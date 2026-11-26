import { Fragment, useState, useEffect, useRef } from 'react';
import Select from 'react-select';

// URL da API, injetada em build time via VITE_API_URL (no GitHub Pages aponta
// para o Cloud Run). O fallback mantem o desenvolvimento local funcionando.
const API = (import.meta.env.VITE_API_URL || 'http://127.0.0.1:33857').replace(/\/$/, '');

// A entrada padrao da aba Gene List e um gene por linha, opcionalmente seguido
// do log2FC ("SOX10  2.41"). Com log2FC a API troca o ssGSEA por um score
// ponderado; sem ele o comportamento continua o de sempre.
const GENES_EXEMPLO = ['SOX10', 'TSPAN6', 'SCYL3', 'GABRB3', 'GAD1', 'GAD2', 'SLC32A1'];

const ASSINATURA_EXEMPLO = [
  ['SOX10', 2.41], ['TSPAN6', -1.87], ['SCYL3', 0.92], ['GABRB3', -2.15],
  ['GAD1', 3.02], ['GAD2', 2.58], ['SLC32A1', -0.55],
];

// Espelha parse_assinatura() do plumber.R, mas so para exibicao (rotulo do
// historico); quem pontua de fato e a API. Se as duas regras divergirem, a
// tela mente sobre o que foi enviado.
const parseAssinatura = (texto) => {
  const divide = (linha) => linha.split(/[ \t,;]+/).filter(Boolean);
  const linhas = (texto || '').split(/\r?\n/).map((l) => l.trim()).filter(Boolean);

  // Cabecalho de export do DESeq2/edgeR/limma ("gene,log2FoldChange")
  const primeiro = linhas.length ? divide(linhas[0]) : [];
  if (primeiro.length >= 2 && /^(gene|genes|symbol|id)$/i.test(primeiro[0])
      && !Number.isFinite(Number(primeiro[1]))) {
    linhas.shift();
  }

  return linhas.flatMap((linha) => {
    const tokens = divide(linha);
    if (tokens.length === 2 && Number.isFinite(Number(tokens[1]))) {
      return [{ gene: tokens[0], peso: Number(tokens[1]), tinhaPeso: true }];
    }
    // Um simbolo so, ou varios sem numero nenhum: lista simples.
    return tokens.map((gene) => ({ gene, peso: 1, tinhaPeso: false }));
  });
};

const separarGenes = (texto) => parseAssinatura(texto).map((g) => g.gene);

// Per-area significance tests offered by the API (significance= on the Gene
// List endpoints). The text under the control is the only place the screen
// says what a mark means before the map arrives; the map caption repeats it.
const TESTES_SIGNIFICANCIA = [
  { value: 'none', label: 'Off', help: 'No significance test: the map shows the score only.' },
  { value: 'score', label: 'Score (*)',
    help: "Asterisk: the list's score in this area beats 4,000 random gene sets of the same size and expression level (q < 0.05), and still does when any single donor is left out." },
  { value: 'markers', label: 'Markers (◆)',
    help: "Diamond: genes specific to this area are over-represented in the list (Fisher's exact test, q < 0.05, at least 2 shared genes). With log2FC values, triangles mark up- (▲) and down-regulated (▼) genes separately." },
];

const formataQ = (q) => (q === null || q === undefined ? '—' : q < 0.001 ? q.toExponential(1) : q.toFixed(3));

// Stage test (/stage_enrichment): one row per (list, region, window). z is
// centred over the 5 windows, so blue = lower than random genes there, red =
// higher. The panel only draws it; the numbers go to Download Data (CSV).
const JANELAS_ESTAGIO = ['1st trimester', '2nd trimester', '3rd trimester', 'Child', 'Adult'];
const JANELAS_CURTAS = { '1st trimester': '1st', '2nd trimester': '2nd', '3rd trimester': '3rd', Child: 'Child', Adult: 'Adult' };
const ROTULOS_ESTAGIO = { all: 'All genes', up: '▲ Up-regulated', down: '▼ Down-regulated' };
const ROTULOS_CURTOS = { all: 'All', up: '▲ Up', down: '▼ Down' };
// Same navy-to-maroon palette as the paper figures, capped at |z| = 10.
const ESCALA_Z = [[-10, [26, 49, 139]], [-5, [79, 113, 190]], [0, [242, 242, 242]], [5, [209, 73, 140]], [10, [122, 8, 69]]];
const corZ = (z, alfa = 1) => {
  const v = Math.max(-10, Math.min(10, Number(z) || 0));
  for (let i = 1; i < ESCALA_Z.length; i++) {
    const [z0, c0] = ESCALA_Z[i - 1];
    const [z1, c1] = ESCALA_Z[i];
    if (v <= z1) {
      const t = (v - z0) / (z1 - z0);
      return `rgba(${c0.map((c, j) => Math.round(c + t * (c1[j] - c))).join(',')},${alfa})`;
    }
  }
  return `rgba(122,8,69,${alfa})`;
};
// Whole brain first, then PsychENCODE's six regions (the API's `region` field).
const REGIOES_ESTAGIO = [
  { id: 'brain', rotulo: 'Brain', nome: 'Whole brain: the six regions below, equally weighted' },
  { id: 'NCX', rotulo: 'NCX', nome: 'Neocortex (11 areas, equally weighted)' },
  { id: 'HIP', rotulo: 'HIP', nome: 'Hippocampus' },
  { id: 'AMY', rotulo: 'AMY', nome: 'Amygdala' },
  { id: 'STR', rotulo: 'STR', nome: 'Striatum' },
  { id: 'MD', rotulo: 'MD', nome: 'Mediodorsal nucleus of the thalamus (dorsal thalamus in the 1st trimester)' },
  { id: 'CBC', rotulo: 'CBC', nome: 'Cerebellar cortex (whole cerebellum in early samples)' },
];
// Short line per list: the windows where the whole brain peaks.
const resumoEstagio = (linhas) => {
  const picos = linhas.filter((l) => l.region === 'brain' && l.enriched && l.direction === 'higher');
  return picos.length ? `peaks in ${picos.map((l) => JANELAS_CURTAS[l.window]).join(', ')}` : 'no whole-brain peak';
};
// Cells with too few donors keep a pale tint of their colour under a grey hatch:
// shown, not tested.
const HACHURA = 'repeating-linear-gradient(135deg, rgba(113,113,122,0.22) 0 2px, transparent 2px 6px)';

export default function Home() {
  const [abaAtual, setAbaAtual] = useState('genelist');
  const [escalaRegiao, setEscalaRegiao] = useState('micro');
  const [opcoes, setOpcoes] = useState([]);
  const [selecionada, setSelecionada] = useState(null);
  const [carregandoLista, setCarregandoLista] = useState(false);
  const [carregandoImagem, setCarregandoImagem] = useState(false);
  const [inputValue, setInputValue] = useState('');

  const [textoListaGenes, setTextoListaGenes] = useState(''); 
  const [imagemCustomizada, setImagemCustomizada] = useState(null);

  const [historicoGene, setHistoricoGene] = useState([]);
  const [historicoOntology, setHistoricoOntology] = useState([]);
  const [historicoGenelist, setHistoricoGenelist] = useState([]);
  const [genesInvalidos, setGenesInvalidos] = useState([]);
  // Errors show as a dismissible line above the map instead of an alert() popup.
  const [aviso, setAviso] = useState(null);
  const [infoAssinatura, setInfoAssinatura] = useState(null);
  const [mostrarModalGenes, setMostrarModalGenes] = useState(false);

  // Which test marks the map, and the per-cell results of both tests (table).
  const [testeSignificancia, setTesteSignificancia] = useState('score');
  // Which developmental window the list peaks in (stage test).
  const [resultadoEstagio, setResultadoEstagio] = useState(null);
  const [carregandoEstagio, setCarregandoEstagio] = useState(false);

  const [genesDaVia, setGenesDaVia] = useState([]);
  const [mostrarGenesDaVia, setMostrarGenesDaVia] = useState(false);
  const [carregandoGenesVia, setCarregandoGenesVia] = useState(false);

  const [baixandoPlot, setBaixandoPlot] = useState(null);

  const [modalAberto, setModalAberto] = useState(false);
  const [abaMapaRef, setAbaMapaRef] = useState('lateral');
  const [escalaMapaRef, setEscalaMapaRef] = useState('micro');

  const inputArquivoRef = useRef(null);

  // Le o arquivo (.txt/.csv/.tsv) como texto puro e joga no mesmo estado do
  // textarea. O parseAssinatura ja aceita virgula/tab/espaco e cabecalho de
  // DESeq2/edgeR/limma, entao nao ha parsing extra para fazer aqui.
  const handleUploadArquivo = (evento) => {
    const arquivo = evento.target.files?.[0];
    if (!arquivo) return;

    const leitor = new FileReader();
    leitor.onload = (e) => {
      setTextoListaGenes(e.target.result);
    };
    leitor.onerror = () => {
      setAviso('Could not read the file.');
    };
    leitor.readAsText(arquivo);

    // Permite selecionar o mesmo arquivo de novo depois (o <input> nao dispara
    // onChange se o value nao mudar).
    evento.target.value = '';
  };

  // Buscando lista de Genes ou Ontologias no R (com Retry automático)
  useEffect(() => {
    setSelecionada(null);
    if (abaAtual === 'genelist') return;

    let montado = true;
    setCarregandoLista(true);
    setOpcoes([]);
    setImagemCustomizada(null);

    setSelecionada(null);

    const endpoint = abaAtual === 'gene' ? '/list_genes' : '/list_ontologies';

    const buscarComRetry = () => {
      fetch(`${API}${endpoint}`)
        .then((res) => {
          if (!res.ok) throw new Error("API starting up...");
          return res.json();
        })
        .then((dados) => {
          if (!montado) return;
          const arrayDados = Array.isArray(dados) ? dados : dados[0] || [];
          const formatadas = arrayDados.map(item => ({ value: item, label: item }));
          
          setOpcoes(formatadas);
          setCarregandoLista(false);
          
          // Se não houver nada selecionado, seleciona o primeiro
          if (formatadas.length > 0) {
            setSelecionada(formatadas[0]);
            atualizarHistorico(formatadas[0].value, abaAtual);
            setCarregandoImagem(true);
          }
        })
        .catch(() => { if (montado) setTimeout(buscarComRetry, 2000); });
    };

    buscarComRetry();
    return () => { montado = false; };
  }, [abaAtual]);

  // Atualiza os históricos limitando a 3 itens
  const atualizarHistorico = (termo, aba) => {
    if (!termo || termo.trim() === '') return;
    
    if (aba === 'gene') {
      setHistoricoGene(prev => [termo, ...prev.filter(t => t !== termo)].slice(0, 3));
    } else if (aba === 'ontology') {
      setHistoricoOntology(prev => [termo, ...prev.filter(t => t !== termo)].slice(0, 3));
    } else if (aba === 'genelist') {
      setHistoricoGenelist(prev => [termo, ...prev.filter(t => t !== termo)].slice(0, 3));
    }
  };

  // FUNÇÃO: Processar Genelist Nova
  const gerarMapaLista = async () => {
    if (!textoListaGenes.trim()) return; 
    
    setCarregandoImagem(true);
    setSelecionada({ label: "Custom List" }); 
    atualizarHistorico(textoListaGenes, 'genelist');

    // 1. Limpa o aviso da pesquisa anterior para não acumular
    setGenesInvalidos([]);
    setAviso(null);
    setInfoAssinatura(null);
    setMostrarModalGenes(false);
    setResultadoEstagio(null);

    try {
      const params = new URLSearchParams();
      params.append('gene_string', textoListaGenes);
      params.append('scale', escalaRegiao);
      params.append('significance', testeSignificancia);

      // 2. CHECAGEM RÁPIDA: Quem existe e quem não existe?
      const resValidacao = await fetch(`${API}/validate_genelist`, {
        method: 'POST', body: params
      });
      const dadosValidacao = await resValidacao.json();

      // Formata a resposta do R com segurança
      const listaInvalidos = Array.isArray(dadosValidacao.invalid) ? dadosValidacao.invalid : dadosValidacao.invalid?.[0] || [];
      const listaValidos = Array.isArray(dadosValidacao.valid) ? dadosValidacao.valid : dadosValidacao.valid?.[0] || [];

      // Escalares vem desembrulhados da API atual, mas uma versao mais antiga
      // no Cloud Run ainda manda [false] / [3].
      const desembrulha = (v) => (Array.isArray(v) ? v[0] : v);
      const ponderado = Boolean(desembrulha(dadosValidacao.weighted));
      const nUp = Number(desembrulha(dadosValidacao.n_up)) || 0;
      const nDown = Number(desembrulha(dadosValidacao.n_down)) || 0;
      const nSemPeso = Number(desembrulha(dadosValidacao.n_unweighted)) || 0;
      setInfoAssinatura(ponderado ? { nUp, nDown, nSemPeso } : null);

      // Se houver algum erro, salva no Estado (Isso faz a caixa amarela aparecer!)
      if (listaInvalidos.length > 0) {
        setGenesInvalidos(listaInvalidos); 
      }

      // Sem nenhum log2FC diferente de zero nao ha o que ponderar, e o erro da
      // API voltaria como um 500 generico.
      if (ponderado && nUp === 0 && nDown === 0) {
        setAviso('All submitted log2FC values are zero — there is nothing to weight.');
        setCarregandoImagem(false);
        return;
      }

      // Se TUDO estiver errado, para por aqui e avisa
      if (listaValidos.length === 0) {
        setAviso('None of the submitted genes was found in the expression data.');
        setCarregandoImagem(false);
        return; 
      }

      // 3. GERAÇÃO DA IMAGEM
      const resposta = await fetch(`${API}/plot_genelist_png`, {
        method: 'POST', body: params
      });
      if (!resposta.ok) throw new Error("Request failed");

      const blob = await resposta.blob();
      setImagemCustomizada(URL.createObjectURL(blob));

      // The stage panel is fetched after the image so it never delays the map.
      if (testeSignificancia !== 'none') buscarEstagio();
    } catch (erro) {
      console.error(erro);
      setAviso('Error processing the list.');
    } finally {
      setCarregandoImagem(false);
    }
  };

  const buscarEstagio = async () => {
    setCarregandoEstagio(true);
    try {
      const params = new URLSearchParams();
      params.append('gene_string', textoListaGenes);
      const res = await fetch(`${API}/stage_enrichment`, { method: 'POST', body: params });
      if (!res.ok) throw new Error('Request failed');
      setResultadoEstagio(await res.json());
    } catch (erro) {
      console.error(erro);
      setResultadoEstagio({ erro: true });
    } finally {
      setCarregandoEstagio(false);
    }
  };

  // Historico
  const carregarDoHistorico = (termo) => {
    if (abaAtual !== 'genelist') {
      setSelecionada({ value: termo, label: termo });
      atualizarHistorico(termo, abaAtual);
      setCarregandoImagem(true);
    } else {
      setTextoListaGenes(termo);
    }
  };

  // Os endpoints _png devolvem o mesmo mapa em bitmap. A tela usa o PNG (o SVG
  // do ggseg traz milhares de poligonos e pesa no navegador); o SVG continua
  // inteiro no botao de download, que e o que vai para a publicacao.
  const endpointPlot = (formato) => {
    const sufixo = formato === 'png' ? '_png' : '';
    if (abaAtual === 'gene') {
      return `${API}/plot_brain${sufixo}?gene=${encodeURIComponent(selecionada.value)}&scale=${escalaRegiao}`;
    }
    if (abaAtual === 'ontology') {
      return `${API}/plot_ontology${sufixo}?geneset=${encodeURIComponent(selecionada.value)}&scale=${escalaRegiao}`;
    }
    return `${API}/plot_genelist${sufixo}`;
  };

  let urlImagem = '';
  if (abaAtual === 'genelist' && imagemCustomizada) {
    urlImagem = imagemCustomizada;
  } else if (selecionada && abaAtual !== 'genelist') {
    urlImagem = endpointPlot('png');
  }


  const baixarPlot = async (formato) => {
    if (!urlImagem || baixandoPlot) return;
    setBaixandoPlot(formato);
    // O PNG da aba Gene List ja esta na tela: refazer o POST custaria outro
    // ssGSEA inteiro so para produzir o mesmo arquivo.
    const jaBaixado = formato === 'png' && abaAtual === 'genelist' && imagemCustomizada;
    try {
      let url = imagemCustomizada;

      if (!jaBaixado) {
        const opcoesFetch = {};
        if (abaAtual === 'genelist') {
          const params = new URLSearchParams();
          params.append('gene_string', textoListaGenes);
          params.append('scale', escalaRegiao);
          params.append('significance', testeSignificancia);
          opcoesFetch.method = 'POST';
          opcoesFetch.body = params;
        }

        const response = await fetch(endpointPlot(formato), opcoesFetch);
        if (!response.ok) throw new Error("Request failed");
        url = URL.createObjectURL(await response.blob());
      }

      const a = document.createElement('a');
      a.href = url;
      a.download = `Map_${abaAtual}.${formato}`;
      a.click();
      // Revogar o blob da aba Gene List apagaria a imagem que esta na tela.
      if (!jaBaixado) URL.revokeObjectURL(url);
    } catch (err) {
      setAviso(`Error downloading the ${formato.toUpperCase()}.`);
    } finally {
      setBaixandoPlot(null);
    }
  };

  const baixarDados = async () => {
    try {
      let url = '';
      let opcoesFetch = {};

    if (abaAtual === 'gene') {
        url = `${API}/data_brain?gene=${selecionada.value}&scale=${escalaRegiao}`;
      } else if (abaAtual === 'ontology') {
        url = `${API}/data_ontology?geneset=${selecionada.value}&scale=${escalaRegiao}`;
      } else if (abaAtual === 'genelist') {
        url = `${API}/data_genelist`;
        const params = new URLSearchParams();
        params.append('gene_string', textoListaGenes);
        params.append('scale', escalaRegiao);
        // With any test on, the CSV carries all of them: both area tests and the
        // stage test, as columns of each (window, region) row.
        params.append('significance', testeSignificancia === 'none' ? 'none' : 'all');
        opcoesFetch = { method: 'POST', body: params };
      }

      const response = await fetch(url, opcoesFetch);
      const blob = await response.blob();
      const objUrl = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = objUrl;
      a.download = `Data_${abaAtual}_${escalaRegiao}.csv`;
      a.click();
      URL.revokeObjectURL(objUrl);
    } catch (err) {
      setAviso('Error downloading the data table.');
    }
  };

  // The stage strip sits beside the map on the Gene List tab once a test ran.
  const mostrarEstagio = abaAtual === 'genelist' && imagemCustomizada && (carregandoEstagio || resultadoEstagio);

  // Histórico ativo atual
  const historicoAtivo = abaAtual === 'gene' ? historicoGene : abaAtual === 'ontology' ? historicoOntology : historicoGenelist;

  const buscarGenesDaVia = async (nomeOntologia) => {
    if (!nomeOntologia) return;
    
    setCarregandoGenesVia(true);
    setMostrarGenesDaVia(true);
    setGenesDaVia([]); 
    
    try {
      // encodeURIComponent protege nomes complexos na hora de virar URL
      const urlFetch = `${API}/geneset_genes?geneset=${encodeURIComponent(nomeOntologia)}`;
      
      const res = await fetch(urlFetch);
      const dados = await res.json();

      let listaLimpa = [];
      if (dados.genes && Array.isArray(dados.genes)) {
        listaLimpa = dados.genes;
      } else if (Array.isArray(dados)) {
        listaLimpa = dados; 
      } else if (dados.genes && typeof dados.genes === 'string') {
        listaLimpa = [dados.genes];
      }

      setGenesDaVia(listaLimpa);

    } catch (err) {
      console.error("Error fetching gene set members:", err);
      setGenesDaVia(["Error loading the genes."]);
    } finally {
      setCarregandoGenesVia(false);
    }
  };

  // Three columns on wide screens -- controls | map | stage. On lg+ the page is
  // exactly the height of <main>: every column is capped at it, the map image
  // scales to fit, and a column that still cannot fit scrolls on its own.
  return (
    <section className="max-w-[1800px] mx-auto px-6 py-4 lg:h-full">
      <div className="flex flex-col lg:flex-row gap-6 items-start lg:h-full">
        
        {/* DIVISÃO ENTRE ABAS MICRO - MACRO REGIÃO */}
        <div className="w-full lg:w-[330px] shrink-0 bg-[#f2f4f4] rounded-xl p-6 border border-zinc-100 h-fit lg:max-h-full lg:overflow-y-auto">

          <div className="flex bg-zinc-300 p-1 rounded-lg mb-4">
            {['micro', 'macro'].map((escala) => (
              <button 
                key={escala}
                onClick={() => {
                  setEscalaRegiao(escala);
                  setCarregandoImagem(true);
                }}
                className={`flex-1 py-2 text-[11px] font-black uppercase tracking-widest rounded-md transition-all ${escalaRegiao === escala ? 'bg-[#58614c] text-white shadow-md' : 'text-zinc-500 hover:text-zinc-700'}`}
              >
                {escala === 'micro' ? 'MicroRegions' : 'MacroRegions'}
              </button>
            ))}
          </div>

          <div className="flex bg-[#e4e9ea] p-1 rounded-lg mb-6">
            {['genelist','gene','ontology'].map((tab) => (
              <button 
                key={tab}
                onClick={() => setAbaAtual(tab)}
                className={`flex-1 py-2 text-[10px] font-bold uppercase tracking-widest rounded-md transition-all ${abaAtual === tab ? 'bg-white text-[#58614c]' : 'text-[#596061]'}`}
              >
                {tab === 'genelist' ? 'Gene List' : tab}
              </button>
            ))}
          </div>

          <div className="space-y-4">
            
            {/* CABEÇALHO COM O BOTÃO DE EXEMPLO (Lado a lado) */}
            <div className="flex justify-between items-center">
              <label className="block text-[10px] font-bold uppercase tracking-widest text-[#58614c]">
                {abaAtual === 'gene' ? 'Search Identifier' : abaAtual === 'ontology' ? 'Search Ontology Term' : 'Input List'}
              </label>

              {/* SÓ MOSTRA O BOTÃO SE ESTIVER NA ABA GENELIST */}
              {abaAtual === 'genelist' && (
                <div className="flex gap-3">
                  <button
                    onClick={() => setTextoListaGenes(GENES_EXEMPLO.join("\n"))}
                    className="text-[10px] text-zinc-500 hover:text-[#58614c] font-bold underline transition-colors"
                  >
                    Example list
                  </button>
                  <button
                    onClick={() => setTextoListaGenes(
                      ASSINATURA_EXEMPLO.map(([gene, fc]) => `${gene}\t${fc}`).join("\n"))}
                    className="text-[10px] text-zinc-500 hover:text-[#58614c] font-bold underline transition-colors"
                  >
                    Example with log2FC
                  </button>
                  <button
                    onClick={() => inputArquivoRef.current?.click()}
                    className="text-[10px] text-zinc-500 hover:text-[#58614c] font-bold underline transition-colors"
                  >
                    Upload file
                  </button>
                  <input
                    ref={inputArquivoRef}
                    type="file"
                    accept=".txt,.csv,.tsv"
                    onChange={handleUploadArquivo}
                    className="hidden"
                  />
                </div>
              )}
            </div>

            {abaAtual !== 'genelist' ? (
              <div className="flex flex-col gap-3">
                <Select
                  value={selecionada}
                  onChange={(opt) => { 
                    setSelecionada(opt); 
                    atualizarHistorico(opt.value, abaAtual);
                    setCarregandoImagem(true); 
                    setMostrarGenesDaVia(false); // Esconde a lista ao trocar de via
                  }}
                  onInputChange={(val) => setInputValue(val)}
                  options={opcoes.filter(o => o.label.toLowerCase().includes(inputValue.toLowerCase())).slice(0, 100)}
                  filterOption={null}
                  styles={{ menu: (base) => ({ ...base, zIndex: 9999 }) }}
                  isLoading={carregandoLista}
                  placeholder="Type to search..."
                  className="text-sm"
                />
                
                {/* BOTÃO E CAIXA DE GENES DA VIA (Só aparece na aba Ontology e se algo estiver selecionado) */}
                {abaAtual === 'ontology' && selecionada && (
                  <div className="flex flex-col gap-2">
                    <button
                      onClick={() => mostrarGenesDaVia ? setMostrarGenesDaVia(false) : buscarGenesDaVia(selecionada.value)}
                      className="text-[10px] font-bold uppercase tracking-widest text-[#58614c] bg-[#e4e9ea] hover:bg-zinc-300 py-2 rounded-lg transition-colors w-full"
                    >
                      {mostrarGenesDaVia ? 'Hide Genes in this Set' : '👁️ Show Genes in this Set'}
                    </button>

                    {mostrarGenesDaVia && (
                      <div className="bg-white p-3 rounded-lg border border-zinc-200 shadow-inner animate-fade-in">
                        {carregandoGenesVia ? (
                          <p className="text-xs text-zinc-400 text-center py-4">Loading genes...</p>
                        ) : (
                          <div className="flex flex-col gap-2">
                            <div className="flex justify-between items-center">
                              <span className="text-[10px] font-bold text-zinc-500">{genesDaVia.length} genes present in the data:</span>
                              <button
                                onClick={() => {
                                  navigator.clipboard.writeText(genesDaVia.join('\n') + '\n');
                                }}
                                className="text-[10px] bg-zinc-100 hover:bg-zinc-200 text-zinc-700 px-2 py-1 rounded"
                              >
                                Copy
                              </button>
                            </div>
                            <textarea
                              readOnly
                              value={genesDaVia.join('\n')}
                              className="w-full bg-zinc-50 text-zinc-600 font-mono text-[10px] p-2 rounded border border-zinc-200 outline-none resize-y h-24"
                            />
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                )}
              </div>
            ) : (
              <div>
                <textarea 
                  rows="6" 
                  value={textoListaGenes}
                  onChange={(e) => setTextoListaGenes(e.target.value)}
                  placeholder={"SOX10\nTSPAN6\nGAD1"} 
                  className="w-full px-4 py-3 bg-white border border-transparent focus:border-[#58614c] rounded-lg text-sm font-mono outline-none resize-y"
                ></textarea>

                {/* Applies to plain lists (ssGSEA) and log2FC signatures alike. */}
                <div className="mt-3">
                  <span className="block text-[10px] font-bold uppercase tracking-widest text-[#58614c] mb-1.5">
                    Mark significant areas
                  </span>
                  <div className="flex bg-[#e4e9ea] p-1 rounded-lg">
                    {TESTES_SIGNIFICANCIA.map((t) => (
                      <button
                        key={t.value}
                        onClick={() => setTesteSignificancia(t.value)}
                        className={`flex-1 py-1.5 text-[10px] font-bold uppercase tracking-widest rounded-md transition-all ${testeSignificancia === t.value ? 'bg-white text-[#58614c]' : 'text-[#596061]'}`}
                      >
                        {t.label}
                      </button>
                    ))}
                  </div>
                  <p className="mt-1.5 text-[11px] leading-snug text-zinc-600">
                    {TESTES_SIGNIFICANCIA.find((t) => t.value === testeSignificancia).help}
                    {testeSignificancia !== 'none' && ' Needs at least 5 genes found in the data.'}
                  </p>
                </div>

                <button 
                  onClick={gerarMapaLista}
                  className="w-full mt-4 bg-[#58614c] text-white py-3 rounded-lg font-bold text-xs uppercase tracking-widest hover:opacity-90 active:scale-[0.98] transition-all"
                >
                  Generate Plot
                </button>

                {/* Unica confirmacao na tela de que os log2FC foram lidos: sem
                    isso um erro de formato vira um mapa nao ponderado silencioso. */}
                {infoAssinatura && (
                  <div className="mt-3 flex flex-wrap items-center gap-x-2 gap-y-1 text-[11px] bg-[#e9edf5] border border-[#c9d4e8] text-[#2f3f66] rounded-lg px-3 py-2">
                    <span className="text-[10px] font-bold uppercase tracking-widest">Weighted</span>
                    <span>{infoAssinatura.nUp} up / {infoAssinatura.nDown} down</span>
                    {infoAssinatura.nSemPeso > 0 && (
                      <span className="text-[#8c6d1f]">
                        · {infoAssinatura.nSemPeso} without log2FC (counted as +1)
                      </span>
                    )}
                  </div>
                )}
              </div>
            )}
          </div>

          {/* HISTÓRICO RECENTE */}
          {historicoAtivo.length > 0 && (
              <div className="mt-8 pt-6 border-t border-zinc-200">
                <span className="text-[10px] font-bold uppercase tracking-widest text-zinc-400 block mb-3">Recent:</span>
                <div className="flex flex-wrap gap-2">
                  {historicoAtivo.map((termo, idx) => (
                    <button 
                      key={idx}
                      onClick={() => carregarDoHistorico(termo)}
                      className="bg-white hover:bg-zinc-100 border border-zinc-200 text-zinc-600 text-[11px] font-medium px-3 py-1.5 rounded-full transition-colors max-w-full truncate"
                      title={termo}
                    >
                      {(() => {
                        const rotulo = separarGenes(termo).join(', ');
                        return rotulo.length > 25 ? rotulo.substring(0, 25) + '...' : rotulo;
                      })()}
                    </button>
                  ))}
                </div>
              </div>
            )}

        </div>

        {/* COLUNA DIREITA: mapa + botoes, e o bloco de estagio ao lado (xl) ou abaixo */}
        <div className="flex-1 min-w-0 w-full flex flex-col xl:flex-row gap-6 items-start lg:max-h-full lg:overflow-y-auto xl:h-full xl:overflow-visible">
        <div className="flex-1 min-w-0 w-full flex flex-col gap-3">
          {aviso && (
            <div className="flex items-center justify-between gap-3 bg-[#fff8e6] border border-[#f0dca5] rounded-xl px-4 py-2.5 text-sm text-[#8c6d1f]">
              <span>⚠️ {aviso}</span>
              <button onClick={() => setAviso(null)} className="text-[#8c6d1f] hover:opacity-70 text-base leading-none" aria-label="Dismiss">×</button>
            </div>
          )}
          
        {genesInvalidos.length > 0 && (
          <div className="bg-[#fff8e6] border border-[#f0dca5] rounded-xl p-5 flex flex-col gap-3 shadow-sm">
            <div className="flex justify-between items-center">
              <p className="text-sm text-[#8c6d1f] font-medium">
                ⚠️ <strong>{genesInvalidos.length} gene(s)</strong> not found in the data and ignored.
              </p>
              <div className="flex gap-2">
                <button
                  onClick={() => {
                    navigator.clipboard.writeText(genesInvalidos.join('\n') + '\n');
                  }}
                  className="text-[11px] font-bold uppercase tracking-widest bg-white border border-[#f0dca5] text-[#8c6d1f] px-4 py-2 rounded-lg hover:bg-[#fcf4dc] transition-colors"
                >
                  Copy Genes
                </button>
          
          <button
            onClick={() => setMostrarModalGenes(!mostrarModalGenes)}
            className="text-[11px] font-bold uppercase tracking-widest bg-white border border-[#f0dca5] text-[#8c6d1f] px-4 py-2 rounded-lg hover:bg-[#fcf4dc] transition-colors"
          >
            {mostrarModalGenes ? 'Hide List' : 'Show List'}
          </button>
        </div>
      </div>

      {mostrarModalGenes && (
        <textarea
          readOnly
          value={genesInvalidos.join('\n')}
          className="w-full mt-2 bg-white text-[#8c6d1f] font-mono text-xs p-3 rounded-lg border border-[#f0dca5] outline-none resize-y h-32 cursor-text focus:border-[#8c6d1f]"
        />
      )}
    </div>
  )}

          {/* Caixa do Mapa: the image is capped at the window height (minus the
              header, page padding and download buttons) so it never scrolls. */}
          <div className="bg-white border border-zinc-100 rounded-xl overflow-hidden shadow-sm p-4 flex items-center justify-center min-h-[400px] relative">
            {carregandoImagem && <p className="absolute font-bold text-[#58614c] animate-pulse z-10 bg-white/80 px-4 py-2 rounded">⏳ Rendering map...</p>}
            {urlImagem && (
              <img 
                src={urlImagem} 
                className="max-w-full w-auto h-auto max-h-[calc(100vh-230px)] transition-opacity duration-300" 
                style={{ opacity: carregandoImagem ? 0.2 : 1 }}
                onLoad={() => setCarregandoImagem(false)}
                alt="Brain Map"
              />
            )}
          </div>

          {/* Botões de Download */}
          {urlImagem && (
            <div className="flex flex-wrap gap-3 justify-end">
              <button 
                onClick={baixarDados}
                className="flex items-center gap-2 px-4 py-2.5 whitespace-nowrap bg-white border border-zinc-200 text-zinc-700 rounded-lg text-xs font-bold uppercase tracking-widest hover:border-zinc-300 hover:bg-zinc-50 transition-all shadow-sm"
              >
                📊 Download Data (CSV)
              </button>
              <button 
                onClick={() => baixarPlot('svg')}
                disabled={baixandoPlot !== null}
                className="flex items-center gap-2 px-4 py-2.5 whitespace-nowrap bg-[#58614c] text-white rounded-lg text-xs font-bold uppercase tracking-widest hover:opacity-90 shadow-sm transition-all disabled:opacity-60 disabled:cursor-wait"
              >
                {baixandoPlot === 'svg' ? '⏳ Generating SVG...' : '🖼️ Download Plot (SVG)'}
              </button>
              <button 
                onClick={() => baixarPlot('png')}
                disabled={baixandoPlot !== null}
                className="flex items-center gap-2 px-4 py-2.5 whitespace-nowrap bg-[#58614c] text-white rounded-lg text-xs font-bold uppercase tracking-widest hover:opacity-90 shadow-sm transition-all disabled:opacity-60 disabled:cursor-wait"
              >
                {baixandoPlot === 'png' ? '⏳ Generating PNG...' : '🖼️ Download Plot (PNG)'}
              </button>
            </div>
          )}

        </div>

        {/* Stage test: its own card, right of the map on wide screens */}
        {mostrarEstagio && (
          <div className="w-full xl:w-[240px] shrink-0 bg-white border border-zinc-100 rounded-xl shadow-sm p-4 flex flex-col gap-3 xl:max-h-full xl:overflow-y-auto">
              <span className="text-[10px] font-bold uppercase tracking-widest text-[#58614c]">Developmental stage</span>
              {carregandoEstagio && <p className="text-[11px] text-[#58614c] animate-pulse">⏳ Testing stages...</p>}
              {resultadoEstagio?.erro && <p className="text-[11px] text-[#8c6d1f]">Could not compute the stage test.</p>}
              {resultadoEstagio && !resultadoEstagio.erro && resultadoEstagio.reason && (
                <p className="text-[11px] text-zinc-500">Not tested — {resultadoEstagio.reason}.</p>
              )}
              {resultadoEstagio && !resultadoEstagio.erro && !resultadoEstagio.reason && (() => {
                const linhas = resultadoEstagio.rows || [];
                const listas = ['all', 'up', 'down'].filter((l) => linhas.some((r) => r.list === l));
                const celula = (lista, regiao, janela) =>
                  linhas.find((r) => r.list === lista && r.region === regiao && r.window === janela);
                const doadores = (janela) => celula('all', 'brain', janela)?.n_donors;
                const skipped = resultadoEstagio.skipped || [];
                return (
                  <>
                    <table className="text-[10px] self-start" style={{ borderCollapse: 'separate', borderSpacing: '2px' }}>
                      <thead>
                        <tr>
                          <th />
                          {JANELAS_ESTAGIO.map((w) => (
                            <th key={w} className="font-bold text-zinc-500 text-center w-[28px]" title={`${w}, ${doadores(w)} donors`}>
                              {JANELAS_CURTAS[w]}
                            </th>
                          ))}
                        </tr>
                      </thead>
                      <tbody>
                        {listas.map((l) => (
                          <Fragment key={l}>
                            {listas.length > 1 && (
                              <tr>
                                <td colSpan={JANELAS_ESTAGIO.length + 1} className="pt-1 font-bold text-zinc-700 whitespace-nowrap">
                                  {ROTULOS_ESTAGIO[l]}
                                </td>
                              </tr>
                            )}
                            {REGIOES_ESTAGIO.map((r) => (
                              <tr key={r.id}>
                                <td className={`pr-1.5 py-0 leading-[14px] whitespace-nowrap ${r.id === 'brain' ? 'font-bold text-zinc-700' : 'text-zinc-500'}`} title={r.nome}>
                                  {r.rotulo}
                                </td>
                                {JANELAS_ESTAGIO.map((w) => {
                                  const c = celula(l, r.id, w);
                                  if (!c) return <td key={w} />;
                                  const marca = !c.tested ? '' : c.enriched ? '★' : c.q < 0.05 ? '○' : '';
                                  const estilo = c.tested
                                    ? { background: corZ(c.z), color: Math.abs(c.z) > 4 ? 'white' : '#27272a' }
                                    : { backgroundColor: corZ(c.z, 0.3), backgroundImage: HACHURA };
                                  const dica = c.tested
                                    ? `${r.nome}, ${w}: z = ${c.z}; q = ${formataQ(c.q)}; ${c.n_donors} donors; worst p leaving one donor out: ${c.p_loo ?? '—'}`
                                    : `${r.nome}, ${w}: not tested, ${c.n_donors} donors (needs 3); z = ${c.z}`;
                                  return (
                                    <td key={w} className="h-[16px] p-0 rounded-sm text-center leading-none text-[9px]" style={estilo} title={dica}>
                                      {marca}
                                    </td>
                                  );
                                })}
                              </tr>
                            ))}
                          </Fragment>
                        ))}
                      </tbody>
                    </table>
                    <div className="flex flex-col gap-0.5 text-[10px] text-zinc-600 w-full max-w-[240px]">
                      {listas.map((l) => (
                        <span key={l}><strong>{listas.length > 1 ? ROTULOS_CURTOS[l] : 'Whole brain'}:</strong> {resumoEstagio(linhas.filter((r) => r.list === l))}</span>
                      ))}
                      {skipped.length > 0 && (
                        <span>{skipped.map((s) => (s === 'up' ? 'Up' : 'Down')).join(' and ')}-regulated genes not tested: fewer than 5.</span>
                      )}
                    </div>
                    <div className="flex flex-col gap-1 text-[10px] text-zinc-500 w-full max-w-[240px]">
                      <div className="h-2 rounded-sm" style={{ background: `linear-gradient(to right, ${[-10, -5, 0, 5, 10].map((z) => corZ(z)).join(', ')})` }} />
                      <div className="flex justify-between"><span>lower</span><span>z vs random genes</span><span>higher</span></div>
                      <span>★ q &lt; 0.05, robust to any one donor</span>
                      <span>○ q &lt; 0.05, rests on one donor</span>
                      <span className="flex items-center gap-1">
                        <span className="inline-block w-3 h-3 rounded-sm" style={{ backgroundColor: 'rgba(161,161,170,0.25)', backgroundImage: HACHURA }} />
                        fewer than 3 donors, not tested
                      </span>
                      <span className="text-zinc-400 leading-snug">
                        Brain: six regions, equal weight. Hover for z and q; values in the CSV.
                      </span>
                    </div>
                  </>
                );
              })()}
          </div>
        )}

        </div>

      </div>

      {/* BOTÃO FLUTUANTE */}
      <button
        onClick={() => setModalAberto(true)}
        style={{
          position: 'fixed',
          bottom: '30px',
          right: '30px',
          backgroundColor: '#58614c',
          color: 'white',
          border: 'none',
          borderRadius: '50%',
          width: '60px',
          height: '60px',
          fontSize: '28px',
          fontWeight: 'bold',
          cursor: 'pointer',
          boxShadow: '0px 4px 10px rgba(0,0,0,0.3)',
          zIndex: 999,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          transition: 'transform 0.2s'
        }}
        onMouseOver={(e) => e.currentTarget.style.transform = 'scale(1.1)'}
        onMouseOut={(e) => e.currentTarget.style.transform = 'scale(1)'}
        title="Reference Map"
      >
        ?
      </button>

      {/* MODAL DE REFERÊNCIA */}
      {modalAberto && (
        <div style={{
          position: 'fixed',
          top: 0, left: 0, width: '100vw', height: '100vh',
          backgroundColor: 'rgba(0,0,0,0.7)',
          display: 'flex', justifyContent: 'center', alignItems: 'center',
          zIndex: 999999
        }}>
          
          <div style={{
            backgroundColor: 'white',
            padding: '30px',
            borderRadius: '12px',
            position: 'relative',
            width: '90%',
            maxWidth: '900px',
            maxHeight: '90vh',
            overflowY: 'auto',
            display: 'block', 
            boxShadow: '0px 15px 40px rgba(0,0,0,0.6)'
          }}>
            
            {/* Botão Fechar */}
            <button
              onClick={() => setModalAberto(false)}
              style={{
                position: 'absolute', top: '15px', right: '15px', 
                backgroundColor: '#e74c3c', color: 'white', border: 'none', 
                borderRadius: '5px', padding: '8px 15px', cursor: 'pointer', fontWeight: 'bold'
              }}
            >
              Close
            </button>

            {/* Título */}
            <h2 style={{ marginTop: 0, marginBottom: '20px', color: '#58614c', fontWeight: 'bold', fontSize: '20px' }}>
              Reference Map (Regions)
            </h2>

              {/*TROCA DE VISUALIZAÇÃO*/}
            <div style={{ 
              display: 'flex', 
              justifyContent: 'center', 
              gap: '10px', 
              marginBottom: '20px' 
            }}>
              {['micro', 'macro'].map((nivel) => (
                <button
                  key={nivel}
                  onClick={() => setEscalaMapaRef(nivel)}
                  style={{
                    padding: '8px 20px',
                    borderRadius: '20px',
                    border: '2px solid #58614c',
                    cursor: 'pointer',
                    fontSize: '12px',
                    fontWeight: 'bold',
                    textTransform: 'uppercase',
                    backgroundColor: escalaMapaRef === nivel ? '#58614c' : 'white',
                    color: escalaMapaRef === nivel ? 'white' : '#58614c',
                    transition: 'all 0.3s'
                  }}
                >
                  {nivel === 'micro' ? 'Microregions' : 'Macroregions'}
                </button>
              ))}
            </div>
            
            {/* MENU DE ABAS */}
            <div style={{ 
              display: 'block', 
              backgroundColor: '#f3f4f6', 
              padding: '10px', 
              borderRadius: '8px', 
              marginBottom: '20px',
              textAlign: 'center' 
            }}>
              {[
                { id: 'lateral', nome: 'Lateral Cortex' },
                { id: 'medial', nome: 'Medial Cortex' },
                { id: 'sagittal', nome: 'Subcortical (Sagittal)' },
                { id: 'coronal', nome: 'Subcortical (Coronal)' }
              ].map((corte) => (
                <button
                  key={corte.id}
                  onClick={() => setAbaMapaRef(corte.id)}
                  style={{
                    display: 'inline-block', 
                    margin: '5px', 
                    padding: '10px 16px',
                    fontSize: '11px',
                    fontWeight: '900',
                    textTransform: 'uppercase',
                    letterSpacing: '1px',
                    borderRadius: '6px',
                    border: 'none',
                    cursor: 'pointer',
                    backgroundColor: abaMapaRef === corte.id ? '#58614c' : '#e4e9ea',
                    color: abaMapaRef === corte.id ? 'white' : '#596061',
                    transition: 'all 0.2s'
                  }}
                >
                  {corte.nome}
                </button>
              ))}
            </div>

            {/* ÁREA DA IMAGEM */}
            <div style={{
              backgroundColor: '#f9fafb',
              border: '1px solid #e5e7eb',
              borderRadius: '8px',
              padding: '15px',
              textAlign: 'center' 
            }}>
              <img
                key={`${abaMapaRef}-${escalaMapaRef}`}
                src={`${API}/plot_reference?view=${abaMapaRef}&scale=${escalaMapaRef}`}
                alt={`Reference map, ${abaMapaRef} view`}
                style={{ 
                  maxWidth: '70%', 
                  height: 'auto', 
                  display: 'inline-block' 
                }}
              />
            </div>
            
          </div>
        </div>
      )}

    </section>
  );
}