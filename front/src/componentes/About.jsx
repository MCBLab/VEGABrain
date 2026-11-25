// Shared text styles, so every section reads the same.
const H2 = ({ children }) => (
  <h2 className="text-sm font-bold text-[#2d3435] uppercase tracking-[0.2em] mb-4">{children}</h2>
);
const P = ({ children, className = '' }) => (
  <p className={`text-base text-[#596061] leading-relaxed ${className}`}>{children}</p>
);
const Strong = ({ children }) => <strong className="font-semibold text-[#2d3435]">{children}</strong>;

// Score (*) vs Markers (◆): one row per question a user is likely to ask.
const COMPARACAO = [
  ['Question it answers',
    'Is my list, as a whole, more active in this area than in the other areas sampled at the same age?',
    'Does my list contain the genes that set this area apart from the other areas at the same age?'],
  ['What it measures',
    "The list's ssGSEA score: every gene counts, through its expression rank in each sample.",
    "The overlap between the list and the area's marker genes. Only area-specific genes count."],
  ['Compared with',
    '4,000 random gene sets of the same size and expression level, scored the same way.',
    "All genes expressed in that window (Fisher's exact test)."],
  ['Areas tested',
    'Each BrainSpan structure (e.g. striatum, hippocampus, V1) in each window.',
    'The six region classes against each other; then each neocortical area against the rest of the neocortex.'],
  ['A mark means',
    'q < 0.05 (Benjamini-Hochberg over all tested cells) and it stays significant when any single donor is left out.',
    'q < 0.05 (Benjamini-Hochberg over all tested cells) and at least 2 genes shared with the list.'],
  ['What you get',
    'Which areas. It does not say which genes drive the result.',
    'Which areas, and which of your genes are behind each one (markers_genes in the CSV).'],
  ['Best for',
    'Pathways, GO terms and gene programs whose genes shift together, including modest, coordinated changes.',
    'Differential-expression lists and signatures: which region’s identity genes your list carries, e.g. a cell-composition or dissection question.'],
  ['Not tested when',
    'Fewer than 3 donors in that structure and window.',
    'Fewer than 10 marker genes or 3 donors. Common after birth, when cortical areas become alike.'],
  ['With log2FC values',
    'Scores the signed signature: each gene weighted by its log2FC.',
    'Tests up-regulated (▲) and down-regulated (▼) genes separately.'],
];

// Columns of Download Data (CSV) on the Gene List tab.
const COLUNAS_CSV = [
  ['broad_age, region', 'Developmental window and atlas parcel of the row.'],
  ['Avg_ssGSEA_Score / Avg_Weighted_Score', 'The value painted on the map (plain list / list with log2FC).'],
  ['tested_score, enriched_score, q_score, p_loo', 'Score test (*) of the structure(s) feeding the parcel. p_loo is the worst p when one donor is left out.'],
  ['enriched_markers, q_markers, markers_level, markers_genes', 'Marker test (◆). Level is class or area; genes are the list genes behind the hit. Suffixes _up / _down for log2FC lists.'],
  ['stage_region, stage_tested, stage_z, stage_q, stage_p_loo, stage_enriched', "Stage test of the parcel's PsychENCODE region in that window."],
  ['stage_brain_*', 'The same for the whole brain (the six regions, equally weighted).'],
  ['_up / _down', 'Copies of a column for the up- and down-regulated genes of a log2FC list.'],
  ['*_note', 'Why a test did not run (e.g. fewer than 5 genes found).'],
];

const REGIOES = [
  ['NCX', 'Neocortex (its 11 areas, equally weighted)'],
  ['HIP', 'Hippocampus'],
  ['AMY', 'Amygdala'],
  ['STR', 'Striatum'],
  ['MD', 'Mediodorsal nucleus of the thalamus (dorsal thalamus in the 1st trimester)'],
  ['CBC', 'Cerebellar cortex (whole cerebellum in early samples)'],
];

// Same palette as the stage card on the Home page.
const GRADIENTE = 'linear-gradient(to right, rgb(26,49,139), rgb(79,113,190), rgb(242,242,242), rgb(209,73,140), rgb(122,8,69))';
const HACHURA = 'repeating-linear-gradient(135deg, rgba(113,113,122,0.22) 0 2px, transparent 2px 6px)';

export default function About() {
  return (
    <section className="max-w-7xl mx-auto px-6 sm:px-8 lg:px-12 py-16">
      <div className="max-w-4xl mx-auto">
        <p className="text-[10px] font-bold text-zinc-400 uppercase tracking-[0.3em] mb-4">
          About This App
        </p>
        <h1 className="text-5xl font-bold text-[#2d3435] leading-tight tracking-tighter mb-8 whitespace-nowrap">
          Human Developmental Brain RNA-Seq Dataset
        </h1>

        <div className="text-xl text-[#596061] leading-relaxed max-w-2xl mb-16 space-y-4">
          <p>
            This application lets you explore gene expression across human
            brain development, using RNA-Seq data from the BrainSpan Atlas.
            It's built for researchers who want to move from a gene, an
            ontology term, or a custom gene list straight to an expression
            pattern, without setting up their own analysis pipeline.
          </p>
          <p>
            The underlying data comes from the{' '}
            <a
              href="https://www.brainspan.org/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
            >
              BrainSpan Atlas of the Developing Human Brain
            </a>
            , which profiles gene expression across brain regions and
            developmental stages, from early fetal periods through adulthood.
          </p>
        </div>

        <div className="max-w-3xl mb-16">
          <H2>What you can do here</H2>
          <div className="grid grid-cols-2 gap-x-12 border-t border-zinc-200">
            <div className="py-6 border-b border-zinc-200">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Search by gene
              </h3>
              <P>
                Look up a single gene and see how its expression changes
                across brain regions and developmental time points.
              </P>
            </div>
            <div className="py-6 border-b border-l border-zinc-200 pl-12">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Score custom gene lists with ssGSEA
              </h3>
              <P>
                Upload your own list of genes and get single-sample GSEA
                enrichment scores computed against the BrainSpan expression
                data, so you can test whether your set of interest shifts
                across regions or stages.
              </P>
            </div>
            <div className="py-6 border-b border-zinc-200">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Browse by ontology
              </h3>
              <P>
                Explore expression patterns for genes grouped under
                biological ontology terms, to see how whole processes or
                pathways trend across development.
              </P>
            </div>
            <div className="py-6 border-b border-l border-zinc-200 pl-12">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Plot custom log2FC values
              </h3>
              <P>
                Bring in your own log2 fold-change values from an external
                comparison — a different disease model, treatment, or
                condition — and plot them against the BrainSpan expression
                data.
              </P>
            </div>
            <div className="py-6">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Test where a list is enriched
              </h3>
              <P>
                Mark the brain areas where your list stands out, with two
                complementary tests: Score (*) and Markers (◆).
              </P>
            </div>
            <div className="py-6 border-l border-zinc-200 pl-12">
              <h3 className="text-lg font-semibold text-[#2d3435] mb-1">
                Test when a list peaks
              </h3>
              <P>
                The Developmental stage card shows in which window, and in
                which of six major regions, your list is higher or lower than
                random genes.
              </P>
            </div>
          </div>
        </div>

        {/* Gene List tab: the significance switch */}
        <div className="max-w-3xl mb-16">
          <H2>Testing a gene list</H2>
          <div className="space-y-4">
            <P>
              On the <Strong>Gene List</Strong> tab, paste one gene per line, or a
              gene followed by its log2FC. The map paints the list's score in every
              area and window. The <Strong>Mark significant areas</Strong> switch
              decides which test, if any, marks the map:
            </P>
            <ul className="space-y-2 text-base text-[#596061] leading-relaxed list-none">
              <li><Strong>Off</Strong>: no test. The map shows the score only, and it is the fastest option.</li>
              <li><Strong>Score (*)</Strong>: an asterisk on each area where the list's score is significantly higher than random genes.</li>
              <li>
                <Strong>Markers (◆)</Strong>: a diamond on each area whose specific genes are
                over-represented in the list. With log2FC values, ▲ marks up-regulated and
                ▼ down-regulated genes.
              </li>
            </ul>
            <P>
              With either test on, the <Strong>Developmental stage</Strong> card appears
              to the right of the map, and <Strong>Download Data (CSV)</Strong> carries
              the results of all three tests. Each test needs at least 5 of your genes
              to be found in the data.
            </P>
          </div>
        </div>

        {/* Score vs Markers */}
        <div className="mb-16">
          <H2>Score (*) or Markers (◆)?</H2>
          <div className="max-w-3xl space-y-4 mb-6">
            <P>
              Both tests ask <em>where</em> your list stands out, and both compare
              areas only with other areas sampled at the same age. They differ in what
              counts as standing out.
            </P>
            <P>
              <Strong>Score (*) is about activity.</Strong> It asks whether the genes of
              your list, taken together, run higher in an area than in its neighbours.
              Every gene contributes, so it picks up programs whose genes move together,
              even when no single gene is specific to that area.
            </P>
            <P>
              <Strong>Markers (◆) is about identity.</Strong> First it finds, for each
              area and window, the genes that distinguish it from the other areas
              (limma, donor as a fixed effect; q &lt; 0.05, at least 2-fold, expressed).
              Then it asks whether your list contains more of them than chance. Only
              area-specific genes count, and it tells you which ones are shared.
            </P>
          </div>
          <div className="overflow-x-auto border-t border-zinc-200">
            <table className="w-full text-sm text-left">
              <thead>
                <tr className="border-b border-zinc-200">
                  <th className="py-3 pr-6 w-[22%]" />
                  <th className="py-3 pr-6 w-[39%] text-[#2d3435] font-semibold text-base">Score (*)</th>
                  <th className="py-3 w-[39%] text-[#2d3435] font-semibold text-base">Markers (◆)</th>
                </tr>
              </thead>
              <tbody>
                {COMPARACAO.map(([rotulo, score, markers]) => (
                  <tr key={rotulo} className="border-b border-zinc-200 align-top">
                    <td className="py-3 pr-6 text-[11px] font-bold uppercase tracking-widest text-zinc-500">{rotulo}</td>
                    <td className="py-3 pr-6 text-[#596061] leading-relaxed">{score}</td>
                    <td className="py-3 text-[#596061] leading-relaxed">{markers}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="max-w-3xl space-y-4 mt-6">
            <P>
              <Strong>They often agree, and the differences are informative.</Strong> A
              set of broadly expressed genes that rise together in the striatum can earn
              an asterisk with no diamond, because none of them is striatum-specific. A
              long list that contains a handful of strong cerebellar markers among many
              unrelated genes can earn a diamond with no asterisk, because those few
              genes barely move the list's overall score.
            </P>
            <P>
              <Strong>Rule of thumb:</Strong> use Score for curated pathways and gene
              programs; use Markers for lists from an experiment, when you want to know
              which region's identity genes are in it and which genes they are.
            </P>
          </div>
        </div>

        {/* Developmental stage */}
        <div className="max-w-3xl mb-16">
          <H2>Developmental stage</H2>
          <div className="space-y-4">
            <P>
              Score and Markers compare areas <em>within</em> a window, so they answer
              <em> where</em> but not <em>when</em>. A gene set that rises in every area at
              once leaves no mark on the map, because no area stands out from the
              others of its own window. Astrocyte development is an example: it is high
              all over the brain in the 3rd trimester and gets no area mark there.
            </P>
            <P>
              The <Strong>Developmental stage</Strong> card answers <em>when</em>. For each
              window it takes the list's score, centres it on the list's average over the
              five windows, and compares it with 4,000 random gene sets of the same
              expression level. The first row is the whole brain; the other rows repeat
              the test inside each of the six regions used by PsychENCODE:
            </P>
            <div className="grid grid-cols-[auto_1fr] gap-x-4 gap-y-1 text-sm text-[#596061]">
              {REGIOES.map(([sigla, nome]) => [
                <span key={`${sigla}-s`} className="font-bold text-[#2d3435]">{sigla}</span>,
                <span key={`${sigla}-n`}>{nome}</span>,
              ])}
            </div>
            <div className="flex flex-col gap-2 text-sm text-[#596061] pt-2">
              <div className="flex items-center gap-3">
                <span className="text-[11px] text-zinc-500 w-16 text-right">lower</span>
                <span className="inline-block h-3 w-48 rounded-sm" style={{ background: GRADIENTE }} />
                <span className="text-[11px] text-zinc-500 w-16">higher</span>
                <span>Colour: z, the list against random genes in that window.</span>
              </div>
              <div className="flex items-center gap-3"><span className="w-6 text-center">★</span>q &lt; 0.05 and it holds when any single donor is left out.</div>
              <div className="flex items-center gap-3"><span className="w-6 text-center">○</span>q &lt; 0.05, but the result rests on one donor.</div>
              <div className="flex items-center gap-3">
                <span className="inline-block w-6 h-4 rounded-sm" style={{ backgroundColor: 'rgba(161,161,170,0.25)', backgroundImage: HACHURA }} />
                Fewer than 3 donors: shown, not tested.
              </div>
            </div>
            <P>
              With log2FC values the card shows three blocks: all genes, the
              up-regulated genes and the down-regulated genes. They can point to
              different stages, which a single score would average away.
            </P>
            <P>
              Development is the main source of variation in brain expression, so
              about 4 in 10 curated gene sets peak somewhere. Compare z values, not only
              the stars: z above 5 places a list among the top ~12% of curated gene sets.
              Hover a cell for its z and q.
            </P>
          </div>
        </div>

        {/* Caveats */}
        <div className="max-w-3xl mb-16">
          <H2>Reading the results</H2>
          <ul className="space-y-3 text-base text-[#596061] leading-relaxed list-disc pl-5">
            <li>
              <Strong>Donors per window</Strong> are 5, 10, 5, 14 and 8 (1st, 2nd and 3rd
              trimester, Infant, Adult). The 3rd trimester is mostly two donors, at 24 and
              37 post-conception weeks, so a 3rd-trimester result often means "near term".
            </li>
            <li>
              <Strong>Untested is not the same as not significant.</Strong> Cells with too
              few donors or marker genes are left out of the tests; the CSV says which.
            </li>
            <li>
              <Strong>Tests run on BrainSpan structures, then are painted on the atlas.</Strong>{' '}
              One structure can colour several parcels (the striatum colours caudate,
              putamen and accumbens), so marks on all of them are one result, not three.
            </li>
            <li>
              <Strong>Every test was checked on random gene lists</Strong>, matched for
              expression level, which rarely produce a mark. The validation scripts are in
              the stats/ folder of the repository.
            </li>
          </ul>
        </div>

        {/* CSV glossary */}
        <div className="max-w-3xl mb-16">
          <H2>Download Data (CSV)</H2>
          <P className="mb-4">
            On the Gene List tab, with a test switched on, the CSV has one row per window
            and atlas parcel, with the map value and the results of all three tests:
          </P>
          <div className="overflow-x-auto border-t border-zinc-200">
            <table className="w-full text-sm text-left">
              <tbody>
                {COLUNAS_CSV.map(([coluna, significado]) => (
                  <tr key={coluna} className="border-b border-zinc-200 align-top">
                    <td className="py-2.5 pr-6 font-mono text-[12px] text-[#2d3435] w-[45%]">{coluna}</td>
                    <td className="py-2.5 text-[#596061] leading-relaxed">{significado}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

    {/*}
        <div className="max-w-2xl mb-16">
            <h2 className="text-sm font-bold text-[#2d3435] uppercase tracking-[0.2em] mb-4">
              Citation
            </h2>
            <p className="text-base text-[#596061] leading-relaxed">
              Vamo usar quando tiver paper{' '}
              <a
                href="link do paper publicado"
                target="_blank"
                rel="noopener noreferrer"
                className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
              >
                Titulo do Paper
              </a>{' '}
              resto da frase.
            </p>
        </div>
    {*/}

        <div className="max-w-2xl">
          <H2>Built by</H2>
          <P>
            This application is developed and maintained by the{' '}
            <a
              href="https://mcblab.github.io/"
              target="_blank"
              rel="noopener noreferrer"
              className="text-[#2d3435] underline decoration-zinc-300 hover:decoration-[#2d3435] transition-colors"
            >
              Marques-Coelho Bioinformatics Lab
            </a>
            , based at the Brain Institute (ICe) of the Federal University of Rio Grande do Norte (UFRN), Natal, Brazil.
          </P>
        </div>

      </div>
    </section>
  );
}
