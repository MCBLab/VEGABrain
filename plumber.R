library(plumber)
library(dplyr)
library(ggplot2)
library(ggseg)
library(GSVA)
library(cowplot)

#* @filter cors
cors <- function(res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  plumber::forward()
}

dados_app <- readRDS("dados_otimizados.rds")
# Médias já agregadas por preparaDados.R (o CSV bruto de ~9.7M linhas não é
# mais lido no boot, o que mantém o startup dentro do limite do Cloud Run)
ontologia_micro <- readRDS("ontologia_micro.rds")
ontologia_macro <- readRDS("ontologia_macro.rds")
genesets_list <- readRDS("genesets_list.rds")

matrix_dados <- dados_app$expression_matrix

meta_limpo <- dados_app$col_meta %>%
  rename(region = structure_mapped) %>%
  filter(!is.na(region), !is.na(broad_age)) %>%
  select(column_num, region, broad_age, macro_region)

# Combinações (idade, região) efetivamente amostradas: usadas para espalhar o
# score da macro-região de volta nas regiões que o ggseg desenha
regioes_por_macro <- meta_limpo %>%
  filter(!is.na(macro_region)) %>%
  distinct(broad_age, region, macro_region)

# 2. LISTAS OTIMIZADAS (À prova de falhas)
setgene <- rownames(matrix_dados)
setontologies <- sort(unique(ontologia_micro$GeneSet))

# Devolve (broad_age, region, Avg_ssGSEA_Score) para uma ontologia. Na escala
# macro cada região recebe a média da macro-região a que pertence.
scores_ontologia <- function(geneset, escala) {
  if (escala == "macro") {
    ontologia_macro %>%
      filter(GeneSet == geneset) %>%
      inner_join(regioes_por_macro,
                 by = c("broad_age", "macro_region"),
                 relationship = "many-to-many") %>%
      select(broad_age, region, Avg_ssGSEA_Score = Score_Medio_ssGSEA)
  } else {
    ontologia_micro %>%
      filter(GeneSet == geneset) %>%
      select(broad_age, region, Avg_ssGSEA_Score = Score_Medio_ssGSEA)
  }
}

# --- Aba Gene List: assinatura simples ou ponderada por log2FC --------------

# Redeclarada dentro de cada endpoint; aqui em cima porque agrega_por_regiao()
# resolve o nome pelo ambiente onde foi definida, nao pelo de quem a chama.
ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)",
                  "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")

# Le a caixa de texto da aba Gene List. Cada linha e "GENE" ou "GENE <log2FC>",
# com espaco, tab, virgula ou ponto-e-virgula entre os dois. Uma linha unica
# com varios simbolos separados por virgula continua valendo (formato antigo).
# tinha_peso marca quem realmente trouxe um numero: e o que decide, mais na
# frente, entre o ssGSEA de sempre e o score ponderado.
parse_assinatura <- function(gene_string) {
  vazio <- data.frame(gene = character(0), peso = numeric(0), tinha_peso = logical(0))
  if (is.null(gene_string) || !nzchar(trimws(gene_string))) return(vazio)

  linhas <- trimws(unlist(strsplit(gene_string, "\r?\n")))
  linhas <- linhas[nzchar(linhas)]
  if (length(linhas) == 0) return(vazio)

  divide <- function(linha) {
    tokens <- unlist(strsplit(linha, "[ \t,;]+"))
    tokens[nzchar(tokens)]
  }

  # Cabecalho de export do DESeq2/edgeR/limma ("gene,log2FoldChange")
  primeiro <- divide(linhas[1])
  if (length(primeiro) >= 2 &&
      grepl("^(gene|genes|symbol|id)$", primeiro[1], ignore.case = TRUE) &&
      !is.finite(suppressWarnings(as.numeric(primeiro[2])))) {
    linhas <- linhas[-1]
  }

  partes <- lapply(linhas, function(linha) {
    tokens <- divide(linha)
    if (length(tokens) == 0) return(NULL)

    if (length(tokens) == 2) {
      valor <- suppressWarnings(as.numeric(tokens[2]))
      if (is.finite(valor)) {
        return(data.frame(gene = tokens[1], peso = valor, tinha_peso = TRUE))
      }
    }
    # Um simbolo so, ou varios sem numero nenhum: lista simples.
    data.frame(gene = tokens, peso = 1, tinha_peso = FALSE)
  })

  do.call(rbind, c(list(vazio), partes))
}

# Resolve a caixa digitada para o nome exato da matriz. Um toupper cego nao
# serve: 1.168 dos 47.808 simbolos tem minuscula (a familia C1orf112), e
# maiusculizar transformaria 1.167 genes hoje validos em "nao encontrado".
# Um unico par colide ignorando a caixa (5S_rRNA / 5S_RRNA); fica com o
# primeiro, e quem digitar a grafia exata cai nela pelo teste de igualdade.
indice_genes <- setgene[!duplicated(toupper(setgene))]
names(indice_genes) <- toupper(indice_genes)

# Parse + confronto com a matriz. Concentra tudo o que os tres endpoints da
# aba Gene List precisam saber sobre a lista submetida.
assinatura_valida <- function(gene_string) {
  bruta <- parse_assinatura(gene_string)

  # Grafia exata primeiro; so o que nao bater cai no casamento sem caixa.
  exato <- bruta$gene %in% setgene
  canonico <- as.character(ifelse(exato, bruta$gene,
                                  unname(indice_genes[toupper(bruta$gene)])))
  no_banco <- !is.na(canonico)

  validos <- bruta[no_banco, , drop = FALSE]
  validos$gene <- canonico[no_banco]
  validos <- validos[!duplicated(validos$gene), , drop = FALSE]

  ponderado <- any(validos$tinha_peso)

  list(
    validos    = validos,
    # O texto original, para a pessoa achar o simbolo na propria lista.
    invalidos  = unique(bruta$gene[!no_banco]),
    ponderado  = ponderado,
    # Sem log2FC nenhum a assinatura nao tem direcao: nao ha up nem down.
    n_up       = if (ponderado) sum(validos$peso > 0) else 0L,
    n_down     = if (ponderado) sum(validos$peso < 0) else 0L,
    n_sem_peso = if (ponderado) sum(!validos$tinha_peso) else 0L
  )
}

# score_j = sum_i(w_i * z_ij) / sum_i(|w_i|), com z = expressao padronizada
# entre as amostras. Genes sem variancia nao informam nada e saem dos dois
# lados da fracao: mante-los no denominador achataria o mapa inteiro em direcao
# ao zero.
score_ponderado <- function(genes, pesos) {
  sub <- matrix_dados[genes, , drop = FALSE]
  desv <- apply(sub, 1, sd)
  ok <- is.finite(desv) & desv > 0
  if (!any(ok)) stop("Error: none of the submitted genes varies across samples.")

  sub <- sub[ok, , drop = FALSE]
  pesos <- pesos[ok]
  if (sum(abs(pesos)) == 0) stop("Error: all submitted log2FC values are zero.")

  z <- (sub - rowMeans(sub)) / desv[ok]
  as.numeric(crossprod(z, pesos) / sum(abs(pesos)))
}

# Recebe um score por amostra e devolve a media por (janela, regiao). Na escala
# macro cada regiao desenhada recebe a media da macro-regiao a que pertence.
agrega_por_regiao <- function(score_por_amostra, escala, nome_col,
                              amostras = colnames(matrix_dados)) {
  dados_base <- data.frame(
    column_num = as.numeric(amostras),
    Score = score_por_amostra
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  if (escala == "macro") {
    dados_prontos <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Score_Medio = mean(Score, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Score_Medio)
  } else {
    dados_prontos <- dados_base %>%
      group_by(broad_age, region) %>%
      summarise(Score_Medio = mean(Score, na.rm = TRUE), .groups = "drop")
  }

  names(dados_prontos)[names(dados_prontos) == "Score_Medio"] <- nome_col
  dados_prontos
}

# Score ponderado quando a lista trouxe log2FC, ssGSEA de sempre quando nao.
scores_assinatura <- function(assinatura, escala) {
  genes <- assinatura$validos$gene
  if (length(genes) == 0) {
    stop("Error: none of the submitted genes was found in the expression matrix.")
  }

  if (assinatura$ponderado) {
    coluna <- "Avg_Weighted_Score"
    dados <- agrega_por_regiao(
      score_ponderado(genes, assinatura$validos$peso), escala, coluna)
  } else {
    coluna <- "Avg_ssGSEA_Score"
    ssgsea_parametros <- GSVA::ssgseaParam(
      exprData = as.matrix(matrix_dados),
      geneSets = list(UserSet = genes)
    )
    ssgsea_resultado <- GSVA::gsva(ssgsea_parametros, verbose = FALSE)
    dados <- agrega_por_regiao(
      as.numeric(ssgsea_resultado["UserSet", ]), escala, coluna,
      colnames(ssgsea_resultado))
  }

  list(dados = dados, coluna = coluna)
}
#* All gene symbols in the expression matrix
#* @get /list_genes
function() {
  setgene
}

#* All precomputed gene set names
#* @get /list_ontologies
function() {
  setontologies
}

# Helper function to create the lateral and medial brain plots grid
build_brain_grid <- function(data, fill_var, fill_scale, caption_text = NULL) {
  p_lateral <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  p_medial <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::dk(),
      position = position_brain(c("right medial")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

    p_sagittal <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::aseg(),
      position = position_brain(c("sagittal")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

    p_coronal <- data %>%
    ggplot() +
    geom_brain(
      atlas = ggseg::aseg(),
      position = position_brain(c("coronal_1")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    facet_wrap(~broad_age, ncol = 5) +
    labs(caption = caption_text) +
    theme_void() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(color = "firebrick", face = "bold", size = 11, hjust = 0.5, margin = margin(t = 15))
    )

  plot_grid(p_lateral, p_medial, p_sagittal, p_coronal, nrow = 4, align = "v", rel_heights = c(1, 0.93, 1.12, 1.43))  #1, 0.90, 1.10, 1.40
}

# 3. ENDPOINTS DE PLOTAGEM (SINTAXE OFICIAL MODERNA DO GGSEG)

# Cada mapa e servido em dois formatos: SVG (vetorial, para publicacao) e PNG
# (o que a tela exibe - o SVG do ggseg carrega milhares de poligonos e pesa no
# navegador). O calculo fica nas funcoes grafico_*, e cada endpoint so escolhe
# o device do serializer, para que os dois formatos nunca divirjam.

grafico_gene <- function(gene = "SOX10", escala = "micro") {
  if (!(gene %in% rownames(matrix_dados))) stop("Error: gene not found.")

  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")
  matrix_dados <- as.matrix(dados_app$expression_matrix)

  dados_base <- data.frame(
    Amostra = as.numeric(colnames(matrix_dados)),
    Expressao = as.numeric(matrix_dados[gene, ])
  ) %>%
    left_join(dados_app$col_meta, by = c("Amostra" = "column_num")) %>%
    rename(region = structure_mapped) %>%
    filter(!is.na(region), !is.na(broad_age)) %>%
    mutate(Expressao = ifelse(Expressao > 6, 6, Expressao)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  if(escala == "macro"){
    dados_gene_plot <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      mutate(Avg_Expression = mean(Expressao, na.rm = TRUE)) %>%
      ungroup() %>%
      distinct(broad_age, region, Avg_Expression)
  } else {
    dados_gene_plot <- dados_base %>%
    group_by(broad_age, region) %>%
    summarise(Avg_Expression = mean(Expressao, na.rm = TRUE), .groups = "drop")
  }

  escala <- scale_fill_gradientn( colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"), 
    values = scales::rescale(c(0, 1.5, 3, 4.5, 6)),limits = c(0, 6), 
    name = "Log2 Expr", na.value = "darkgray")  

  build_brain_grid(dados_gene_plot, "Avg_Expression", escala)
}

#* Brain map of one gene's mean log2(RPKM+1) expression (SVG)
#* @param gene Gene symbol
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer svg list(width = 11, height = 7)
#* @get /plot_brain
function(gene = "SOX10", scale = "micro") {
  print(grafico_gene(gene, scale))
}

#* Same map as /plot_brain, as PNG
#* @param gene Gene symbol
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer png list(width = 11, height = 7, units = "in", res = 150)
#* @get /plot_brain_png
function(gene = "SOX10", scale = "micro") {
  print(grafico_gene(gene, scale))
}

grafico_ontologia <- function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", escala = "micro") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")

  dados_prontos <- scores_ontologia(geneset, escala)

  if (nrow(dados_prontos) == 0) stop("Error: gene set not found or has no data to plot.")

  dados_prontos <- dados_prontos %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))

  # Ontology scores span a narrow range (e.g. 0.29-0.39), where the default
  # 0.025-step breaks need three decimals and the labels run together.
  # breaks_pretty only uses 1/2/5 steps, so labels get two decimals whenever
  # the range allows and stay exact (0.325 is never printed as "0.32"); the
  # wider bar fits up to five labels.
  escala_cor <- list(
    scale_fill_gradientn( colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"),
      breaks = scales::breaks_pretty(n = 4),
      name = "ssGSEA Score", na.value = "darkgray"),
    guides(fill = guide_colorbar(barwidth = unit(5, "cm")))
  )

  build_brain_grid(dados_prontos, "Avg_ssGSEA_Score", escala_cor)
}

#* Brain map of a precomputed gene set's mean ssGSEA score (SVG)
#* @param geneset Gene set name, as listed by /list_ontologies
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer svg list(width = 11, height = 7)
#* @get /plot_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", scale = "micro") {
  print(grafico_ontologia(geneset, scale))
}

#* Same map as /plot_ontology, as PNG
#* @param geneset Gene set name, as listed by /list_ontologies
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer png list(width = 11, height = 7, units = "in", res = 150)
#* @get /plot_ontology_png
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", scale = "micro") {
  print(grafico_ontologia(geneset, scale))
}

grafico_genelist <- function(gene_string = "", escala = "micro") {
  assinatura <- assinatura_valida(gene_string)
  resultado <- scores_assinatura(assinatura, escala)

  if (assinatura$ponderado) {
    # O score soma zero entre as amostras por construcao, entao limites
    # simetricos deixam o tom claro da paleta exatamente no zero: azul e down,
    # vinho e up. Sem isso o mapa nao distingue as duas direcoes.
    limite <- max(abs(resultado$dados[[resultado$coluna]]), na.rm = TRUE)
    if (!is.finite(limite) || limite == 0) limite <- 1

    escala_cor <- list(
      scale_fill_gradientn(
        colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"),
        limits = c(-limite, limite),
        breaks = c(-limite, 0, limite),
        labels = scales::label_number(accuracy = 0.01),
        name = "Weighted Score", na.value = "darkgray"),
      # Com as cinco marcas padrao o sinal de menos faz os rotulos colidirem
      # na barra estreita; tres marcas numa barra mais larga cabem.
      guides(fill = guide_colorbar(barwidth = unit(5, "cm")))
    )

    # Um mapa ponderado e um nao ponderado sao indistinguiveis depois de
    # baixados; a legenda carimba o metodo e a composicao da assinatura.
    legenda <- sprintf("Weighted signature - %d genes (%d up / %d down)",
                       nrow(assinatura$validos), assinatura$n_up, assinatura$n_down)
    n_fora <- length(assinatura$invalidos)
    if (n_fora > 0) {
      legenda <- paste0(legenda, sprintf(" - %d symbol%s not found",
                                         n_fora, if (n_fora == 1) "" else "s"))
    }
  } else {
    escala_cor <- scale_fill_gradientn( colors = c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845"), 
      #values = scales::rescale(c(0, 2, 4, 6, 8)),limits = c(0, 8), 
      name = "Custom ssGSEA", na.value = "darkgray")
    legenda <- NULL
  }

  build_brain_grid(resultado$dados, resultado$coluna, escala_cor, legenda)
}

#* Brain map of a user gene list, scored on demand: ssGSEA, or a log2FC-weighted score when values are given (SVG)
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer svg list(width = 11, height = 7)
#* @post /plot_genelist
function(gene_string = "", scale = "micro") {
  print(grafico_genelist(gene_string, scale))
}

#* Same map as /plot_genelist, as PNG
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer png list(width = 11, height = 7, units = "in", res = 150)
#* @post /plot_genelist_png
function(gene_string = "", scale = "micro") {
  print(grafico_genelist(gene_string, scale))
}

#* Table behind /plot_brain: mean expression per window and region (CSV)
#* @param gene Gene symbol
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer csv
#* @get /data_brain
function(gene = "SOX10", scale = "micro") {
  
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")
  
  dados_base <- data.frame(
    column_num = as.numeric(colnames(matrix_dados)), 
    Expressao = as.numeric(matrix_dados[gene, ])
  ) %>%
    inner_join(meta_limpo, by = "column_num") %>%
    filter(!is.na(broad_age)) %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades))
    
  if (scale == "macro") {
    dados_finais <- dados_base %>%
      filter(!is.na(macro_region)) %>%
      group_by(broad_age, macro_region) %>%
      summarise(Avg_Expression = mean(Expressao, na.rm = TRUE), .groups = "drop")
  } else {
    dados_finais <- dados_base %>%
      group_by(broad_age, region) %>%
      summarise(Avg_Expression = mean(Expressao, na.rm = TRUE), .groups = "drop")
  }

  return(dados_finais)
}

#* Table behind /plot_ontology: mean ssGSEA score per window and region (CSV)
#* @param geneset Gene set name, as listed by /list_ontologies
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer csv
#* @get /data_ontology
function(geneset = "GOBP_FOREBRAIN_GENERATION_OF_NEURONS", scale = "micro") {
  ordem_idades <- c("1st trimester (n = 5)", "2nd trimester (n = 10)", "3rd trimester (n = 5)", "Infant (n = 14)", "Adult (n = 8)")

  dados_prontos <- scores_ontologia(geneset, scale)

  if (nrow(dados_prontos) == 0) stop("Error: gene set not found or has no data to plot.")

  dados_prontos <- dados_prontos %>%
    mutate(broad_age = factor(broad_age, levels = ordem_idades)) %>%
    arrange(broad_age, region)

  return(dados_prontos)
}

#* Table behind /plot_genelist: mean score per window and region (CSV)
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer csv
#* @post /data_genelist
function(gene_string = "", scale = "micro") {
  assinatura <- assinatura_valida(gene_string)

  if (nrow(assinatura$validos) == 0) {
    return(data.frame(Error = "No valid genes"))
  }

  return(scores_assinatura(assinatura, scale)$dados)
}

#* Checks which submitted symbols exist in the expression matrix, and whether the list carries log2FC weights
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @post /validate_genelist
function(gene_string = "") {
  assinatura <- assinatura_valida(gene_string)

  # unbox: sem isso o serializer manda [false] / [3] no lugar de escalares,
  # e o front teria de desembrulhar campo a campo.
  return(list(
    valid        = assinatura$validos$gene,
    invalid      = assinatura$invalidos,
    weighted     = jsonlite::unbox(assinatura$ponderado),
    n_up         = jsonlite::unbox(assinatura$n_up),
    n_down       = jsonlite::unbox(assinatura$n_down),
    n_unweighted = jsonlite::unbox(assinatura$n_sem_peso)
  ))
}

#* Member genes of a precomputed gene set that are present in the expression matrix
#* @param geneset Gene set name, as listed by /list_ontologies
#* @get /geneset_genes
function(geneset = "") {
  if (geneset == "" || is.null(genesets_list[[geneset]])) {
    return(list(genes = character(0)))
  }
  
  genes_da_via <- genesets_list[[geneset]]
  genes_presentes <- intersect(genes_da_via, rownames(matrix_dados))
  
  return(list(genes = as.character(genes_presentes)))
}

#* Reference map naming the atlas regions (micro) or anatomical groups (macro) of one view (SVG)
#* @param view Which view to draw: lateral, medial, sagittal or coronal
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @serializer svg list(width = 9, height = 7)
#* @get /plot_reference
function(view = "lateral", scale = "micro") {
  
  if (scale == "macro") {
    mapa_referencia <- dados_app$mapping_info$region_to_macro_structure
    coluna <- "macro_region"
  } else {
    mapa_referencia <- dados_app$mapping_info$region_to_structure
    coluna <- "structure_name"
  }

  mapa_referencia <- mapa_referencia %>%
    filter(!is.na(.data[[coluna]]))
  
  todas_areas <- sort(unique(mapa_referencia[[coluna]]))
  
  mapa_referencia[[coluna]] <- factor(mapa_referencia[[coluna]], levels = todas_areas)
  
  escala_sem_NA <- scale_fill_discrete(drop = TRUE, na.translate = FALSE, name = "")

  if (view == "lateral") {
    df_dk <- as.data.frame(ggseg::dk())
    areas <- df_dk$region[df_dk$hemi == "right" & df_dk$view == "lateral"]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::dk(), position = position_brain("right lateral"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Lateral Cortex -", ifelse(scale == "macro", "Macro", "Micro")))
    
  } else if (view == "medial") {
    df_dk <- as.data.frame(ggseg::dk())
    areas <- df_dk$region[df_dk$hemi == "right" & df_dk$view == "medial"]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::dk(), position = position_brain("right medial"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Medial Cortex -", ifelse(scale == "macro", "Macro", "Micro")))
    
  } else if (view == "sagittal") {
    df_aseg <- as.data.frame(ggseg::aseg())
    areas <- df_aseg$region[grepl("sagittal", df_aseg$view, ignore.case = TRUE)]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::aseg(), position = position_brain("sagittal"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Subcortical (Sagittal) -", ifelse(scale == "macro", "Macro", "Micro")))
    
  } else if (view == "coronal") {
    df_aseg <- as.data.frame(ggseg::aseg())
    areas <- df_aseg$region[grepl("coronal_1", df_aseg$view, ignore.case = TRUE)]
    dados_plot <- mapa_referencia %>% filter(region %in% areas)
    
    p <- ggplot(dados_plot) + geom_brain(atlas = ggseg::aseg(), position = position_brain("coronal_1"), mapping = aes(fill = .data[[coluna]]), color = "black", size = 0.3) + ggtitle(paste("Subcortical (Coronal) -", ifelse(scale == "macro", "Macro", "Micro")))
    
  } else {
    stop("Error: unknown view; use lateral, medial, sagittal or coronal.")
  }

  # 4. Finalização do Plot
  p_final <- p + 
    escala_sem_NA + 
    theme_void() + 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    guides(fill = guide_legend(ncol = 2, keywidth = 1, keyheight = 1))

  print(p_final)
}