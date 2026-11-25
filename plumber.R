library(plumber)
library(dplyr)
library(ggplot2)
library(ggseg)
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

# --- Per-area significance (Gene List tab) -----------------------------------
# Two complementary tests, one row per (BrainSpan structure or class, window):
#   score    does the list's ssGSEA score peak here beyond expression-matched
#            random gene sets? (stats/enrichment_stats.R)
#   markers  are this area's specific genes over-represented in the list?
#            (stats/marker_stats.R)
# Nothing is precomputed: everything comes from the expression matrix and its
# sample metadata. The large tables are built on first use, not at boot, because
# Cloud Run times the cold start.
source("stats/enrichment_stats.R")
source("stats/marker_stats.R")
# Plus a stage-level test: which developmental window the list peaks in
# (stats/stage_stats.R). The area tests cannot see a rise shared by every
# structure of a window; this one asks exactly that.
source("stats/stage_stats.R")

score_ref <- cells_from_meta(dados_app$col_meta)
stage_ref <- stage_meta(dados_app$col_meta)
# Null genes and marker backgrounds come from the MSigDB symbols: the kind of
# gene a curated or user list is made of, not the ~30k unannotated matrix rows.
universe_genes <- intersect(unique(unlist(genesets_list, use.names = FALSE)), setgene)

stats_cache <- new.env()
# Rank tables (~300 MB, ~6 s): ssGSEA of any list in ~1 s, identical to GSVA.
rank_tables_cached <- function() {
  if (is.null(stats_cache$rt)) stats_cache$rt <- rank_tables(matrix_dados, universe = universe_genes)
  stats_cache$rt
}
# z-scores (~200 MB) for the null of log2FC-weighted signatures.
weighted_reference_cached <- function() {
  if (is.null(stats_cache$wref)) stats_cache$wref <- weighted_reference(matrix_dados, universe = universe_genes)
  stats_cache$wref
}
# Marker genes per (area, window) (~2 s, < 50 MB).
marker_reference <- function() {
  if (is.null(stats_cache$mk)) stats_cache$mk <- area_markers(matrix_dados, dados_app$col_meta, universe_genes)
  stats_cache$mk
}

# The map, the table and the CSV of one list each ask for the same test; the
# results are deterministic (fixed seeds), so the last few are kept.
results_cache <- new.env()
remember <- function(key, value) {
  if (!is.null(results_cache[[key]])) return(results_cache[[key]])
  if (length(ls(results_cache)) >= 32) rm(list = ls(results_cache), envir = results_cache)
  assign(key, value, envir = results_cache)  # `value` is evaluated here, once
  value
}

parse_significance <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (length(x) == 0 || x %in% c("", "none", "false", "off", "0")) return("none")
  if (x %in% c("score", "true", "1")) return("score")
  if (x == "markers") return("markers")
  stop("Error: significance must be none, score or markers.")
}

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
    # GSVA's ssGSEA computed algebraically from the rank tables (identical to
    # GSVA::gsva() to 1e-6, see stats/validate_stats.R check 0) in ~1 s instead
    # of ~45 s. GSVA normalises by the score's own range; done here the same way.
    rt <- rank_tables_cached()
    per_sample <- ssgsea_from_ranks(rt, membership(match(genes, rt$genes), rt$n))[1, ]
    dados <- agrega_por_regiao(per_sample / diff(range(per_sample)), escala, coluna,
                               names(per_sample))
  }

  list(dados = dados, coluna = coluna)
}

# Runs one significance test on a validated signature. Returns
# list(test, cells, reason): `cells` is the per-cell table of enrich(),
# enrich_weighted() or area_enrichment(); `reason` says why nothing was tested.
run_significance <- function(assinatura, test) {
  genes <- assinatura$validos$gene
  weights <- if (assinatura$ponderado) assinatura$validos$peso else NULL
  key <- paste(test, paste(genes, if (is.null(weights)) "" else weights, sep = ":", collapse = ","))
  remember(key, tryCatch({
    if (test == "score") {
      if (length(genes) < 5) {
        stop(sprintf("needs at least 5 genes found in the data (got %d)", length(genes)))
      }
      cells <- if (is.null(weights)) {
        enrich(score_ref, rank_tables_cached(), genes)
      } else {
        enrich_weighted(score_ref, weighted_reference_cached(), genes, weights)
      }
    } else {
      cells <- area_enrichment(marker_reference(), genes, weights)
    }
    list(test = test, cells = cells, reason = NULL)
  }, error = function(e) list(test = test, cells = NULL, reason = conditionMessage(e))))
}

# Region-level verdicts: a parcel inherits the verdict of the structure(s) or
# class that feed it. Score: one row per (window, region). Markers: one row per
# (window, region, direction), direction being all, up or down.
painted_significance <- function(sig) {
  mapping <- dados_app$mapping_info$region_to_structure
  if (sig$test == "score") paint_outward(sig$cells, mapping) else paint_markers(sig$cells, mapping)
}

# One row per marked parcel and window; `kind` picks the glyph: score (*), or
# all / up / down for the marker test.
significance_marks <- function(sig) {
  if (is.null(sig) || is.null(sig$cells)) return(NULL)
  painted <- painted_significance(sig) %>% filter(enriched)
  marks <- if (sig$test == "score") {
    painted %>% distinct(broad_age, region) %>% mutate(kind = "score")
  } else {
    painted %>% distinct(broad_age, region, kind = direction)
  }
  marks %>% mutate(broad_age = factor(broad_age, levels = ordem_idades))
}

# Stamped on the map, so a downloaded image says on its own what a mark means.
significance_caption <- function(sig) {
  if (is.null(sig)) return(NULL)
  if (!is.null(sig$reason)) return(paste0("Significance not tested: ", sig$reason))
  cells <- sig$cells
  if (sig$test == "score") {
    return(sprintf(paste0("* q < 0.05 vs 4,000 expression-matched random gene sets (BH over %d structure x window cells),\n",
                          "robust to leaving any single donor out; %d enriched cells"),
                   sum(cells$tested), sum(cells$enriched)))
  }
  weighted <- any(cells$direction != "all")
  skipped <- attr(cells, "skipped")
  sprintf(paste0("%s over-represented among area-specific genes (Fisher, BH over %d cells, q < 0.05, >= %d genes);\n",
                 "region classes vs each other, neocortical areas vs the rest of neocortex; %d enriched cells%s"),
          if (weighted) "Triangles up / down: up- / down-regulated genes" else "Diamonds: list",
          sum(cells$tested), MIN_OVERLAP, sum(cells$enriched_marker),
          if (length(skipped)) sprintf(" (%s genes not tested: fewer than %d)", paste(skipped, collapse = ", "),
                                       MIN_LIST_GENES) else "")
}

# Region-level columns appended to /data_genelist.
significance_columns <- function(sig) {
  if (is.null(sig$cells)) return(NULL)
  painted <- painted_significance(sig)
  if (sig$test == "score") {
    out <- painted %>%
      transmute(broad_age, region, tested_score = tested, enriched_score = enriched,
                q_score = q_competitive, p_loo)
  } else {
    out <- painted %>%
      transmute(broad_age, region, direction,
                enriched_markers = enriched, q_markers = q_marker, markers_level = level,
                markers_genes = genes) %>%
      tidyr::pivot_wider(names_from = direction,
                         values_from = c(enriched_markers, q_markers, markers_level, markers_genes),
                         names_glue = "{.value}_{direction}")
    names(out) <- sub("_all$", "", names(out))
  }
  out %>% mutate(broad_age = factor(broad_age, levels = ordem_idades))
}

# Stage test of a validated signature (stats/stage_stats.R), cached like the
# area tests: the panel and the CSV of one list ask for it twice.
run_stage <- function(assinatura) {
  genes <- assinatura$validos$gene
  weights <- if (assinatura$ponderado) assinatura$validos$peso else NULL
  key <- paste("stage", paste(genes, if (is.null(weights)) "" else weights, sep = ":", collapse = ","))
  remember(key, tryCatch(list(rows = stage_enrichment(rank_tables_cached(), stage_ref, genes, weights),
                              reason = NULL),
                         error = function(e) list(rows = NULL, reason = conditionMessage(e))))
}

# Drawn parcel -> PsychENCODE stage region (NCX, HIP, ...), through the region
# class of the structure(s) feeding it. Every parcel has exactly one class.
stage_region_of <- dados_app$mapping_info$region_to_structure %>%
  mutate(stage_region = names(STAGE_REGIONS)[match(unname(STRUCTURE_CLASS[structure_name]), STAGE_REGIONS)]) %>%
  distinct(region, stage_region)

# Stage columns for /data_genelist, one set per row of the map data: the test of
# the parcel's PsychENCODE region (stage_*) and of the whole brain (stage_brain_*)
# in that window. log2FC lists add _up / _down copies for the split lists.
stage_columns <- function(dados, rows) {
  sfx <- c(all = "", up = "_up", down = "_down")
  widen <- function(d, keys, prefix) {
    vals <- setdiff(names(d), c(keys, "list"))
    d %>% mutate(list = sfx[list]) %>%
      tidyr::pivot_wider(id_cols = all_of(keys), names_from = list, values_from = all_of(vals),
                         names_glue = paste0(prefix, "{.value}{list}"))
  }
  regional <- rows %>% filter(level == "region") %>%
    transmute(list, stage_region = region, window, tested, z = round(z, 3), q, p_loo, enriched) %>%
    widen(c("stage_region", "window"), "stage_")
  brain <- rows %>% filter(level == "brain") %>%
    transmute(list, window, z = round(z, 3), q, p_loo, enriched) %>%
    widen("window", "stage_brain_")
  dados %>%
    mutate(window = window_short(as.character(broad_age))) %>%
    left_join(stage_region_of, by = "region") %>%
    left_join(regional, by = c("stage_region", "window")) %>%
    left_join(brain, by = "window") %>%
    select(-window)
}

# Tested cells, most significant first, for /area_enrichment.
significance_rows <- function(sig) {
  cells <- sig$cells
  if (is.null(cells)) return(data.frame())
  if (sig$test == "score") {
    cells %>%
      filter(tested) %>%
      arrange(q_competitive) %>%
      transmute(window = broad_age, area = structure_original, level = "structure",
                n_donors, p = signif(p_competitive, 3), q = signif(q_competitive, 3),
                p_loo = signif(p_loo, 3), enriched, underpowered_window)
  } else {
    cells %>%
      filter(tested) %>%
      arrange(q_marker, p_marker) %>%
      transmute(window = broad_age, area = cell, level, direction, n_donors, n_markers,
                n_list, overlap, expected = round(expected, 2),
                # an infinite odds ratio (every list gene is a marker) has no JSON number
                odds_ratio = ifelse(is.finite(odds_ratio), signif(odds_ratio, 3), NA),
                or_ci_low = signif(or_ci_low, 3), p = signif(p_marker, 3), q = signif(q_marker, 3),
                enriched = enriched_marker, genes = overlap_genes, structures)
  }
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

# Where a mark goes: a point guaranteed to lie INSIDE the largest polygon of each
# parcel (sf::st_point_on_surface; the centroid falls outside thin, curved
# parcels such as superior parietal in the medial view), computed on the
# geometry ggseg actually draws in that view -- one branch per ggseg version
# (Docker pins 2.1.1; 2.2 changed the internals). Hemispheres drawn on top of
# each other (the cerebellum in the sagittal view) collapse to one point.
parcel_centroids <- function(atlas, pos) {
  measure <- function(pg, region, hemi, view) {
    bb <- sf::st_bbox(pg)
    pt <- tryCatch(sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(pg))),
                   error = function(e) sf::st_coordinates(sf::st_centroid(pg)))
    tibble(region = region, hemi = hemi, view = view, area = as.numeric(sf::st_area(pg)),
           x = pt[1, 1], y = pt[1, 2],
           xmin = bb[["xmin"]], xmax = bb[["xmax"]], ymin = bb[["ymin"]], ymax = bb[["ymax"]])
  }
  if (exists("prepare_polygon_atlas", envir = asNamespace("ggseg"), inherits = FALSE)) {
    # ggseg >= 2.2 draws flat polygons (x, y)
    flat <- ggseg:::prepare_polygon_atlas(atlas, position = position_brain(pos)) %>%
      filter(!is.na(region)) %>%
      mutate(key = paste(region, hemi, .feature_id, subgroup, sep = "|"))
    rings <- bind_rows(lapply(split(flat, flat$key), function(d) {
      xy <- cbind(d$x, d$y)
      if (any(xy[1, ] != xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
      measure(sf::st_sfc(sf::st_polygon(list(xy))), d$region[1], d$hemi[1], d$view[1])
    }))
  } else {
    # ggseg 2.1.x draws repositioned sf geometry
    rep <- ggseg::reposition_brain(as.data.frame(atlas), position = pos)
    rep <- rep[!is.na(rep$region), ]
    geo <- sf::st_geometry(rep)
    rings <- bind_rows(lapply(seq_len(nrow(rep)), function(i) {
      polys <- suppressWarnings(sf::st_cast(geo[i], "POLYGON"))
      bind_rows(lapply(seq_along(polys), function(k) measure(polys[k], rep$region[i], rep$hemi[i], rep$view[i])))
    }))
  }
  # a parcel that is only a sliver in this view (precentral seen from the medial
  # side: 0.05% of the panel) is marked where it is actually drawn, not here
  view_area <- (max(rings$xmax) - min(rings$xmin)) * (max(rings$ymax) - min(rings$ymin))
  rings %>%
    group_by(region, hemi, view) %>%
    slice_max(area, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(area >= 0.0025 * view_area) %>%
    group_by(region, view) %>%
    # keep both hemispheres only when their bounding boxes are disjoint
    filter(n() == 1 | max(xmin) > min(xmax) | max(ymin) > min(ymax) | area == max(area)) %>%
    ungroup() %>%
    select(region, hemi, view, x, y, xmin, xmax)
}

# Significance marks for one view: `marks` has (broad_age, region, kind), one row
# per marked parcel and window. Score test: a white asterisk on a dark halo.
# Marker test: a white diamond, or triangles pointing up / down for the two halves
# of a log2FC signature, nudged apart so both fit. White on a dark outline reads
# over any tone of the palette. Fill is a constant, never an aesthetic, so the
# layer cannot interfere with the map's fill scale. NULL when there is nothing to
# mark (adding NULL to a ggplot does nothing).
significance_layer <- function(marks, atlas, pos, size = 9) {
  if (is.null(marks) || nrow(marks) == 0) return(NULL)
  centroids <- parcel_centroids(atlas, pos)
  points <- marks %>% inner_join(centroids, by = "region", relationship = "many-to-many")
  if (nrow(points) == 0) return(NULL)

  layers <- list()
  stars <- points %>% filter(kind == "score")
  if (nrow(stars) > 0) {
    layers <- c(layers, list(
      geom_text(data = stars, aes(x = x, y = y), inherit.aes = FALSE, label = "*",
                size = size * 1.2, colour = "grey20", fontface = "bold", vjust = 0.3),
      geom_text(data = stars, aes(x = x, y = y), inherit.aes = FALSE, label = "*",
                size = size, colour = "white", fontface = "bold", vjust = 0.3)))
  }
  nudge <- 0.012 * (max(centroids$xmax) - min(centroids$xmin))
  for (k in c("all", "up", "down")) {
    pts <- points %>% filter(kind == k)
    if (nrow(pts) == 0) next
    pts$x <- pts$x + switch(k, all = 0, up = -nudge, down = nudge)
    layers <- c(layers, list(
      geom_point(data = pts, aes(x = x, y = y), inherit.aes = FALSE,
                 shape = switch(k, all = 23, up = 24, down = 25),
                 size = size / 3, fill = "white", colour = "grey15", stroke = 0.9)))
  }
  layers
}

# Lateral and medial cortex (dk), sagittal and coronal subcortex (aseg), one
# column per window. The data go to geom_brain explicitly, grouped by the facet:
# ggseg >= 2.2 no longer reads the plot's data, and on 2.1.x the result is the same.
build_brain_grid <- function(data, fill_var, fill_scale, caption_text = NULL, marks = NULL,
                             coronal_height = 1.43) {
  data <- group_by(data, broad_age)
  p_lateral <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::dk(),
      position = position_brain(c("right lateral")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    significance_layer(marks, ggseg::dk(), "right lateral") +
    facet_wrap(~broad_age, ncol = 5) +
    theme_void() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 5, b = 0, l = 5),
      legend.key.height = unit(2, "cm"),
      plot.background = element_rect(fill = "white", colour = NA)
    )

  p_medial <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::dk(),
      position = position_brain(c("right medial")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    significance_layer(marks, ggseg::dk(), "right medial") +
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

    p_sagittal <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::aseg(),
      position = position_brain(c("sagittal")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    significance_layer(marks, ggseg::aseg(), "sagittal") +
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

    p_coronal <- ggplot() +
    geom_brain(
      data = data,
      atlas = ggseg::aseg(),
      position = position_brain(c("coronal_1")),
      mapping = aes(fill = .data[[fill_var]]),
      color = 'black',
      size = 0.50
    ) +
    fill_scale +
    # small, adjacent parcels: a smaller glyph so marks do not overlap
    significance_layer(marks, ggseg::aseg(), "coronal_1", size = 6) +
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

  # The coronal row also holds the legend and caption. The Gene List map carries
  # a two-line caption, so it asks for a taller row (1.9, with a 7.6 in image) to
  # keep its slices as large as the Gene and Ontology ones at 1.43 / 7 in.
  plot_grid(p_lateral, p_medial, p_sagittal, p_coronal, nrow = 4, align = "v",
            rel_heights = c(1, 0.93, 1.12, coronal_height))
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

grafico_genelist <- function(gene_string = "", escala = "micro", significance = "none") {
  assinatura <- assinatura_valida(gene_string)
  resultado <- scores_assinatura(assinatura, escala)
  significance <- parse_significance(significance)
  sig <- if (significance == "none") NULL else run_significance(assinatura, significance)

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

  # The significance note joins the same caption, so a downloaded image explains
  # its own marks.
  legenda <- c(legenda, significance_caption(sig))
  legenda <- if (length(legenda) == 0) NULL else paste(legenda, collapse = "\n")

  build_brain_grid(resultado$dados, resultado$coluna, escala_cor, legenda,
                   marks = significance_marks(sig), coronal_height = 1.9)
}

#* Brain map of a user gene list, scored on demand: ssGSEA, or a log2FC-weighted score when values are given (SVG)
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @param significance Mark significant areas: none, score (asterisk; list score vs random gene sets) or markers (diamond / triangles; list vs area-specific genes)
#* @serializer svg list(width = 11, height = 7.6)
#* @post /plot_genelist
function(gene_string = "", scale = "micro", significance = "none") {
  print(grafico_genelist(gene_string, scale, significance))
}

#* Same map as /plot_genelist, as PNG
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @param scale Region granularity: micro (atlas regions) or macro (anatomical groups)
#* @param significance Mark significant areas: none, score (asterisk; list score vs random gene sets) or markers (diamond / triangles; list vs area-specific genes)
#* @serializer png list(width = 11, height = 7.6, units = "in", res = 150)
#* @post /plot_genelist_png
function(gene_string = "", scale = "micro", significance = "none") {
  print(grafico_genelist(gene_string, scale, significance))
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
#* @param significance Add per-region significance columns: none, score (tested_score, enriched_score, q_score, p_loo), markers (enriched_markers, q_markers, markers_level, markers_genes; _up / _down for log2FC lists) or all (both, plus the stage test of the parcel's PsychENCODE region and of the whole brain: stage_region, stage_tested, stage_z, stage_q, stage_p_loo, stage_enriched, stage_brain_*)
#* @serializer csv
#* @post /data_genelist
function(gene_string = "", scale = "micro", significance = "none") {
  assinatura <- assinatura_valida(gene_string)

  if (nrow(assinatura$validos) == 0) {
    return(data.frame(Error = "No valid genes"))
  }

  dados <- scores_assinatura(assinatura, scale)$dados
  every <- identical(tolower(trimws(as.character(significance))), "all")
  significance <- if (every) "all" else parse_significance(significance)
  if (significance == "none") return(dados)

  # Same rule as the map: a region inherits the verdict of the structure(s) or
  # class feeding it, with the smallest q among them.
  for (t in if (every) c("score", "markers") else significance) {
    sig <- run_significance(assinatura, t)
    dados <- if (is.null(sig$cells)) {
      mutate(dados, "{t}_note" := paste("Not tested:", sig$reason))
    } else {
      left_join(dados, significance_columns(sig), by = c("broad_age", "region"))
    }
  }
  if (!every) return(dados)
  st <- run_stage(assinatura)
  if (is.null(st$rows)) return(mutate(dados, stage_note = paste("Not tested:", st$reason)))
  stage_columns(dados, st$rows)
}

#* Per-area significance of a user gene list, one row per tested (area, window) cell, most significant first (JSON). score: the list's ssGSEA deviation vs 4,000 expression-matched random gene sets, with a leave-one-donor-out check. markers: Fisher test of the list against each area's specific genes (region classes vs each other, neocortical areas vs the rest of neocortex), with the overlapping genes.
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @param test Which test to run: score, markers or both
#* @serializer json list(digits = NA)
#* @post /area_enrichment
function(gene_string = "", test = "both") {
  assinatura <- assinatura_valida(gene_string)
  if (nrow(assinatura$validos) == 0) {
    stop("Error: none of the submitted genes was found in the expression matrix.")
  }
  tests <- if (identical(tolower(test), "both")) c("score", "markers") else parse_significance(test)
  if (identical(tests, "none")) stop("Error: test must be score, markers or both.")

  out <- list(weighted = jsonlite::unbox(assinatura$ponderado))
  for (t in tests) {
    sig <- run_significance(assinatura, t)
    n_tested <- if (is.null(sig$cells)) 0L else sum(sig$cells$tested)
    out[[t]] <- list(
      reason     = jsonlite::unbox(if (is.null(sig$reason)) NA_character_ else sig$reason),
      n_tested   = jsonlite::unbox(n_tested),
      n_untested = jsonlite::unbox(if (is.null(sig$cells)) 0L else nrow(sig$cells) - n_tested),
      n_enriched = jsonlite::unbox(if (is.null(sig$cells)) 0L else
        sum(if (t == "score") sig$cells$enriched else sig$cells$enriched_marker)),
      cells      = significance_rows(sig)
    )
  }
  out
}

#* Developmental stage of a user gene list, one row per (region, window) (JSON). The list's ssGSEA score per window, for the whole brain (region classes equally weighted) and within each PsychENCODE region (NCX, HIP, AMY, STR, MD, CBC), centred over the 5 windows, vs 4,000 expression-matched random gene sets; two-sided, BH over the tested windows (whole brain and regions as separate families), with a leave-one-donor-out check. Region x window cells with fewer than 3 donors are reported but not tested. A log2FC list is also split into up- and down-regulated genes, tested as a second family.
#* @param gene_string One gene per line, optionally followed by its log2FC (weighted signature)
#* @serializer json list(digits = NA)
#* @post /stage_enrichment
function(gene_string = "") {
  assinatura <- assinatura_valida(gene_string)
  if (nrow(assinatura$validos) == 0) {
    stop("Error: none of the submitted genes was found in the expression matrix.")
  }
  res <- run_stage(assinatura)
  rows <- if (is.null(res$rows)) data.frame() else res$rows %>%
    transmute(family, list, region, level, window, n_donors, tested, n_genes, z = round(z, 2), p = signif(p, 3), q = signif(q, 3),
              p_loo = signif(p_loo, 3), direction, enriched)
  list(weighted = jsonlite::unbox(assinatura$ponderado),
       reason   = jsonlite::unbox(if (is.null(res$reason)) NA_character_ else res$reason),
       skipped  = if (is.null(res$rows)) character(0) else attr(res$rows, "skipped"),
       rows     = rows)
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