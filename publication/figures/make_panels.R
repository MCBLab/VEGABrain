#!/usr/bin/env Rscript
# Figures 2-5 of main.tex, drawn with the app's own code at print resolution.
#
# The brain maps come from plumber.R's grafico_*() functions -- the code behind
# the API's *_png endpoints -- so the paper shows exactly what the tool serves,
# only at 300 dpi instead of the endpoints' 150. The developmental-stage panels
# come from stats/stage_stats.R, the code behind POST /stage_enrichment, and the
# marker genes from stats/marker_stats.R (POST /area_enrichment).
#
#   Rscript publication/figures/make_panels.R          # all figures
#   Rscript publication/figures/make_panels.R fig5     # one figure
#
# Run from anywhere inside the repository; needs the prepared data objects in
# its root (dados_otimizados.rds, ontologia_*.rds, genesets_list.rds), built by
# preparaDados.R.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(cowplot) })

find_root <- function() {
  for (p in c(".", "..", "../..", "../../..")) {
    if (file.exists(file.path(p, "dados_otimizados.rds")) && file.exists(file.path(p, "plumber.R"))) {
      return(normalizePath(p))
    }
  }
  stop("run from inside the repository (dados_otimizados.rds and plumber.R not found)")
}
ROOT <- find_root()
OUT  <- file.path(ROOT, "publication", "figures")
only <- commandArgs(trailingOnly = TRUE)
want <- function(id) length(only) == 0 || id %in% only

# plumber.R reads its data relative to the working directory.
setwd(ROOT)
message("loading the app (plumber.R) from ", ROOT)
suppressMessages(suppressWarnings(source("plumber.R", local = globalenv())))

PAL    <- c("#1A318B", "#4F71BE", "#C2B4D6", "#D1498C", "#7A0845")
REGION <- c(brain = "Whole brain", NCX = "NCX", HIP = "HIP", AMY = "AMY", STR = "STR", MD = "MD", CBC = "CBC")
WINDOW <- c("1st trimester" = "1st", "2nd trimester" = "2nd", "3rd trimester" = "3rd",
            Child = "Child", Adult = "Adult")

# ---- developmental-stage panel ------------------------------------------
# One list: rows = whole brain + the six PsychENCODE regions, columns = the
# five windows. Colour = z against random genes; star = q < 0.05 and robust to
# leaving any donor out; ring = q < 0.05 resting on one donor; hatched = fewer
# than 3 donors, not tested.
stage_panel <- function(genes, weights = NULL, title = "Developmental stage") {
  # A log2FC list is drawn as the weighted list the map shows (list "all"); its
  # up and down halves stay in the table only.
  st <- stage_enrichment(rank_tables_cached(), stage_ref, genes, weights)
  attr(st, "skipped") <- NULL
  df <- st %>% filter(list == "all")
  df <- df %>%
    mutate(region = factor(REGION[region], rev(REGION)), window = factor(WINDOW[window], WINDOW),
           zc = pmax(-8, pmin(8, z)),
           mark = case_when(!tested ~ "", enriched ~ "★", !is.na(q) & q < 0.05 ~ "○", TRUE ~ ""),
           ink = ifelse(abs(zc) > 4.5, "white", "grey15"))
  hatch <- df %>% filter(!tested) %>%
    mutate(x = as.numeric(window), y = as.numeric(region)) %>%
    tidyr::crossing(k = c(-0.3, 0, 0.3)) %>%
    mutate(x0 = x - 0.45 + pmax(0, k), y0 = y - 0.45 + pmax(0, -k),
           x1 = x + 0.45 - pmax(0, -k), y1 = y + 0.45 - pmax(0, k))
  p <- ggplot(df, aes(window, region)) +
    geom_tile(aes(fill = zc, alpha = tested), colour = "white", linewidth = 0.7) +
    geom_segment(data = hatch, aes(x = x0, y = y0, xend = x1, yend = y1),
                 colour = "grey55", linewidth = 0.3, inherit.aes = FALSE) +
    geom_text(aes(label = mark, colour = ink), size = 5, show.legend = FALSE) +
    scale_colour_identity() +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    scale_fill_gradientn(colours = c(PAL[1], PAL[2], "grey95", PAL[4], PAL[5]), limits = c(-8, 8),
                         breaks = c(-8, 0, 8), labels = c("lower", "0", "higher"),
                         name = "z vs random genes") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14, hjust = 0),
          axis.text.x = element_text(face = "bold", size = 11.5, colour = "grey20"),
          axis.text.y = element_text(size = 12.5, colour = "grey20"),
          legend.position = "bottom", legend.title.position = "top",
          legend.title = element_text(size = 12, hjust = 0.5),
          legend.text = element_text(size = 11),
          legend.key.width = unit(0.85, "cm"), legend.key.height = unit(0.28, "cm"),
          plot.margin = margin(18, 8, 6, 4),
          plot.background = element_rect(fill = "white", colour = NA))
  attr(p, "table") <- st
  p
}

# Text key under the stage panel, shared by Figs 3-5.
stage_key <- function() {
  ggdraw() +
    draw_label(paste("★     q < 0.05, robust to any one donor",
                     "○     q < 0.05, rests on one donor",
                     "hatched: < 3 donors, not tested", sep = "\n"),
               x = 0.04, hjust = 0, size = 11.5, colour = "grey25", lineheight = 1.15)
}

panel <- function(map, side, side_heights = c(1, 0.2), widths = c(8.2, 3.7)) {
  right <- plot_grid(plotlist = side, ncol = 1, rel_heights = side_heights)
  plot_grid(map, right, nrow = 1, rel_widths = widths, labels = c("a", "b"),
            label_size = 15, label_fontface = "bold")
}

save_fig <- function(p, name, w, h) {
  # The cairo PNG device, as the API's png serializer uses: same fonts and the
  # same asterisk halos as the maps the site serves.
  ggsave(file.path(OUT, name), p, width = w, height = h, dpi = 300, bg = "white",
         device = grDevices::png, type = "cairo")
  message("wrote ", name, sprintf(" (%.1f x %.1f in, 300 dpi)", w, h))
}

# ---- Fig. 2: single gene ---------------------------------------------------
if (want("fig2")) {
  save_fig(grafico_gene("SOX10", "micro"), "fig2_gene_sox10.png", 11, 7)
}

# ---- Fig. 3: ontology term + its stage test --------------------------------
# Astrocyte development peaks in the 3rd trimester everywhere at once, so the two
# area tests (which compare areas within a window) mark nothing there; the stage
# test is the one that sees it.
if (want("fig3")) {
  gs <- "GOBP_ASTROCYTE_DEVELOPMENT"
  sp <- stage_panel(intersect(genesets_list[[gs]], setgene))
  save_fig(panel(grafico_ontologia(gs, "micro"), list(sp, stage_key())), "fig3_ontology.png", 11.9, 6.8)
  write.csv(attr(sp, "table"), file.path(OUT, "genes", "fig3_stage.csv"), row.names = FALSE)
}

# ---- Fig. 4: autism genes (HP_AUTISM): marker test, stage, marker genes ----
if (want("fig4")) {
  aut <- intersect(genesets_list[["HP_AUTISM"]], setgene)
  sp  <- stage_panel(aut)
  map <- grafico_genelist(paste(aut, collapse = "\n"), "micro", "markers")
  mk  <- area_enrichment(marker_reference(), aut) %>% filter(enriched_marker) %>% arrange(q_marker)
  wrap_genes <- function(g, width = 32) paste(strwrap(gsub(",", "", g), width = width), collapse = "\n")
  txt <- paste0("Marker test (", "◆", ")\n\n",
                paste(sprintf("%s, %s\n%d of %d markers, q = %s\n%s",
                              sub(" \\(n = .*", "", mk$broad_age),
                              gsub(" \\([^)]*\\)", "", mk$cell), mk$overlap, mk$n_markers,
                              formatC(mk$q_marker, digits = 2, format = "g"),
                              vapply(mk$overlap_genes, wrap_genes, character(1))),
                      collapse = "\n\n"))
  genes_panel <- ggdraw() + draw_label(txt, x = 0.09, y = 0.97, hjust = 0, vjust = 1, size = 12,
                                       colour = "grey20", lineheight = 1.12)
  right <- plot_grid(sp, stage_key(), genes_panel, ncol = 1, rel_heights = c(1, 0.19, 0.78),
                     labels = c("b", "", "c"), label_size = 15, label_fontface = "bold")
  fig <- plot_grid(map, right, nrow = 1, rel_widths = c(8.2, 3.7), labels = c("a", ""),
                   label_size = 15, label_fontface = "bold")
  save_fig(fig, "fig4_autism.png", 11.9, 9.2)
  write.csv(attr(sp, "table"), file.path(OUT, "genes", "fig4_stage.csv"), row.names = FALSE)
  write.csv(mk %>% transmute(window = broad_age, area = cell, level, overlap, n_markers,
                             odds_ratio, q = q_marker, genes = overlap_genes),
            file.path(OUT, "genes", "fig4_markers.csv"), row.names = FALSE)
}

# ---- Fig. 5: ethanol in cerebral organoids (Arzua et al. 2020), log2FC ------
# genes/fig5_organoid_arzua2020.csv is Supplementary Table 1 of Arzua et al.
# (2020, Transl Psychiatry 10:347): 199 genes changed > 2-fold (P < 0.05) after
# 6 h of 50 mM ethanol in two-month-old iPSC cerebral organoids. Signed log2 fold
# change = log2(fold change), negative for "down". The whole table is submitted,
# as a user would paste it; the app keeps the 173 symbols BrainSpan has.
if (want("fig5")) {
  org <- read.csv(file.path(OUT, "genes", "fig5_organoid_arzua2020.csv"))
  txt <- paste(sprintf("%s %.3f", org$gene, org$log2fc), collapse = "\n")
  sig <- assinatura_valida(txt)
  sp  <- stage_panel(sig$validos$gene, sig$validos$peso)
  map <- grafico_genelist(txt, "micro", "score")
  save_fig(panel(map, list(sp, stage_key())), "fig5_organoid.png", 11.9, 7.7)
  write.csv(attr(sp, "table"), file.path(OUT, "genes", "fig5_stage.csv"), row.names = FALSE)
  score <- run_significance(sig, "score")$cells
  write.csv(score %>% filter(tested) %>%
              transmute(window = broad_age, structure = structure_original, n_donors,
                        delta = delta_within_window, p = p_competitive, q = q_competitive,
                        p_loo, enriched),
            file.path(OUT, "genes", "fig5_score.csv"), row.names = FALSE)
}
