#!/usr/bin/env Rscript
# Calibration, power and controls for the marker-based area test.
#
#   Rscript stats/validate_markers.R [n_null] [n_sets]
#
# Run from the repository root. Needs dados_otimizados.rds and genesets_list.rds
# there (preparaDados.R). The numbers it prints are what README and the API
# caption rely on: a random list must rarely produce a hit, curated anatomical
# sets must find their own structure, and the result must not hinge on how the
# donor effect is modelled. About 10 min at the defaults.

suppressPackageStartupMessages({ library(dplyr) })

find_root <- function() {
  for (p in c(".", "..", "../..")) {
    if (file.exists(file.path(p, "dados_otimizados.rds"))) return(normalizePath(p))
  }
  stop("dados_otimizados.rds not found; run preparaDados.R first")
}
ROOT <- find_root()
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "enrichment_stats.R"))
source(file.path(HERE, "marker_stats.R"))

args   <- commandArgs(trailingOnly = TRUE)
N_NULL <- if (length(args) >= 1) as.integer(args[1]) else 500
N_SETS <- if (length(args) >= 2) as.integer(args[2]) else 300
set.seed(1)

dados_app <- readRDS(file.path(ROOT, "dados_otimizados.rds"))
m   <- dados_app$expression_matrix
gs  <- readRDS(file.path(ROOT, "genesets_list.rds"))
universe <- intersect(unique(unlist(gs, use.names = FALSE)), rownames(m))
message("MSigDB universe: ", length(universe), " symbols")

all_cells <- function(mk) bind_rows(lapply(mk, `[[`, "cells"))
short <- function(x, n = 28) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "~"), x)
win   <- function(x) sub(" \\(n = .*", "", x)

# ----------------------------------------------- 0. design and tested cells
message("\n==> 0. markers (hierarchical contrast, donor fixed effect)")
t0 <- Sys.time()
mk <- area_markers(m, dados_app$col_meta, universe)
cat(sprintf("   built in %.1f s; background per window: %s\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            paste(sapply(mk, function(w) length(w$background)), collapse = " / ")))
cl <- all_cells(mk)
cat(sprintf("   cells: %d, testable (>= %d markers, >= %d informative donors): %d (%d class, %d area)\n",
            nrow(cl), MIN_MARKERS, MARKER_MIN_DONORS, sum(cl$testable),
            sum(cl$testable & cl$level == "class"), sum(cl$testable & cl$level == "area")))
print(cl %>% filter(testable) %>%
        transmute(window = win(broad_age), level, cell = short(cell), n_donors, n_markers) %>%
        as.data.frame(), row.names = FALSE)

# ------------------------------------- 1. why the contrast is hierarchical
# Flat: every structure vs the mean of all the others. Neocortex is most of the
# reference, so cortical areas lose their markers and the rest gains.
message("\n==> 1. marker counts: flat contrast (structure vs all others) vs hierarchical")
flat <- all_cells(area_markers(m, dados_app$col_meta, universe, contrast = "flat"))
flat$class <- unname(STRUCTURE_CLASS[flat$cell])
print(flat %>% group_by(window = win(broad_age), kind = ifelse(class == "neocortex", "cortical area", "non-cortical")) %>%
        summarise(median_markers = median(n_markers), max = max(n_markers), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("   hierarchical, same summary:\n")
print(cl %>% group_by(window = win(broad_age), level) %>%
        summarise(median_markers = median(n_markers), max = max(n_markers), testable = sum(testable),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)

# ---------------------------------------------- 2. stringency sensitivity
message("\n==> 2. testable cells and median markers at log2FC 0.5 / 1 / 2")
for (l in c(0.5, 1, 2)) {
  c2 <- bind_rows(lapply(mk, function(w) select_markers(w, lfc = l)$cells))
  cat(sprintf("   log2FC >= %.1f: testable %2d (class %2d, area %2d); median markers class %4.0f, area %4.0f\n",
              l, sum(c2$testable), sum(c2$testable & c2$level == "class"),
              sum(c2$testable & c2$level == "area"),
              median(c2$n_markers[c2$level == "class"]), median(c2$n_markers[c2$level == "area"])))
}

# ------------------------------------------- 3. donor model: fixed vs block
message("\n==> 3. donor as fixed effect vs duplicateCorrelation (random effect)")
mk_b <- area_markers(m, dados_app$col_meta, universe, donor = "block")
jac <- bind_rows(lapply(names(mk), function(w) {
  a <- mk[[w]]; b <- mk_b[[w]]
  keep <- a$cells$testable | b$cells$testable
  data.frame(window = win(w), cell = short(a$cells$cell[keep]),
             fixed = lengths(a$markers[keep]), block = lengths(b$markers[keep]),
             shared = mapply(function(x, y) length(intersect(x, y)), a$markers[keep], b$markers[keep]))
})) %>% mutate(jaccard = round(shared / pmax(1, fixed + block - shared), 2))
cat(sprintf("   Jaccard of marker sets over testable cells: median %.2f (expected >= 0.9), min %.2f\n",
            median(jac$jaccard), min(jac$jaccard)))
cat("   cells below 0.9 (in a small set one gene moves the index a lot):\n")
print(jac %>% filter(jaccard < 0.9) %>% as.data.frame(), row.names = FALSE)

# ---------------------------------------------- 4. calibration, random lists
# A list drawn uniformly from the universe is the null. Fisher's test is discrete
# and conservative, so per-cell rates sit below nominal; what matters is that a
# null list rarely draws a single mark (family rate <= 0.05).
message("\n==> 4. calibration: ", N_NULL, " uniform random lists per size")
for (n in c(15, 50, 200)) {
  cellp <- c(); fam <- c()
  for (b in seq_len(N_NULL)) {
    r <- area_enrichment(mk, sample(universe, n), details = FALSE)
    cellp <- c(cellp, r$p_marker[r$tested]); fam <- c(fam, any(r$enriched_marker))
  }
  cat(sprintf("   size %3d  per cell: p<=0.01 %.4f  p<=0.05 %.4f  |  lists with >= 1 enriched cell %.3f\n",
              n, mean(cellp <= 0.01), mean(cellp <= 0.05), mean(fam)))
}

# ----------------------------- 5. calibration, expression-matched real sets
# Each gene of a real set swapped for one of the same mean-expression bin keeps
# the set's expression profile and destroys its biology. If markers were biased
# toward an expression range, this is where false hits would appear.
message("\n==> 5. calibration: ", N_SETS, " real sets with every gene swapped for an expression-matched one")
rt <- rank_tables(m, universe = universe)
eff  <- vapply(gs, function(g) length(intersect(unique(g), rownames(m))), integer(1))
pick <- sample(names(gs)[eff >= 15 & eff <= 300], N_SETS)
fam <- c(); cellp <- c()
for (i in seq_along(pick)) {
  idx <- match(intersect(unique(gs[[pick[i]]]), rownames(m)), rt$genes)
  sw  <- rt$genes[draw_matched(rt$bins, idx, k = 1, seed = i, pool = rt$pool)[1, ]]
  r <- area_enrichment(mk, sw, details = FALSE)
  cellp <- c(cellp, r$p_marker[r$tested]); fam <- c(fam, any(r$enriched_marker))
}
cat(sprintf("   per cell p<=0.05 %.4f  |  sets with >= 1 enriched cell %.3f (expected <= 0.05)\n",
            mean(cellp <= 0.05), mean(fam)))

# ------------------------------------------------------ 6. weighted lists
message("\n==> 6. weighted lists: ", N_NULL, " random lists of 50 genes, log2FC ~ N(0, 1)")
cellp <- c(); fam <- c()
for (b in seq_len(N_NULL)) {
  r <- area_enrichment(mk, sample(universe, 50), weights = rnorm(50), details = FALSE)
  cellp <- c(cellp, r$p_marker[r$tested]); fam <- c(fam, any(r$enriched_marker))
}
cat(sprintf("   per cell p<=0.05 %.4f  |  lists with >= 1 enriched cell %.3f\n",
            mean(cellp <= 0.05), mean(fam)))

# --------------------------------------------------------- 7. controls
show_hits <- function(label, r) {
  h <- r %>% filter(enriched_marker)
  cat(sprintf("   %-44s %2d cells  %s\n", short(label, 44), nrow(h),
              if (nrow(h)) paste0(win(h$broad_age), ":", short(h$cell, 18),
                                  ifelse(h$direction == "all", "", paste0("(", h$direction, ")")),
                                  collapse = ", ") else "-"))
}
message("\n==> 7. controls")
show_hits("8 striatal genes (PPP1R1B DRD1 DRD2 ...)",
          area_enrichment(mk, c("PPP1R1B", "DRD1", "DRD2", "PENK", "TAC1", "GPR88", "RGS9", "ADORA2A")))
for (nm in c("GOBP_STRIATUM_DEVELOPMENT", "GOBP_CEREBELLAR_CORTEX_DEVELOPMENT",
             "GOBP_THALAMUS_DEVELOPMENT", "GOBP_HIPPOCAMPUS_DEVELOPMENT",
             "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION", "GOBP_MYELIN_ASSEMBLY",
             "HALLMARK_OXIDATIVE_PHOSPHORYLATION")) {
  if (is.null(gs[[nm]])) next
  show_hits(nm, area_enrichment(mk, intersect(gs[[nm]], rownames(m))))
}

# ------------------------------------------------------------ 8. power
message("\n==> 8. power on the same ", N_SETS, " real MSigDB sets (effective size 15-300)")
res <- bind_rows(lapply(pick, function(nm)
  area_enrichment(mk, intersect(gs[[nm]], rownames(m)), details = FALSE) %>% mutate(set = nm)))
hits <- res %>% filter(enriched_marker)
cat(sprintf("   sets with >= 1 enriched cell: %d / %d (%.0f%%); median cells per hit set %s\n",
            n_distinct(hits$set), N_SETS, 100 * n_distinct(hits$set) / N_SETS,
            if (nrow(hits)) median(table(hits$set)) else "-"))
print(hits %>% count(level, cell, sort = TRUE) %>% mutate(cell = short(cell, 40)) %>% as.data.frame(),
      row.names = FALSE)

# --------------------------------------------- 9. agreement with the score test
# Per (set, window, region class): a score-test hit on any structure of the class
# vs a marker hit on the class (or, for neocortex, on any area).
n_agree <- min(100, N_SETS)
message("\n==> 9. agreement with the score test on ", n_agree, " of those sets")
ref <- cells_from_meta(dados_app$col_meta)
score_hits <- bind_rows(lapply(pick[seq_len(n_agree)], function(nm)
  enrich(ref, rt, gs[[nm]]) %>% filter(enriched) %>%
    transmute(set = nm, broad_age, class = unname(STRUCTURE_CLASS[structure_original])))) %>%
  distinct()
marker_hits <- hits %>% filter(set %in% pick[seq_len(n_agree)]) %>%
  transmute(set, broad_age, class = ifelse(level == "area", "neocortex", cell)) %>% distinct()
grid <- expand.grid(set = pick[seq_len(n_agree)], broad_age = names(mk),
                    class = unique(STRUCTURE_CLASS), stringsAsFactors = FALSE)
key <- function(d) paste(d$set, d$broad_age, d$class)
grid$score  <- key(grid) %in% key(score_hits)
grid$marker <- key(grid) %in% key(marker_hits)
tab <- table(score = grid$score, marker = grid$marker)
print(tab)
cat(sprintf("   both %d, score only %d, markers only %d; agreement odds ratio %.1f\n",
            tab["TRUE", "TRUE"], tab["TRUE", "FALSE"], tab["FALSE", "TRUE"],
            fisher.test(tab)$estimate))
message("\ndone.")
