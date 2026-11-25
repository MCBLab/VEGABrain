#!/usr/bin/env Rscript
# Stage-level test for VEGABrain: in which developmental window does a gene list
# peak (or dip), beyond what random genes of the same expression level do?
#
# The two area tests (enrichment_stats.R, marker_stats.R) compare structures
# within one window, so a rise that happens in every structure at once is
# invisible to them by design. This test asks the complementary "when" question.
#
#   window value  per-sample score -> cell mean (window x structure) -> region
#                 class mean (6 classes, equal weight) -> window mean. Every class
#                 is sampled in all 5 windows, so the regional make-up of a window
#                 is constant even though the structures standing for a class
#                 change (dorsal thalamus vs mediodorsal nucleus, cerebellum vs
#                 cerebellar cortex, fetal lobes vs adult areas).
#   regions       the same test is also run inside each of the six PsychENCODE
#                 regions (NCX, HIP, AMY, STR, MD, CBC): window value = mean over
#                 that region's structures only. A region x window cell is tested
#                 only with >= STAGE_MIN_DONORS donors, as in the area tests; in
#                 the 3rd trimester only NCX and CBC qualify.
#   statistic     the window value minus the mean over the 5 windows.
#   null          the same statistic for k expression-matched random gene sets.
#   p             two-sided competitive p; BH over the tested windows of one
#                 family, whole brain and regions as separate families.
#   guard         drop each donor of the window in turn; a hit must keep p < alpha
#                 every time. It matters most for the 3rd trimester, where two
#                 donors (24 and 37 pcw) give 32 of 38 samples.
#
# Calibration (k = 1000). Whole brain, 300 lists each: random lists of 15 / 50 /
# 200 genes had >= 1 guarded hit in 3.0 / 3.0 / 4.3% of cases; expression-matched
# swaps of real sets 2.3%; random log2FC-weighted lists 3.0%. Regions, 200 lists
# each: 2.5 / 2.5 / 3.5%, swaps 2.5%; per-cell p <= 0.05 at 0.04-0.05 in both.
# 44% (whole brain) and 46% (regions) of real MSigDB sets have a stage hit
# somewhere -- development is the main axis of brain expression -- so read the
# effect size z, not only the mark.
#
# Needs enrichment_stats.R and marker_stats.R sourced first.

STAGE_WINDOWS <- c("1st trimester", "2nd trimester", "3rd trimester", "Infant", "Adult")
STAGE_MIN_GENES <- 5
# PsychENCODE's six regions (Li et al. 2018). MD and CBC stand for the thalamus
# and cerebellum classes: in early windows BrainSpan sampled dorsal thalamus and
# whole cerebellum instead.
STAGE_REGIONS <- c(NCX = "neocortex", HIP = "hippocampus", AMY = "amygdala",
                   STR = "striatum", MD = "thalamus", CBC = "cerebellum")
# Fewer donors than this in a region x window cell: shown, not tested.
STAGE_MIN_DONORS <- 3

window_short <- function(x) sub(" \\(n = .*", "", x)

# One row per sample, with its window and region class.
stage_meta <- function(col_meta) {
  meta <- col_meta %>%
    filter(!is.na(structure_mapped), !is.na(broad_age)) %>%
    distinct(column_num, donor_id, structure_original, broad_age) %>%
    mutate(window = window_short(broad_age),
           class  = unname(STRUCTURE_CLASS[structure_original]))
  if (anyNA(meta$class)) stop("structure without a region class: ",
                              paste(unique(meta$structure_original[is.na(meta$class)]), collapse = ", "))
  meta
}

# [samples x windows] weights such that S %*% M gives the window values.
stage_weights <- function(meta, sample_ids) {
  meta <- meta %>%
    group_by(window, structure_original) %>% mutate(n_cell = n()) %>%
    group_by(window, class) %>% mutate(n_struct = n_distinct(structure_original)) %>%
    group_by(window) %>% mutate(n_class = n_distinct(class)) %>% ungroup() %>%
    mutate(wt = 1 / (n_cell * n_struct * n_class))
  M <- matrix(0, length(sample_ids), length(STAGE_WINDOWS), dimnames = list(sample_ids, STAGE_WINDOWS))
  i <- match(as.character(meta$column_num), sample_ids)
  ok <- !is.na(i)
  M[cbind(i[ok], match(meta$window[ok], STAGE_WINDOWS))] <- meta$wt[ok]
  M
}

# S: [1 + k, samples], row 1 the query. One row per window, unadjusted.
stage_stat <- function(S, meta) {
  ids  <- colnames(S)
  meta <- meta[as.character(meta$column_num) %in% ids, ]
  V  <- S %*% stage_weights(meta, ids)
  D  <- V - rowMeans(V)
  k  <- nrow(S) - 1
  Dq <- D[1, ]; Dn <- D[-1, , drop = FALSE]
  p_up   <- (1 + colSums(Dn >= rep(Dq, each = k))) / (k + 1)
  p_down <- (1 + colSums(Dn <= rep(Dq, each = k))) / (k + 1)
  donors <- tapply(meta$donor_id, factor(meta$window, STAGE_WINDOWS), function(x) length(unique(x)))
  data.frame(window    = STAGE_WINDOWS,
             n_donors  = as.integer(donors[STAGE_WINDOWS]),
             z         = (Dq - colMeans(Dn)) / apply(Dn, 2, sd),
             p         = pmin(1, 2 * pmin(p_up, p_down)),
             direction = ifelse(Dq >= colMeans(Dn), "higher", "lower"),
             row.names = NULL, stringsAsFactors = FALSE)
}

# Leave one donor out, for rows `rows` of a stage_stat() table. Returns the worst
# two-sided p, in the hit's direction, over the donors of each window.
stage_guard <- function(S, meta, st, rows) {
  ids  <- colnames(S)
  meta <- meta[as.character(meta$column_num) %in% ids, ]
  V <- S %*% stage_weights(meta, ids)
  k <- nrow(S) - 1
  vapply(rows, function(h) {
    w <- st$window[h]
    max(vapply(unique(meta$donor_id[meta$window == w]), function(dn) {
      Md <- stage_weights(meta[!(meta$donor_id == dn & meta$window == w), ], ids)
      Vd <- V; Vd[, h] <- S %*% Md[, h]
      Dd <- Vd[, h] - rowMeans(Vd)
      p1 <- if (st$direction[h] == "higher") (1 + sum(Dd[-1] >= Dd[1])) / (k + 1)
            else (1 + sum(Dd[-1] <= Dd[1])) / (k + 1)
      min(1, 2 * p1)
    }, numeric(1)))
  }, numeric(1))
}

# ssGSEA of a plain list and its k expression-matched null sets.
stage_scores <- function(rt, genes, k = 4000, seed = 1) {
  idx <- match(unique(genes), rt$genes)
  idx <- idx[!is.na(idx)]
  if (length(idx) < STAGE_MIN_GENES) stop(sprintf("needs at least %d genes; got %d", STAGE_MIN_GENES, length(idx)))
  ssgsea_from_ranks(rt, rbind(membership(idx, rt$n), null_sets(rt, idx, k, seed)))
}

# The units a list is tested in: the whole brain (class-balanced) and each region.
stage_units <- function(meta) {
  c(list(brain = meta), lapply(STAGE_REGIONS, function(cl) meta[meta$class == cl, , drop = FALSE]))
}

# One family per level (whole brain; regions): BH over the tested rows of every
# list in Ss, donor guard on BH hits.
# Ss: named list of score matrices, e.g. list(all = S) or list(up = S1, down = S2).
stage_family <- function(Ss, meta, alpha = 0.05) {
  units <- stage_units(meta)
  st <- bind_rows(lapply(names(Ss), function(nm) bind_rows(lapply(names(units), function(u)
    stage_stat(Ss[[nm]], units[[u]]) %>% mutate(list = nm, region = u))))) %>%
    mutate(level = ifelse(region == "brain", "brain", "region"),
           tested = !is.na(n_donors) & n_donors >= STAGE_MIN_DONORS)
  st$q <- NA_real_
  for (lv in unique(st$level)) {
    i <- which(st$level == lv & st$tested)
    st$q[i] <- p.adjust(st$p[i], "BH")
  }
  st$p_loo <- NA_real_
  hit <- which(!is.na(st$q) & st$q < alpha)
  for (key in unique(paste(st$list[hit], st$region[hit]))) {
    rows <- which(paste(st$list, st$region) == key)
    h <- intersect(hit, rows)
    st$p_loo[h] <- stage_guard(Ss[[st$list[rows[1]]]], units[[st$region[rows[1]]]], st[rows, ], match(h, rows))
  }
  st %>% mutate(enriched = !is.na(q) & q < alpha & !is.na(p_loo) & p_loo < alpha)
}

# Top level. The whole list is one family; for a log2FC list, the up- and
# down-regulated genes are scored separately (ssGSEA, like a plain list) and form
# a second family, so the "all" row reads the same with or without log2FC.
stage_enrichment <- function(rt, meta, genes, weights = NULL, k = 4000, alpha = 0.05) {
  out <- stage_family(list(all = stage_scores(rt, genes, k, seed = 1)), meta, alpha) %>%
    mutate(family = "list", n_genes = length(intersect(genes, rt$genes)))
  skipped <- character(0)
  if (!is.null(weights)) {
    sides <- list(up = genes[weights > 0], down = genes[weights < 0])
    Ss <- list(); n <- c()
    for (s in names(sides)) {
      g <- intersect(sides[[s]], rt$genes)
      if (length(g) >= STAGE_MIN_GENES) {
        Ss[[s]] <- stage_scores(rt, g, k, seed = if (s == "up") 1 else 2); n[s] <- length(g)
      } else skipped <- c(skipped, s)
    }
    if (length(Ss)) out <- bind_rows(out, stage_family(Ss, meta, alpha) %>%
                                       mutate(family = "direction", n_genes = unname(n[list])))
  }
  attr(out, "skipped") <- skipped
  out
}
