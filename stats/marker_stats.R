#!/usr/bin/env Rscript
# Marker-based over-representation test per brain area (VEGABrain).
#
# enrichment_stats.R asks whether a list's aggregate ssGSEA score peaks in an
# area. This file asks a different question: are the genes that are *specific*
# to an area, in a developmental window, over-represented in the list -- and
# which ones? It is the Specific Expression Analysis idea (TSEA / pSI, Xu et al.
# 2014) rebuilt on BrainSpan with donor-aware markers:
#
#   markers  limma per window, ~ 0 + structure + donor, so every contrast is a
#            within-donor comparison. A gene is a marker of a cell if BH q < 0.05,
#            log2FC >= 1 and the cell's mean expression is >= 1 log2(RPKM + 1).
#   test     one-sided Fisher's exact test of the list against each cell's
#            markers, inside that window's background; BH over every tested cell
#            of the request (both directions together for a log2FC signature).
#
# The contrast is hierarchical, because neocortex is 15 of the 22 structures:
#   level "class"  a region class (neocortex, hippocampus, amygdala, striatum,
#                  thalamus, cerebellum) vs the equal-weighted mean of the other
#                  classes present in the window; structures are equal-weighted
#                  within their class.
#   level "area"   one neocortical area vs the mean of the other neocortical areas.
# Comparing each structure with the mean of all the others lets cortex dominate
# the reference (the bias that broke the first design described in
# enrichment_stats.R), and requiring a gene to beat every other structure leaves
# no cortical area with a single marker. Few area-level markers after birth is
# biology (cortical areas converge postnatally), not a defect; every cell reports
# n_markers, and cells with too few markers are not tested.
#
# Pure functions. The only inputs are the expression matrix, its sample metadata
# and the gene universe; nothing is precomputed or read from disk.

suppressPackageStartupMessages({ library(dplyr) })

MARKER_Q          <- 0.05  # BH q of the limma contrast, per cell
MARKER_LFC        <- 1     # log2 fold change vs the reference (2x)
MARKER_MIN_EXPR   <- 1     # cell mean, log2(RPKM + 1)
MIN_MARKERS       <- 10    # a cell with fewer markers is not tested
MARKER_MIN_DONORS <- 3     # informative donors (>= 2 samples in the window) in the cell
MIN_REF_SAMPLES   <- 2     # a structure needs this many samples to enter a window's design
MIN_LIST_GENES    <- 5     # per direction, genes present in the expression matrix
MIN_OVERLAP       <- 2     # an enriched cell must share at least this many genes with the list

# Every BrainSpan structure the atlas can draw, and its region class. Explicit on
# purpose: a structure missing here is an error, not a silent drop.
STRUCTURE_CLASS <- c(
  "anterior (rostral) cingulate (medial prefrontal) cortex" = "neocortex",
  "dorsolateral prefrontal cortex"                          = "neocortex",
  "orbital frontal cortex"                                  = "neocortex",
  "ventrolateral prefrontal cortex"                         = "neocortex",
  "primary motor cortex (area M1, area 4)"                  = "neocortex",
  "primary motor-sensory cortex (samples)"                  = "neocortex",
  "primary somatosensory cortex (area S1, areas 3,1,2)"     = "neocortex",
  "parietal neocortex"                                      = "neocortex",
  "posteroventral (inferior) parietal cortex"               = "neocortex",
  "posterior (caudal) superior temporal cortex (area 22c)"  = "neocortex",
  "temporal neocortex"                                      = "neocortex",
  "inferolateral temporal cortex (area TEv, area 20)"       = "neocortex",
  "primary auditory cortex (core)"                          = "neocortex",
  "primary visual cortex (striate cortex, area V1/17)"      = "neocortex",
  "occipital neocortex"                                     = "neocortex",
  "hippocampus (hippocampal formation)"                     = "hippocampus",
  "amygdaloid complex"                                      = "amygdala",
  "striatum"                                                = "striatum",
  "dorsal thalamus"                                         = "thalamus",
  "mediodorsal nucleus of thalamus"                         = "thalamus",
  "cerebellum"                                              = "cerebellum",
  "cerebellar cortex"                                       = "cerebellum")


# ---- samples --------------------------------------------------------------
# One row per sample. col_meta repeats a sample once per atlas parcel it paints,
# so it is collapsed here; structure_mapped restricts to drawable structures.
marker_samples <- function(col_meta) {
  meta <- col_meta %>%
    filter(!is.na(structure_mapped), !is.na(broad_age)) %>%
    distinct(column_num, donor_id, structure_original, broad_age) %>%
    as.data.frame()
  unknown <- setdiff(unique(meta$structure_original), names(STRUCTURE_CLASS))
  if (length(unknown)) stop("structure without a region class: ", paste(unknown, collapse = ", "))
  meta$class <- unname(STRUCTURE_CLASS[meta$structure_original])
  meta
}

# Samples that can inform a within-donor contrast: a donor with a single sample
# in the window is absorbed by its own coefficient, and a structure with fewer
# than MIN_REF_SAMPLES samples would only add noise to the references. Dropping
# one can expose the other, so repeat until stable.
design_samples <- function(mw) {
  repeat {
    n0 <- nrow(mw)
    mw <- mw[mw$donor_id %in% names(which(table(mw$donor_id) >= 2)), ]
    mw <- mw[mw$structure_original %in% names(which(table(mw$structure_original) >= MIN_REF_SAMPLES)), ]
    if (nrow(mw) == n0) return(mw)
  }
}


# ---- contrasts ------------------------------------------------------------
# Rows: structure coefficients (then donor coefficients, all zero). Columns: one
# per cell. `cls` is the class of each structure, in coefficient order.
# hierarchical: classes vs the mean of the other classes, neocortical areas vs the
#               mean of the other neocortical areas.
# flat:         every structure vs the mean of all other structures (validation
#               only: shows how cortex dominates that reference).
build_contrasts <- function(structures, cls, contrast = c("hierarchical", "flat")) {
  contrast <- match.arg(contrast)
  ns <- length(structures)
  cols <- list(); cells <- list()
  add <- function(v, level, cell, members) {
    cols[[length(cols) + 1]] <<- v
    cells[[length(cells) + 1]] <<- data.frame(level = level, cell = cell,
                                              structures = paste(members, collapse = "; "),
                                              stringsAsFactors = FALSE)
  }
  if (contrast == "hierarchical") {
    classes <- unique(cls)
    if (length(classes) >= 2) {
      for (k in classes) {
        v <- numeric(ns)
        v[cls == k] <- 1 / sum(cls == k)
        for (o in setdiff(classes, k)) v[cls == o] <- -1 / sum(cls == o) / (length(classes) - 1)
        add(v, "class", k, structures[cls == k])
      }
    }
    ncx <- which(cls == "neocortex")
    if (length(ncx) >= 2) {
      for (a in ncx) {
        v <- numeric(ns)
        v[a] <- 1
        v[setdiff(ncx, a)] <- -1 / (length(ncx) - 1)
        add(v, "area", structures[a], structures[a])
      }
    }
  } else {
    for (a in seq_len(ns)) {
      v <- rep(-1 / (ns - 1), ns); v[a] <- 1
      add(v, "structure", structures[a], structures[a])
    }
  }
  list(C = do.call(cbind, cols), cells = bind_rows(cells))
}


# ---- markers --------------------------------------------------------------
# m: expression matrix [genes x samples], log2(RPKM + 1), colnames = column_num.
# universe: genes a marker (and a list gene) may be; NULL = all rows.
# Returns one element per window:
#   $background  genes of the universe expressed (>= min_expr) in some structure
#   $cells       one row per cell: window, level, cell, structures, n_donors,
#                n_markers, testable
#   $logfc, $q, $cell_mean  [background x cells], to re-threshold without refitting
#   $markers     named list, cell_id -> marker genes
# donor = "block" swaps the donor fixed effect for limma's duplicateCorrelation
# (a random donor effect); ~10x slower, used only to validate the default.
area_markers <- function(m, col_meta, universe = NULL, lfc = MARKER_LFC, q = MARKER_Q,
                         min_expr = MARKER_MIN_EXPR, contrast = "hierarchical",
                         donor = c("fixed", "block")) {
  donor <- match.arg(donor)
  if (!requireNamespace("limma", quietly = TRUE)) stop("the marker test needs the limma package")
  meta <- marker_samples(col_meta)
  genes <- if (is.null(universe)) rownames(m) else rownames(m)[rownames(m) %in% universe]

  windows <- sort(unique(meta$broad_age))
  out <- lapply(windows, function(w) {
    mw <- design_samples(meta[meta$broad_age == w, ])
    structures <- sort(unique(mw$structure_original))
    cls <- unname(STRUCTURE_CLASS[structures])
    X <- m[genes, as.character(mw$column_num), drop = FALSE]

    # mean per structure, then per cell (equal weight per structure)
    S <- sapply(structures, function(s) rowMeans(X[, mw$structure_original == s, drop = FALSE]))
    expressed <- apply(S, 1, max) >= min_expr
    X <- X[expressed, , drop = FALSE]
    S <- S[expressed, , drop = FALSE]

    # short coefficient names; the structure names are too long for limma
    code <- setNames(paste0("s", seq_along(structures)), structures)
    sf <- factor(code[mw$structure_original], levels = code)
    design <- if (donor == "fixed") model.matrix(~ 0 + sf + factor(mw$donor_id)) else model.matrix(~ 0 + sf)
    if (qr(design)$rank < ncol(design)) stop("design not full rank in ", w)

    ct <- build_contrasts(structures, cls, contrast)
    C <- rbind(ct$C, matrix(0, ncol(design) - length(structures), ncol(ct$C)))
    rownames(C) <- colnames(design)
    lm <- if (donor == "fixed") limma::lmFit(X, design) else {
      corr <- limma::duplicateCorrelation(X, design, block = mw$donor_id)$consensus.correlation
      limma::lmFit(X, design, block = mw$donor_id, correlation = corr)
    }
    fit <- limma::eBayes(limma::contrasts.fit(lm, C), trend = TRUE)

    cell_mean <- S %*% pmax(ct$C, 0)  # positive weights sum to 1: the cell's own mean
    qv <- apply(fit$p.value, 2, p.adjust, method = "BH")
    cell_id <- paste(w, ct$cells$level, ct$cells$cell, sep = " || ")
    colnames(qv) <- colnames(cell_mean) <- cell_id
    logfc <- fit$coefficients; colnames(logfc) <- cell_id

    members <- strsplit(ct$cells$structures, "; ", fixed = TRUE)
    n_donors <- vapply(members, function(s) n_distinct(mw$donor_id[mw$structure_original %in% s]), integer(1))

    res <- list(window = w, background = rownames(X), logfc = logfc, q = qv,
                cell_mean = cell_mean,
                cells = data.frame(cell_id = cell_id, broad_age = w, ct$cells,
                                   n_donors = n_donors, stringsAsFactors = FALSE))
    select_markers(res, lfc = lfc, q = q, min_expr = min_expr)
  })
  names(out) <- windows
  out
}

# (Re)apply the marker thresholds to one window of area_markers(), filling
# $markers, $cells$n_markers and $cells$testable.
select_markers <- function(win, lfc = MARKER_LFC, q = MARKER_Q, min_expr = MARKER_MIN_EXPR) {
  is_marker <- win$q < q & win$logfc >= lfc & win$cell_mean >= min_expr
  win$markers <- lapply(setNames(seq_len(ncol(is_marker)), colnames(is_marker)),
                        function(j) win$background[is_marker[, j]])
  win$cells$n_markers <- lengths(win$markers)
  win$cells$testable  <- win$cells$n_markers >= MIN_MARKERS &
                         win$cells$n_donors >= MARKER_MIN_DONORS
  win
}


# ---- the test -------------------------------------------------------------
# One row per cell of every window. p_marker is the one-sided Fisher p (computed
# as the hypergeometric tail, which is identical for a 2x2 table); details = TRUE
# adds the conditional-MLE odds ratio and its one-sided lower bound, which only a
# reported result needs.
ora_cells <- function(mk, genes, direction = "all", details = TRUE) {
  genes <- unique(genes)
  bind_rows(lapply(mk, function(win) {
    L  <- genes[genes %in% win$background]
    nb <- length(win$background); nl <- length(L)
    ov <- lapply(win$markers, function(M) L[L %in% M])
    k  <- lengths(ov); nm <- win$cells$n_markers
    out <- win$cells %>%
      select(cell_id, broad_age, level, cell, structures, n_donors, n_markers) %>%
      mutate(direction = direction, n_list = nl, n_background = nb, overlap = k,
             expected = nl * nm / nb,
             p_marker = phyper(k - 1, nm, nb - nm, nl, lower.tail = FALSE),
             tested = win$cells$testable & nl > 0)
    if (details) {
      ft <- lapply(seq_along(k), function(i) {
        if (k[i] == 0) return(c(0, 0))
        f <- fisher.test(matrix(c(k[i], nl - k[i], nm[i] - k[i], nb - nm[i] - nl + k[i]), 2),
                         alternative = "greater")
        c(unname(f$estimate), f$conf.int[1])
      })
      out$odds_ratio <- vapply(ft, `[`, numeric(1), 1)
      out$or_ci_low  <- vapply(ft, `[`, numeric(1), 2)
      out$overlap_genes <- vapply(ov, function(g) paste(sort(g), collapse = ", "), character(1))
    }
    out
  }))
}

# genes: symbols present in the matrix. weights: NULL for a plain list, or the
# log2FC of each gene (lines without a value count as +1, as in score_ponderado),
# in which case up and down genes are tested separately. BH runs over every
# tested cell of the request, both directions together: one error budget per
# figure. Directions with fewer than MIN_LIST_GENES genes are skipped and named in
# attr(, "skipped").
area_enrichment <- function(mk, genes, weights = NULL, alpha = 0.05, details = TRUE) {
  sets <- if (is.null(weights)) list(all = genes) else
    list(up = genes[weights > 0], down = genes[weights < 0])
  sizes <- lengths(sets)
  skipped <- names(sets)[sizes < MIN_LIST_GENES]
  sets <- sets[sizes >= MIN_LIST_GENES]
  if (length(sets) == 0) {
    stop(if (is.null(weights)) {
      sprintf("needs at least %d genes found in the data (got %d)", MIN_LIST_GENES, sizes[["all"]])
    } else {
      sprintf("needs at least %d up- or %d down-regulated genes found in the data (got %d up, %d down)",
              MIN_LIST_GENES, MIN_LIST_GENES, sizes[["up"]], sizes[["down"]])
    })
  }
  res <- bind_rows(lapply(names(sets), function(d) ora_cells(mk, sets[[d]], d, details)))
  res$q_marker <- NA_real_
  res$q_marker[res$tested] <- p.adjust(res$p_marker[res$tested], "BH")
  res$enriched_marker <- res$tested & !is.na(res$q_marker) & res$q_marker < alpha &
                         res$overlap >= MIN_OVERLAP
  attr(res, "skipped") <- skipped
  res
}


# ---- painting outward -----------------------------------------------------
# Cells are classes or neocortical areas; the map draws parcels. A parcel is
# enriched in a window if any tested cell containing one of its structures is (a
# neocortex class hit therefore marks every cortical parcel of that window). q is
# the smallest over those cells; level says which kind of cell carried the call.
paint_markers <- function(res, mapping) {
  if (!"overlap_genes" %in% names(res)) res$overlap_genes <- ""   # details = FALSE
  res %>%
    filter(tested) %>%
    mutate(structure_original = strsplit(structures, "; ", fixed = TRUE)) %>%
    tidyr::unnest(structure_original) %>%
    inner_join(mapping %>% select(structure_name, region),
               by = c("structure_original" = "structure_name"),
               relationship = "many-to-many") %>%
    group_by(broad_age, region, direction) %>%
    summarise(enriched = any(enriched_marker),
              q_marker = min(q_marker),
              level    = paste(sort(unique(level[enriched_marker])), collapse = "+"),
              # the list genes behind the enriched cell(s); empty if none
              genes    = paste(sort(unique(unlist(strsplit(overlap_genes[enriched_marker], ", ", fixed = TRUE)))),
                               collapse = ", "),
              .groups = "drop")
}
