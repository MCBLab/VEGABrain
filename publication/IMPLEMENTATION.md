# VEGABrain manuscript — implementation notes

How the *Neuroinformatics* submission in this folder is put together: what each
file is for, where every number in the text comes from, what is still a
placeholder, and how to build the PDF.

- **Target journal:** *Neuroinformatics* (Springer), <https://link.springer.com/journal/12021>
- **Article type:** **Software Original Article** — the journal's type for a paper
  describing a software product. It asks for more than a description: a short
  account of how to install and run the tool, a validation of its correctness
  against accepted standards, and an evaluation against comparable packages
  including benchmarking where applicable. Section 4 of `main.tex` exists
  specifically to satisfy that requirement (`Kennedy2009` in the bibliography is
  the journal's own statement of the bar for software papers).
- **Working title:** VEGABrain: a Visualizer of Enriched Gene Analysis for the
  developing human brain
- **Status:** full text drafted; all four figures are real (the interface
  screenshot was dropped on 2026-09-14), and most measured values in Section 4 are still marked
  `[TODO]`. See [Second pass](#second-pass).

---

## 1. Files

```
publication/
├── main.tex                  the manuscript
├── main.pdf                  built output (19 pp, A4)
├── references.bib            38 entries
├── build.sh                  build / clean
├── sn-jnl.cls                Springer Nature class (copied from template/)
├── sn-basic.bst              bibliography style used by main.tex
├── bst/                      the other eight SN styles, unused
├── template/                 pristine unpacked template, incl. user-manual.pdf
├── figures/
│   ├── export_figures.sh     regenerates fig1-fig4, numbered as in the PDF
│   ├── fig1_architecture.tex Fig. 1 TikZ source -> fig1_architecture.pdf
│   ├── fig{2,3,4}_*.png      Figs. 2-4, the API's own PNG maps
│   └── make_figures.R        offline vector mirror of Figs. 2-4 (pinned ggseg)
├── supplementary/
│   └── S1_region_mapping.csv   Supplementary File 1, the region mapping
├── references/               20 PDFs of cited works
└── IMPLEMENTATION.md         this file
```

The template package is the official Springer Nature journal article template,
v3.1 (December 2024), from
<https://www.springernature.com/gp/authors/campaigns/latex-author-support>.
`template/` holds the untouched extraction, including `user-manual.pdf` (the
authoritative reference for the class) and `sn-article.pdf` (a rendered
example). `sn-jnl.cls` and `sn-basic.bst` are copied up one level because BibTeX
does not search subdirectories.

**Document class:** `\documentclass[pdflatex,sn-basic]{sn-jnl}`. SN Basic is the
author–year reference style *Neuroinformatics* uses. Add the `lineno` option for
the reviewers' copy:

```latex
\documentclass[lineno,pdflatex,sn-basic]{sn-jnl}
```

`\usepackage[T1]{fontenc}` is deliberately **absent**. This TinyTeX has no Type 1
EC fonts, so T1 triggers Metafont bitmap generation, which bloats the PDF and
breaks `latexmk`'s error detection. The SN template leaves fontenc off for the
same reason; UTF-8 accents in author names compose correctly without it.

## 2. Building

```bash
bash build.sh          # -> main.pdf, keeps main.log / main.blg
bash build.sh clean    # remove LaTeX intermediates
```

LaTeX comes from TinyTeX in `~/.TinyTeX` (TeX Live 2026), installed with
`Rscript -e 'tinytex::install_tinytex()'`. The build goes through
`tinytex::latexmk()` rather than calling `latexmk` directly because it
auto-installs any LaTeX package the template pulls in.

Verifying a build:

```bash
grep -c 'Citation .* undefined'  main.log   # expect 0
grep -c 'Reference .* undefined' main.log   # expect 0
grep -c '^!'                     main.log   # expect 0
grep -c 'Overfull \\hbox'        main.log   # expect 0
grep -oP "You've used \K\d+"     main.blg   # expect 38, = entries in references.bib
```

Last verified build (2026-09-04): 19 pages, A4, 0 undefined citations, 0
undefined references, 0 errors, 0 overfull boxes, 38/38 bibliography entries
cited.

## 3. Where the numbers come from

Every quantitative claim in `main.tex` was read off the code or the data in this
repository, not estimated. Keep this table in sync if the pipeline changes.

| Claim in the text | Value | Verified from |
|---|---|---|
| Expression source | BrainSpan RNA-Seq, Gencode v10 RPKM, gene-summarised | `ETL.qmd`, `genes_matrix_csv/` |
| Matrix dimensions | 52,376 genes × 524 samples | `wc -l genes_matrix_csv/{rows,columns}_metadata.csv` |
| Transformation | `log2(RPKM + 1)`, duplicate symbols dropped | `preparaDados.R` |
| Developmental windows | 5; n = 5, 10, 5, 14, 8 donors | `preparaDados.R` (`age_mapping`) |
| Atlas regions (micro) | **41** | `n_distinct(mapeamento_regioes.csv$region)`, and `GET /data_brain?gene=SOX10` returns 41 regions × 5 windows = 169 rows. `ontologia_micro.rds` still says 42 — it is stale; see §5 |
| Anatomical groups (macro) | 6 | `mapeamento_regioes.csv`, column `macro_regions` |
| Views per map | 4 (DK lateral, DK medial, aseg sagittal, aseg coronal_1) | `plumber.R:build_brain_grid` |
| Gene sets scored | 18,650 | `length(readRDS("genesets_list.rds"))` |
| — GOBP / HP / GOMF / REACTOME / GOCC / BIOCARTA / KEGG / HALLMARK | 7,538 / 5,793 / 1,872 / 1,839 / 1,080 / 292 / 186 / 50 | prefix table over `names(genesets_list)` |
| Aggregated score rows | 3,188,466 | `nrow(readRDS("ontologia_micro.rds"))` |
| Sample-level scores before aggregation | ≈9.8 M | 18,650 × 524 |
| Cloud Run config | 4 GiB, 2 vCPU, concurrency 4, timeout 900 s, min-instances 0 | `deploy.sh` |
| Cold start | 24.8 s | timed `GET /list_ontologies` against the live API |
| Live front end | <https://mcblab.github.io/VEGABrain/> | `.github/workflows/pages.yml` |
| Live API | <https://brainmap-api-5kofy56taq-rj.a.run.app> | extracted from the deployed JS bundle |

Regenerate the gene-set breakdown with:

```r
gs <- readRDS("genesets_list.rds")
sort(table(sub("_.*", "", names(gs))), decreasing = TRUE)
```

## 4. Section plan

| § | Section | ~words | Content source |
|---|---|---|---|
| — | Abstract | 200 | — |
| 1 | Introduction | 850 | The three gaps: portals are non-anatomical and gene-at-a-time; programmatic tools assume coding and local data; nothing takes a *gene set* to a developmental map. Contributions listed explicitly. |
| 2.1 | Expression data | 150 | `preparaDados.R`, `ETL.qmd` |
| 2.2 | Developmental windows | 120 | `preparaDados.R` |
| 2.3 | Atlas mapping | 300 | `mapeamento_regioes.csv`, `plumber.R` |
| 2.4 | Composite enrichment scores | 450 | `plumber.R`, `preparaDados.R`; Eq. 1 from `Barbie2009` |
| 2.5 | Implementation & deployment | 300 | `plumber.R`, `Dockerfile.api`, `deploy.sh`, `pages.yml` |
| 3 | Results / tool description | 700 | `front/src/componentes/Home.jsx`, `plumber.R` endpoints |
| 4 | Validation & comparison | 700 | **partly TODO** — see below |
| 5 | Discussion & limitations | 700 | Donor n, bulk tissue, mapping approximation, rank-based scores, cold start |
| 6 | Conclusion | 130 | — |
| — | Back matter | 250 | Declarations; several `[TODO]` |

## 5. Figures and tables

| Item | Label | State | Regenerate with |
|---|---|---|---|
| Fig. 1 Architecture | `fig:arch` | **real** — TikZ, `figures/fig1_architecture.tex` | `build.sh` (when the `.tex` is newer) |
| Fig. 2 Single gene (*SOX10*) | `fig:gene` | **real** | `figures/export_figures.sh` |
| Fig. 3 Ontology composite score | `fig:ontology` | **real** | `figures/export_figures.sh` |
| Fig. 4 User gene list (FAS signature) | `fig:genelist` | **real** | `figures/export_figures.sh` |
| Supplementary File 1 | — | **real** — copy of the tool's `mapeamento_regioes.csv` | by hand (see below) |
| Table 1 Data summary | `tab:data` | **real** | — |
| Table 2 API endpoints | `tab:api` | **real** | `plumber.R` |
| Table 3 Response times | `tab:perf` | cold start real, warm `[TODO]` | benchmark script (to write) |
| Table 4 Tool comparison | `tab:compare` | **real** | — |

Figure files are named by their number in the compiled PDF (`fig1_…` to
`fig4_…`). Fig. 1 is TikZ compiled standalone, because Springer takes figures as
separate files. Figs. 2-4 are the API's own `*_png` output (11 × 7 in, 150 dpi),
downloaded by `figures/export_figures.sh`, so the paper shows exactly what the
tool serves. Point it at an API built from the same tree:

```bash
docker compose up -d plumber-api
API_URL=http://127.0.0.1:33857 bash figures/export_figures.sh
```

Supplementary File 1 is deliberately not refreshed by that script: the root
`mapeamento_regioes.csv` is gitignored and can differ between machines, so the
published copy in `supplementary/` is replaced by hand, from the copy the
deployed image was built with.

`figures/make_figures.R` redraws Figs. 2-4 offline as vector PDF in
`figures/offline/` (gitignored), mirroring `plumber.R:build_brain_grid` — the
same four views, facetting, colour ramp and legends. It is the cross-check
behind the table below, and the source of vector versions if a journal asks. It needs `dados_otimizados.rds`, `ontologyssGSEA.csv` and
`mapeamento_regioes.csv` in the repository root.

**It must run against the ggseg versions `Dockerfile.api` pins**, not whatever is
in the user library. From ggseg 2.2 on, `geom_brain()` only joins data handed to
it directly, so data supplied through `ggplot(df)` — which is how both
`plumber.R` and this script pass it — is silently dropped and the run aborts with
`Column not found`. Install the pinned pair into a scratch library and put it
first on the path:

```bash
mkdir -p /tmp/vegabrain-rlib
Rscript -e '.libPaths(c("/tmp/vegabrain-rlib", .libPaths()))
  remotes::install_github("ggsegverse/ggseg.formats", ref = "v0.0.2",
                          lib = "/tmp/vegabrain-rlib", upgrade = "never")
  remotes::install_github("ggsegverse/ggseg", ref = "v2.1.1",
                          lib = "/tmp/vegabrain-rlib", upgrade = "never")'
R_LIBS=/tmp/vegabrain-rlib Rscript figures/make_figures.R fig2 fig3
```

The API image already carries that pinned pair, so running the script inside it
(`docker run … vegabrain-api Rscript make_figures.R`, with `figures/` mounted
under `/api/publication/figures`) skips the install.

The script rebuilds the sample-to-region mapping from `mapeamento_regioes.csv`
rather than reading `col_meta$structure_mapped`, and aggregates Figure 3 from the
sample-level `ontologyssGSEA.csv` rather than from `ontologia_micro.rds`. Both
detours exist because the `.rds` objects in the repository root were built on
2026-08-05, before the last mapping fix, and still carry two wrong parcels: one
named `bankssts`, which no longer matches the DK atlas and so renders grey, and
`corpus callosum`, onto which the three cerebellum samples were being painted.
Re-running `preparaDados.R` would let both detours go, but it recomputes ssGSEA
over all 18,650 gene sets. See [Second pass](#second-pass).

Both figures were checked against the deployed API, which is the point of
generating them this way:

| | cells vs `/data_*` | max \|Δ\| | Pearson *r* |
|---|---|---|---|
| Fig. 2 `/data_brain?gene=SOX10` | 169 / 169 | 4.6 × 10⁻² | 0.99998 |
| Fig. 3 `/data_ontology?geneset=GOBP_FOREBRAIN_GENERATION_OF_NEURONS` | 169 / 169 | 0 | 1.0 |
| Fig. 4 `POST /data_genelist` (the 13 FAS symbols) | 169 / 169 | 0 | 1.0 |

Figures 3 and 4 are exact. Figure 2's residual is not an error in the figure: the map
caps expression at 6 *before* averaging over samples (`grafico_gene`), while the
CSV endpoint (`/data_brain`) applies no cap at all, so the two differ for any
region holding a sample above 6. Section 3 of `main.tex` currently claims the
opposite — see [Second pass](#second-pass).

Figure 4's demonstration list is `FAS_DORSAL` in `make_figures.R`: the 13
dorsal/glutamatergic specification genes that `Fischer2021` reports as
upregulated by chronic intermittent ethanol in human pluripotent stem cells
differentiated over 50 days to first-trimester-equivalent cortical neurons
(*EMX2*, *LHX1/2/5/9*, *OTX2*, *FEZF2*, *NEUROD1/2/6*, *NEUROG2*, *TBR1*,
*EOMES*). Only genes whose direction that paper states in its text are
included; all 13 are present in the BrainSpan matrix. It was chosen over a
marker panel because the exposure window is part of the experimental design, so
the map is checkable rather than merely illustrative — and it checks out: the
region-averaged score falls monotonically across the five windows (0.81, 0.74,
0.58, 0.47, 0.44), peaking in the window the *in vitro* model represents. Those
five numbers are quoted in Section 4.2 of `main.tex`; `fig_genelist()` prints
them on every run. The paper's other two directional sets — the
ethanol-downregulated ventral/GABAergic programme and the upregulated WNT genes
— were rendered and rejected: both are much flatter across development.

## 6. Reference ledger

38 entries in `references.bib`; 23 have the full text archived in
`references/<Key>.pdf`. Bib key ⇄ filename is 1:1 (`Zhu2019.pdf` was renamed to
`Zhu2018.pdf` and `mowinckel2020.pdf` to `Mowinckel2020.pdf` to keep it so).

**Archived locally (23):** `Bahl2017`, `Cao2023`, `Chen2013`, `Chopra2023`,
`GOConsortium2023`, `Gillespie2022`, `Hanzelmann2013`, `Huisman2017`,
`Kanehisa2000`, `Kohler2021`, `Kuleshov2016`, `Lau2008`, `Li2018`,
`Liberzon2011`, `Liberzon2015`, `Markello2021`, `Markello2022`,
`Mowinckel2020`, `Subramanian2005`, `Sunkin2013`, `Wilkinson2016`, `Yi2020`,
`Zhu2018`.

**No local PDF (15)** — cited on verified Crossref metadata; no open-access copy
was retrievable:

| Key | Why no PDF |
|---|---|
| `Kang2011`, `Miller2014`, `Barbie2009`, `Ashburner2000`, `Fischer2021` | PMC author manuscripts; PMC and Europe PMC both refuse automated PDF retrieval |
| `Hawrylycz2012`, `Desikan2006`, `Fischl2002`, `Kennedy2009`, `Fleck2021` | paywalled, no OA copy |
| `BrainSpan2011` | web resource, not an article |
| `RCoreTeam2026`, `Schloerke2024plumber`, `Chang2025shiny`, `Dolgalev2025msigdbr` | software citations |

Bibliographic metadata for every DOI-registered entry came from the Crossref
REST API, so titles, journals, volumes and pages are publisher records rather
than hand-typed. Author lists longer than 12 names were truncated with
`and others` (14 entries; the consortium papers would otherwise run to hundreds
of names).

To add a reference:

```bash
curl -s "https://api.crossref.org/works/<DOI>?mailto=<you>" | python3 -c '...'   # or by hand
# then archive the PDF if open access:
curl -sL -o references/<Key>.pdf "https://europepmc.org/articles/<pmcid lowercase>?pdf=render"
```

The lowercase `pmcid` matters — `?pdf=render` returns HTTP 500 for an uppercase
`PMCID`.

## 7. Second pass

Everything below is deliberately unfinished in the current draft.

### Figures and measurements

- [x] ~~Capture an interface screenshot~~ — dropped from the manuscript
      (2026-09-14).
- [x] Figure files renumbered to match the PDF (`fig1`-`fig4`), Figs. 2-4
      exported from the API, Fig. 3 legend at two decimals, Fig. 4 legend bar
      thinned; region mapping added as Supplementary File 1 (2026-11-19).
- [x] Figs. 2 and 3 generated and verified against the deployed API (2026-09-04).
- [x] Fig. 4 generated from the Fischer et al. (2021) fetal-alcohol signature
      and verified against `POST /data_genelist` (2026-09-04).
- [ ] Re-run `preparaDados.R` so `dados_otimizados.rds`, `ontologia_micro.rds`
      and `ontologia_macro.rds` pick up the current `mapeamento_regioes.csv`
      (`bankssts` → `banks of superior temporal sulcus`, `corpus callosum`
      dropped). Until then the served maps and the committed `.rds` disagree on
      those two parcels, and `make_figures.R` works around it.
- [ ] **`main.tex` says 42 atlas regions in four places (lines 219, 441, 475,
      557); the deployed API serves 41.** The 42nd is `corpus callosum`, dropped
      by the last mapping fix — the three cerebellum samples that fed it are now
      mapped to `Cerebellum` only. The stale `ontologia_micro.rds` is where the
      42 came from. Fix the four occurrences when the `.rds` objects are rebuilt.
- [ ] Reconcile the map and the CSV export, or soften the claim in
      Section 3 ("the CSV contains exactly the values that were painted, so a
      figure produced in VEGABrain and a statistic computed from its export
      cannot disagree"). `grafico_gene` caps expression at 6 before averaging;
      `/data_brain` does not cap, so for genes with any sample above 6 the two
      disagree — 13 of 169 cells for *SOX10*, by up to 0.046.
- [ ] §4.1 correctness: re-run `GSVA::gsva`/`ssgseaParam` on a random sample of
      gene sets, aggregate as `preparaDados.R` does, and compare against the
      served values. Fill in *n*, Pearson *r*, max |Δ|. Do the same for the
      on-demand path by submitting a precomputed set's member genes.
- [ ] §4.2 face validity: replace the qualitative statements with per-window
      means from the CSV exports for the myelination, GABAergic and neurogenesis
      cases.
- [ ] §4.3 / Table 3: script a warm-latency benchmark of the four endpoints.

### Author and submission metadata

- [ ] Confirm author **order**; add ORCIDs (`\orcid{}` inside `\author`).
- [ ] Replace the two `TODO@ufrn.br` addresses.
- [ ] Full postal affiliation (`\orgdiv`, `\street`, `\postcode`) for the MCB Lab.
- [ ] Funding statement and grant numbers.
- [ ] Author contributions in CRediT terms.

### Repository

- [ ] Choose and add a licence to `MCBLab/VEGABrain` — a software paper needs one.
- [ ] Archive the release on Zenodo and cite the DOI in *Code availability*.
- [x] Merge `VEGAbrain` into `main` and repoint `.github/workflows/pages.yml`
      (done 2026-11-19, together with the rename `BrainMap` → `VEGABrain`).

### Before submission

- [ ] Re-read the journal's submission guidelines at
      <https://link.springer.com/journal/12021/submission-guidelines> — this
      draft was structured from the Software Original Article definition and the
      SN template, and the live page should be checked for length limits and
      declaration wording, which could not be retrieved automatically (the page
      is behind a bot challenge).
- [ ] Build the reviewers' copy with the `lineno` option.
- [ ] Springer wants a single `.tex`; figures are uploaded separately.
