#!/usr/bin/env bash
# Regenerates every figure of main.tex into this folder, numbered as in the PDF:
#
#   fig1_architecture.pdf   Fig. 1  data flow (TikZ, fig1_architecture.tex)
#   fig2_gene_sox10.png     Fig. 2  single-gene mode, SOX10
#   fig3_ontology.png       Fig. 3  ontology mode (GOBP_ASTROCYTE_DEVELOPMENT) +
#                                   its developmental-stage test
#   fig4_autism.png         Fig. 4  HP_AUTISM: marker test marks, stage test,
#                                   marker genes
#   fig5_organoid.png       Fig. 5  ethanol in cerebral organoids (Arzua et al.
#                                   2020, log2FC list): score test + stage test
#
# Figs. 2-5 are drawn by make_panels.R, which sources ../../plumber.R and calls
# the same grafico_*() functions the API's *_png endpoints call, plus the test
# code in ../../stats/ -- so the paper shows exactly what the tool serves, at
# 300 dpi. The numbers behind the panels are written to genes/fig*_stage.csv,
# genes/fig4_markers.csv and genes/fig5_score.csv; Fig. 5's input list is
# genes/fig5_organoid_arzua2020.csv.
#
#   bash publication/figures/export_figures.sh
#
# Supplementary File 1 (../supplementary/S1_region_mapping.csv) is not
# touched here on purpose: the root mapeamento_regioes.csv is gitignored and
# can differ between machines, so the published copy is only replaced by hand.
set -euo pipefail
cd "$(dirname "$0")"

export PATH="$HOME/.TinyTeX/bin/x86_64-linux:$PATH"
Rscript -e 'tinytex::latexmk("fig1_architecture.tex", engine = "pdflatex", clean = TRUE)'
echo "wrote fig1_architecture.pdf"

Rscript make_panels.R
