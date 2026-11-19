#!/usr/bin/env bash
# Regenerates every figure of main.tex into this folder, numbered as in the PDF:
#
#   fig1_architecture.pdf   Fig. 1  data flow (TikZ, fig1_architecture.tex)
#   fig2_gene_sox10.png     Fig. 2  single-gene mode, SOX10
#   fig3_ontology.png       Fig. 3  ontology mode, GOBP_FOREBRAIN_GENERATION_OF_NEURONS
#   fig4_genelist.png       Fig. 4  gene-list mode, FAS dorsal/glutamatergic signature
#
# Figs. 2-4 are downloaded from the API's own *_png endpoints, so the paper
# shows exactly what the tool serves. Point API_URL at an API built from this
# tree -- the deployed one only matches after plumber.R changes are deployed:
#
#   docker compose up -d plumber-api
#   API_URL=http://127.0.0.1:33857 bash figures/export_figures.sh
#
# Supplementary File 1 (../supplementary/S1_region_mapping.csv) is not
# touched here on purpose: the root mapeamento_regioes.csv is gitignored and
# can differ between machines, so the published copy is only replaced by hand.
set -euo pipefail
cd "$(dirname "$0")"

API="${API_URL:?set API_URL to the plumber API, e.g. API_URL=http://127.0.0.1:33857}"
API="${API%/}"

export PATH="$HOME/.TinyTeX/bin/x86_64-linux:$PATH"
Rscript -e 'tinytex::latexmk("fig1_architecture.tex", engine = "pdflatex", clean = TRUE)'
echo "wrote fig1_architecture.pdf"

# -f: an API error (HTTP 500) must fail the script, not save the JSON error
# body under a .png name.
curl -fsS --max-time 300 -o fig2_gene_sox10.png \
  "$API/plot_brain_png?gene=SOX10&scale=micro"
echo "wrote fig2_gene_sox10.png"

curl -fsS --max-time 300 -o fig3_ontology.png \
  "$API/plot_ontology_png?geneset=GOBP_FOREBRAIN_GENERATION_OF_NEURONS&scale=micro"
echo "wrote fig3_ontology.png"

# Same 13 genes as FAS_DORSAL in make_figures.R (Fischer et al. 2021); no
# log2FC column, so the API scores it with plain ssGSEA.
curl -fsS --max-time 900 -o fig4_genelist.png -X POST \
  -H 'Content-Type: application/json' \
  -d '{"gene_string": "EMX2\nLHX1\nLHX2\nLHX5\nLHX9\nOTX2\nFEZF2\nNEUROD1\nNEUROD2\nNEUROD6\nNEUROG2\nTBR1\nEOMES", "scale": "micro"}' \
  "$API/plot_genelist_png"
echo "wrote fig4_genelist.png"
