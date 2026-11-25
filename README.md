# VEGABrain — Gene Expression in Human Developmental Brain

**Live tool:** <https://mcblab.github.io/VEGABrain/> (React front on GitHub Pages, Plumber API on Google Cloud Run).

A Shiny web application for visualizing gene expression throughout the development of the human brain. This tool leverages the **BrainSpan Atlas of the Developing Human Brain** dataset and provides a spatial interface using the `ggseg` package to map expression levels onto the Desikan-Killiany (DK) atlas.

## 🚀 Features

- **Gene Expression Mapping**: Search for any HGNC gene name to see its spatial expression pattern across 5 broad developmental stages (1st Trimester to Adult).
- **GSVA-Enrichr Mode**: Upload a custom list of genes (gene signature) to calculate and visualize its activity across brain regions using GSVA.
- **Per-area significance**: For a gene list, mark the brain areas it is significantly associated with, by two complementary tests (see below).
- **High-Quality Exports**: Download brain maps in SVG or TIFF formats for publications.
- **Interactive UI**: Built with a modern, responsive interface using `shinythemes` (Cerulean).

## 🛠 Prerequisites

The application requires R (≥ 4.3.0) and several system-level dependencies for spatial data processing.

### R Packages

Install the required packages from CRAN, Bioconductor, and GitHub:

```R
# CRAN
install.packages(c("shiny", "dplyr", "tidyr", "ggplot2", "shinythemes", "svglite", "shinycssloaders", "remotes", "ggseg"))

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GSVA")
```

## 📂 Data Preparation

Before running the application, you must process the raw BrainSpan data and calculate the ontology enrichment scores. This is a two-step process based on the provided scripts:

### 1. Prepare the Expression Matrix
Ensure the `genes_matrix_csv/` directory contains:
- `expression_matrix.csv`
- `rows_metadata.csv`
- `columns_metadata.csv`

Run the preparation script from your terminal:
```bash
Rscript preparaDados.R
```
This script generates `dados_otimizados.rds`, which includes mapped brain regions, filtered gene lists, and pre-calculated age groups.

### 2. Generate Ontology Scores (ssGSEA)
The API also requires the `ontologyssGSEA.csv` file to serve the pathway visualization endpoints. This file is generated via the **`ETL.qmd`** document.

1. Open `ETL.qmd` (e.g., in RStudio, VSCode, or via the command line).
2. Ensure you have the `dados_otimizados.rds` file from the previous step.
3. Run the **"Ontology create dataset"** code chunk inside the document. 

This will load the optimized data, fetch gene sets from `msigdbr`, run the GSVA calculations, and output the final `ontologyssGSEA.csv` file to your root directory.

### 3. Per-area significance (no build step)
The Gene List tab can test, for every BrainSpan structure (or region class) and developmental window, whether the submitted list is significantly associated with it. Two tests answer two different questions:

| | `significance=score` (asterisk) | `significance=markers` (diamond, or up/down triangles for log2FC lists) |
|---|---|---|
| question | does the list's ssGSEA score peak in this area beyond random genes? | are this area's specific genes over-represented in the list, and which ones? |
| unit | 87 structure x window cells, 70 with >= 3 donors | region classes (neocortex, hippocampus, amygdala, striatum, thalamus, cerebellum) vs each other, and neocortical areas vs the rest of neocortex; 36 cells with >= 10 markers and >= 3 informative donors |
| test | within-window deviation ranked against 4,000 random gene sets of the same size and expression level; BH; must survive leaving any single donor out | markers from limma per window (`~ 0 + structure + donor`; q < 0.05, log2FC >= 1, mean >= 1); one-sided Fisher's exact test inside the window's expressed MSigDB genes; BH; >= 2 overlapping genes |
| code | `stats/enrichment_stats.R` | `stats/marker_stats.R` |
| calibration | `Rscript stats/validate_stats.R` | `Rscript stats/validate_markers.R` |

`/plot_genelist`, `/plot_genelist_png` and `/data_genelist` take `significance=none|score|markers` (the CSV gains per-region `enriched_*` / `q_*` columns), and `POST /area_enrichment` (`test=score|markers|both`) returns every tested cell as JSON, with the overlapping genes for the marker test. A test needs at least 5 genes found in the data (per direction for the marker test). Everything is computed from `dados_otimizados.rds` at request time: the first request after a start builds the tables (~20 s), and later ones take ~1–3 s. The marker test needs the Bioconductor package `limma`.

### 4. Developmental stage (no build step)
`POST /stage_enrichment` asks in which developmental window a gene list peaks, which the per-area tests cannot see (they compare areas within one window). The list's ssGSEA score per window, with the six PsychENCODE regions (NCX, HIP, AMY, STR, MD, CBC) weighted equally, is centred over the five windows and ranked against 4,000 expression-matched random gene sets: two-sided, BH over the windows, and a hit must survive leaving any single donor out. It runs for the whole brain and inside each region (region x window cells with fewer than 3 donors are reported, not tested); a log2FC list is also split into up- and down-regulated genes. Code: `stats/stage_stats.R`. `/data_genelist` with `significance=all` returns the columns of all three tests in one CSV.

## 🏃 Running the Application

Once the data is prepared, you can launch the app:

```R
shiny::runApp()
```

## 🐳 Docker Support

The full application stack (React Frontend + Plumber API) can be easily served via Docker Compose.

### Using Docker Compose (Recommended)

The easiest way to build and run the complete application is using Docker Compose:

```bash
docker-compose up --build -d
```

This will automatically start both services:
- **React Frontend**: `http://localhost:5173`
- **Plumber API**: `http://localhost:33857` *(Swagger Docs available at `http://localhost:33857/__docs__/`)*

### Live Reloading during Development

Docker Compose is configured to automatically sync your code changes during development without needing full container rebuilds. 

Instead of `docker-compose up`, use the `watch` command:

```bash
docker compose watch
```

With this running:
- Any edits saved to `plumber.R` will instantly sync and restart the API container.
- Any edits saved in the `front/` directory will instantly trigger Vite's Hot Module Replacement (HMR) and update the frontend live in your browser.

### Using Docker Manually (API Only)

If you wish to run *only* the backend API, you can build and run it directly:

```bash
# Build the API image
docker build -t devbrain-markers-api -f Dockerfile.api .

# Run the API container
docker run -p 33857:33857 devbrain-markers-api
```

## 📊 Data Source

Data is sourced from the [BrainSpan Atlas of the Developing Human Brain](https://www.brainspan.org/), specifically the RNA-Seq RPKM values averaged to genes.

## 🤝 Acknowledgements

Made by the **[MCB Lab at UFRN](https://mcblab.github.io/)**.

If you use this tool in your research, please cite the BrainSpan Atlas and the `ggseg` package.
