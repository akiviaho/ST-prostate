### Club-like cells persist throughout treatment and interact with immunosuppressive myeloid cells in the prostate tumor microenvironment


**Author:** Antti Kiviaho

**Email:** antti.kiviaho@tuni.fi

**Last modified:** 30.7.2024

This repository contains necessary code to reproduce results in Kiviaho et al. 2024.

### Package versions:
- **single_cell_analysis_environment.yml** – Containts package versions used in the analysis of spatial transcriptomics data.
- **spatial_analysis_environment.yml** – Containts package versions used in the analysis of single-cell data.

## Spatial Transcriptomics data

Files related to ST data analysis

### Data preprocessing and computations

- **spatial-qc-and-normalization.ipynb** – Quality control and preprocessing of spatial transcriptomics data.

- **spatial-c2l-cell-type-reference-mapping.py** – Spatial transcriptomics data deconvolution using the cell type reference created from single-cell data.

- **spatial-post-c2l-cell-type-mapping.ipynb** – Division of spatial transcriptomics data into single-cell mapping-based (SCM) regions.

- **spatial-gene-set-scoring.py** – Gene set scoring on spatial data.

- **spatial-region-ligand-receptor-signaling-analysis.py** – Ligand-receptor interaction analysis.

- **spatial-to-pseudobulk.py** – Generating pseudobulk expression data from spatial transcriptomics data.

### Data analysis and plotting results

- **spatial-gene-expression-analysis.ipynb** – Gene expression analysis and plots of spatial transcriptomics data. (Figures 1c, 2d)

- **spatial-gene-set-score-analysis.ipynb** – Gene set scoring-based plotting (Figures 3a, 3b, 3c, 4a, 4b, 4c, 4h).

- **spatial-dotplot-ar-basal-club-markers.py** – Plot 2c generation.

- **spatial-dotplot-chemokine-expression.py** – Plot 3d generation.

- **spatial-mapping-based-clusters-receptor-ligand-analyses.ipynb** – Ligand-receptor interaction analysis-based plotting (Figures 4e, 4f, 4g).

- **spatial-metastatic-tumor-sample-analysis.ipynb.ipynb** – Analysis of metastatic prostate cancer spatial transcriptomics samples (Figures 5d, 5e).

- **spatial-pseudobulk-data-analysis.ipynb** – Pseudobulk spatial transcriptomics data analysis and plots (Figures 5f, 5i).

- **public-sc-and-bulk-data-analysis.ipynb** – Analysis of public data He et al. 2021 + TCGA and SU2C gene expression data analyses (Figures 5a, 5b, 5g, 5h).

## Single-cell data analysis 

Files related to single-cell data used in the article

- **single-cell-preprocessing.ipynb** – Preprocessing of single-cell datasets to attain uniform format.

- **single-cell-quality-control.ipynb** – Gene filtering, doublet removal, and normalization steps carried out on each dataset individually.

- **single-cell-scvi-integrate.py** – scVI-based integration of 7 preprocessed single-cell datasets to find a common embedding.

- **single-cell-post-integration.ipynb** – Gene marker-based cell type annotation on the integrated dataset. Removal of sample-specfic clusters.

- **single-cell-nmf-analysis.py** – Non-negative matrix factorization-based annotation of cell type-specific gene expression modules.

- **single-cell-immune-celltypist.ipynb** – Celltypist-based annotation of cells annotated as immune cells on the previous round (broad annotation).

- **single-cell-gene-module-based-annotation.ipynb** – Non-immune cell subtype annotation based on the NMF gene expression modules.

- **single-cell-generate-c2l-celltype-reference.py** – Annotation-based regression of cell2location-compatible cell type signatures.
