## Club-like cells persist throughout treatment and interact with immunosuppressive myeloid cells in the prostate tumor microenvironment


**Author:** Antti Kiviaho

**Email:** antti.kiviaho@tuni.fi

**Date:** 20.6.2024

This repository contains necessary code to perform data processing and analyses related to the manuscript.

### Package versions:
- **single_cell_analysis_environment.yml** – Containts package versions used in the analysis of spatial transcriptomics data.
- **spatial_analysis_environment.yml** – Containts package versions used in the analysis of single-cell data.

### File formats:
- **Jupyter notebooks (.ipynb)** that mostly contain plotting-related functions and that require less computational resources.
- **Python (.py)** files that mostly contain code related to intensive computational tasks run on more resources.
- **Bash (.sh)** files for submitting jobs to the slurm resource manager.

## Detailed description of file contents: 

### Single-cell preprocessing and annotation 
- **single-cell-preprocessing.ipynb** – Preprocessing of single-cell datasets to attain uniform format.

- **single-cell-quality-control.ipynb** – Gene filtering, doublet removal, and normalization steps carried out on each dataset individually.

- **scvi-integrate-sc-reference.py** – scVI-based integration of 7 preprocessed single-cell datasets to find a common embedding.

- **single-cell-post-integration.ipynb** – Gene marker-based cell type annotation on the integrated dataset. Removal of sample-specfic clusters.

- **nmf-on-celltype.py** – Non-negative matrix factorization-based annotation of cell type-specific gene expression modules.

- **single-cell-immune-celltypist.ipynb** – Celltypist-based annotation of cells annotated as immune cells on the previous round (broad annotation).

- **single-cell-annotation-from-nmf-best-iteration.ipynb** – Non-immune cell subtype annotation based on the NMF gene expression modules.

- **c2l-infer-cell-type-reference.py** – Annotation-based regression of cell2location-compatible cell type signatures.

### Spatial transcriptomics data deconvolution and analyses
- **spatial-qc-and-normalization.ipynb** – Quality control and preprocessing of spatial transcriptomics data.

- **c2l-cell-type-mapping.py** – Spatial transcriptomics data deconvolution using the cell type reference created from single-cell data.

- **spatial-post-c2l-mapping.ipynb** – Division of spatial transcriptomics data into single-cell mapping-based (SCM) regions.

- **spatial_scanpy_enrichment_and_plots.py** – Scanpy score calculation of gene sets.

- **spatial-custom-list-enrichment.ipynb** – Gene set scoring-based plotting.

- **spatial-gene-expression-analysis.ipynb** – Gene expression analysis and plots of spatial transcriptomics data.

- **plot_spatial_progenitor_expression_dotplot.py** – Plot 2c generation.

- **plot_spatial_chemokine_expression_dotplot.py** – Plot 3d generation.

- **spatial_mapping_based_ligand_receptor_analysis.py** – Ligand-receptor interaction analysis.

- **spatial-mapping-based-clusters-receptor-ligand-analyses.ipynb** – Ligand-receptor interaction analysis-based plotting (Figures 4e-4g).

- **spatial_to_pseudobulk.py** – Generating pseudobulk expression data from spatial transcriptomics data.

- **single-cell-club-signature-comparison.ipynb** – Analysis of metastatic prostate cancer single-cell transcriptomics daata (Figures 5a-5b).

- **met-sample-gene-list-scoring.ipynb** – Analysis of metastatic prostate cancer spatial transcriptomics samples (Figures 5c-5e).

- **bulk-rna-survival-analysis.ipynb** – TCGA and SU2C gene expression data analyses (Figures 5g, 5h).

- **spatial_pseudobulk_analysis.ipynb** – Pseudobulk spatial transcriptomics data analysis and plots (Figures 5f, 5i).