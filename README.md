# Zebrafish Gene Expression Deconvolution Using MuSiC

This repository contains R scripts and resources for deconvolving bulk RNA-seq gene expression data in **zebrafish**, using **MuSiC (Multi-subject Single-cell deconvolution)**. The goal is to investigate how **environmental exposures** influence gene expression at the **cell-type level**, using single-cell RNA-seq references from both **pancreas-specific** and **whole-zebrafish** sources.

---

##  Project Objective

To estimate **cell type proportions** in bulk RNA-seq samples of zebrafish and assess how they shift in response to environmental exposures. This is achieved through multimodal analysis by integrating:

- **Bulk RNA-seq gene counts** (exposed vs control)
- **Single-cell RNA-seq reference data** (annotated by cell type)

---


##  Script Details

### `deconvolution_bulk_zebrafish_to_zebrafish_pancreas.R`

- Uses pancreas-specific cell types (e.g., alpha, beta_1, duct)  
- Converts Seurat object (`Control_annotated.rds`) into a `SingleCellExperiment`  
- Performs MuSiC deconvolution and outputs:
  - Jitter plot comparing MuSiC vs NNLS estimates  
  - Heatmap of estimated proportions  
  - Stacked bar plot of cell type distributions  
- **Reference Dataset:**  
  `Control_annotated.rds` downloaded from the Broad Institute's Single Cell Portal ([SCP1549](https://singlecell.broadinstitute.org/single_cell/study/SCP1549))  
  > *Citation:* Singh, S. P., Janjuha, S., Hartmann, T., Kayisoglu, O., Konantz, J., Birke, S., ... & Lickert, H. (2022). A single-cell atlas of de novo β-cell regeneration reveals the contribution of hybrid β/δ cells to diabetes recovery in zebrafish. *Development*, 149(2), dev199853. https://doi.org/10.1242/dev.199853

---

### `deconvolution_bulk_zebrafish_to_whole_zebrafish.R`

- Uses DanioCell 2023 dataset as the scRNA-seq reference (`Daniocell2023_SeuratV4.rds`)  
- Converts both bulk and single-cell data into `ExpressionSet`  
- Performs MuSiC deconvolution using the entire cellular atlas  
- Outputs the same set of visualizations and CSV matrices  
- **Reference Dataset:**  
  `Daniocell2023_SeuratV4.rds` downloaded from the [DanioCell Portal](https://daniocell.nichd.nih.gov/)  
  > *Citation:* White, D. T., Freudenberg, J., Riva, A., Snyder, D. T., Lai, S. L., Ang, K. C., ... & Burgess, S. M. (2023). A zebrafish single-cell transcriptome atlas reveals shared and distinct cellular responses to injury. *bioRxiv*. https://doi.org/10.1101/2023.03.29.534734


---

##  Requirements

- R (version ≥ 4.1)
- CRAN & Bioconductor packages:
```r
install.packages(c("ggplot2", "dplyr", "reshape2", "pheatmap"))
BiocManager::install(c("MuSiC", "Seurat", "SingleCellExperiment", "Matrix", "Biobase", "scuttle"))
