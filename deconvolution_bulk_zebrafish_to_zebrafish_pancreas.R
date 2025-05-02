################################################################################
# ZEBRAFISH BULK TO ZEBRAFISH PANCREAS
################################################################################

# -------------------- PART 0: Load Required Libraries -------------------------
library(Seurat)
library(SingleCellExperiment)
library(MuSiC)
library(Matrix)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(scuttle)

# -------------------- PART 1: Process Bulk RNA-seq Data -----------------------
setwd("~/Desktop/Music/RawGeneCounts")
files <- list.files(pattern = "gene_counts_data_.*", full.names = TRUE)

process_file <- function(file) {
  data <- read.delim(file, header = TRUE)
  data <- data %>% select(-Chr, -Start, -End, -Strand, -Length)
  return(data)
}

raw_gene_counts <- lapply(files, process_file)
names(raw_gene_counts) <- gsub("gene_counts_data_|\\.txt$", "", basename(files))

# -------------------- PART 2: Prepare Bulk RNA-seq Data -----------------------
data_file <- 1  # ****CHANGE as necessary
selected_data <- raw_gene_counts[[data_file]] 
file_name <- names(raw_gene_counts)[data_file]
print(file_name)

selected_data <- selected_data[!is.na(selected_data$Geneid), ]
selected_data <- selected_data[!duplicated(selected_data$Geneid), ]
rownames(selected_data) <- selected_data$Geneid
selected_data <- selected_data[, -which(names(selected_data) == "Geneid")]
bulk_matrix <- as.matrix(selected_data)
print(head(bulk_matrix))

# -------------------- PART 3: Load and Prep Single-Cell Data ------------------
setwd("~/Desktop/Music")
zebrafish.sce <- readRDS("Control_annotated.rds")
DefaultAssay(zebrafish.sce) <- "RNA"
zebrafish.sce <- as.SingleCellExperiment(zebrafish.sce)

clusters_col <- 'ident'
samples_col <- 'orig.ident'

# Use pancreas-specific cell types
selected_cell_types <- c('alpha', 'beta_1', 'delta2', 'delta1', 'beta_2', 
                         'endothelial', 'gip', 'epsilon','duct', 'acinar', 'immune')

# Filter to shared genes
intersect_genes <- intersect(rownames(bulk_matrix), rownames(zebrafish.sce))
bulk_matrix <- bulk_matrix[intersect_genes, ]
zebrafish.sce <- zebrafish.sce[intersect_genes, ]

# Assign pseudo-samples for MuSiC
colData(zebrafish.sce)$pseudo_sample <- sample(1:3, size = ncol(zebrafish.sce), replace = TRUE)
samples_col <- "pseudo_sample"

# -------------------- PART 4: Run MuSiC Deconvolution -------------------------
Est.prop.bulk_matrix <- music_prop(
  bulk.mtx = bulk_matrix, 
  sc.sce = zebrafish.sce, 
  clusters = clusters_col, 
  samples = samples_col,
  select.ct = selected_cell_types,
  verbose = TRUE
)

# -------------------- PART 5: Jitter Plot -------------------------------------
jitter.fig <- Jitter_Est(
  list(
    data.matrix(Est.prop.bulk_matrix$Est.prop.weighted),
    data.matrix(Est.prop.bulk_matrix$Est.prop.allgene)
  ),
  method.name = c('MuSiC', 'NNLS'),
  title = 'Jitter plot of Est Proportions'
)

print(jitter.fig)
ggsave(filename = paste0("jitter_plot_est_prop_pancreas_", file_name, ".png"),
       plot = jitter.fig, width = 8, height = 6, dpi = 300)

# -------------------- PART 6: Heatmap ------------------------------------------
heatmap_data <- data.matrix(Est.prop.bulk_matrix$Est.prop.weighted)
write.csv(heatmap_data, file = paste0("est_prop_pancreas_matrix_", file_name, ".csv"), row.names = TRUE)

heatmap <- pheatmap(heatmap_data,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    color = colorRampPalette(c("navy", "white", "firebrick"))(50),
                    main = "Heatmap of Estimated Cell Type Proportions",
                    fontsize_row = 8,
                    fontsize_col = 8)

print(heatmap)
ggsave(filename = paste0("heatmap_est_prop_pancreas_", file_name, ".png"),
       plot = heatmap, width = 6, height = 6, dpi = 300)

# -------------------- PART 7: Stacked Bar Plot ---------------------------------
heatmap_data_df <- as.data.frame(heatmap_data)
heatmap_data_df$Sample <- rownames(heatmap_data_df)
heatmap_data_long <- melt(heatmap_data_df, id.vars = "Sample",
                          variable.name = "CellType", value.name = "Proportion")

stacked_bar_plot <- ggplot(heatmap_data_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportions of Cell Types Across Samples",
       x = "Samples", y = "Proportion") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()

print(stacked_bar_plot)
ggsave(filename = paste0("stacked_bar_plot_est_prop_pancreas_", file_name, ".png"),
       plot = stacked_bar_plot, width = 6, height = 6, dpi = 300)
