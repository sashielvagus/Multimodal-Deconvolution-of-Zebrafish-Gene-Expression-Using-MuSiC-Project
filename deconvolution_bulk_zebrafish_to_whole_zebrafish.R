################################################################################
# ZEBRAFISH BULK TO WHOLE ZEBRAFISH 
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
library(Biobase)

# -------------------- PART 1: Process Bulk RNA-seq Data ------------------------
setwd("~/Gene_count_files")
files <- list.files(pattern = "gene_counts_data_.*", full.names = TRUE)

process_file <- function(file) {
  data <- read.delim(file, header = TRUE)
  data <- data %>% select(-Chr, -Start, -End, -Strand, -Length)
  return(data)
}

raw_gene_counts <- lapply(files, process_file)
names(raw_gene_counts) <- gsub("gene_counts_data_|\\.txt$", "", basename(files))

# -------------------- PART 2: Prepare Bulk RNA-seq Data ------------------------
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

# -------------------- PART 3: Load Daniocell Single-Cell Data ------------------
setwd("~/")  
daniocell <- readRDS("Daniocell2023_SeuratV4.rds")
daniocell.sce <- as.SingleCellExperiment(daniocell)

clusters_col <- "ident"
samples_col <- "orig.ident"

# -------------------- Filter to Shared Genes ------------------------
intersect_genes <- intersect(rownames(bulk_matrix), rownames(daniocell.sce))
bulk_matrix <- bulk_matrix[intersect_genes, , drop = FALSE]
daniocell.sce <- daniocell.sce[intersect_genes, , drop = FALSE]

# -------------------- Convert Bulk to ExpressionSet ------------------------
bulk_metadata <- data.frame(sampleID = colnames(bulk_matrix), row.names = colnames(bulk_matrix))
bulk_eset <- ExpressionSet(
  assayData = bulk_matrix, 
  phenoData = AnnotatedDataFrame(bulk_metadata)
)
print(bulk_eset)

# -------------------- Convert SC to ExpressionSet ------------------------
sc_counts <- assay(daniocell.sce, "logcounts")
sc_metadata <- colData(daniocell.sce)
rownames(sc_metadata) <- colnames(sc_counts)

sc_eset <- ExpressionSet(
  assayData = as.matrix(sc_counts), 
  phenoData = AnnotatedDataFrame(as.data.frame(sc_metadata))
)
print(sc_eset)

# -------------------- PART 4: Run MuSiC Deconvolution ------------------------
Est.prop.bulk_matrix <- music_prop(
  bulk.eset = bulk_eset,  
  sc.eset = sc_eset,      
  clusters = clusters_col, 
  samples = samples_col,
  verbose = TRUE
)

# View results
print(Est.prop.bulk_matrix$Est.prop.weighted)

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

ggsave(filename = paste0("jitter_plot_est_prop_wholeZebrafish_", file_name, ".png"),
       plot = jitter.fig, width = 8, height = 6, dpi = 300)

# -------------------- PART 6: Heatmap ------------------------------------------
heatmap_data <- data.matrix(Est.prop.bulk_matrix$Est.prop.weighted)
write.csv(heatmap_data, file = paste0("est_prop_wholeZebrafish_matrix_", file_name, ".csv"), row.names = TRUE)

heatmap <- pheatmap(heatmap_data,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    color = colorRampPalette(c("navy", "white", "firebrick"))(50),
                    main = "Heatmap of Estimated Cell Type Proportions",
                    fontsize_row = 8,
                    fontsize_col = 8)
print(heatmap)

ggsave(filename = paste0("heatmap_est_prop_wholeZebrafish_", file_name, ".png"),
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
ggsave(filename = paste0("stacked_bar_plot_est_prop_wholeZebrafish_", file_name, ".png"),
       plot = stacked_bar_plot, width = 6, height = 6, dpi = 300)