# Load necessary libraries
library(monocle3)
library(dplyr)
library(ggplot2)
library(Seurat)

# Load the Seurat object
if (is_h5) {
  seurat_object <- LoadH5Seurat(input_file)
} else {
  seurat_object <- readRDS(input_file)
}

# Ensure the object is a Seurat object
if (!inherits(seurat_object, "Seurat")) {
  stop("The input file does not contain a valid Seurat object.")
}

# Standard pre-processing
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Run UMAP for dimensionality reduction
seurat_object <- RunUMAP(seurat_object, dims = 1:20)

# Convert Seurat object to Monocle3 object
cds <- as.cell_data_set(seurat_object)

# Cluster the cells
cds <- cluster_cells(cds) 

# Read the marker gene list from the TSV file
marker_genes_df <- read.table("human_gene_list.csv", sep = ",", header = FALSE, stringsAsFactors = FALSE)

# Loop through each row in the marker_genes_df to perform differential expression analysis
for (i in 1:nrow(marker_genes_df)) {
  # Extract group name and genes
  group_name <- marker_genes_df[i, 1]
  genes <- unlist(marker_genes_df[i, -1])
  genes <- genes[genes != ""]  # Remove empty elements
  
  # Subset the cds object to include only the genes in the current group
  genes_in_cds <- rowData(cds)$gene_short_name %in% genes
  cds_subset <- cds[genes_in_cds, ]
  
  if (sum(genes_in_cds) == 0) {
    message(paste("No genes found in cds for group", group_name))
    next
  }
  
  # Perform differential expression analysis
  gene_fits <- fit_models(cds_subset, model_formula_str = "~cluster")
  
  # Extract the coefficient table
  fit_coefs <- coefficient_table(gene_fits)
  
  # Filter significant genes (adjust q_value threshold as needed)
  significant_genes <- fit_coefs %>% filter(q_value < 0.05) %>%
    select(gene_short_name, term, q_value, estimate)
  
  if (nrow(significant_genes) == 0) {
    message(paste("No significant genes found for group", group_name))
    next
  }
  
  # Print significant genes for the current group
  print(paste("Significant genes for group:", group_name))
  print(significant_genes)
  
  # Generate and save violin plots for the significant genes
  plot <- plot_genes_violin(cds_subset, group_cells_by = "cluster", genes = significant_genes$gene_short_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Violin Plot for Group:", group_name))
  
  # Save the plot
  ggsave(paste0("ViolinPlot_", group_name, ".png"), plot = plot, width = 10, height = 7)
  
  # Optionally, save the significant genes to a CSV file
  write.csv(significant_genes, paste0("SignificantGenes_", group_name, ".csv"), row.names = FALSE)
}
