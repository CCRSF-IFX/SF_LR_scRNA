# Load necessary libraries
library(Seurat)
library(monocle3)
library(SeuratWrappers)
#library(SeuratDisk)
library(SeuratObject)
library(SingleCellExperiment)
library(optparse)
library(Matrix)
library(dplyr)

# Function to run trajectory analysis
runTrajectory <- function(input_file, is_h5 = FALSE, workdir = "./") {
  # Set working directory
  setwd(workdir)
  
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
  # Preprocess the data in Monocle3
  cds <- preprocess_cds(cds, num_dim = 50)
  plot_pc_variance_explained(cds)
  #plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "pseudotime")
  
  # Reduce dimensionality using UMAP
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  plot_cells(cds)
  plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
  
  # Cluster the cells
  #cds <- cluster_cells(cds)
  
  marker_test_res <- top_markers(cds, group_cells_by="partition", 
                                 reference_cells=1000, cores=8)  
  
  
  top_specific_markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(1, pseudo_R2)
  
  top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
  
  
  plot_genes_by_group(cds,
                      top_specific_marker_ids,
                      group_cells_by="partition",
                      ordering_type="maximal_on_diag",
                      max.size=3)
  
  
  top_specific_markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(3, pseudo_R2)
  
  top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
  
  plot_genes_by_group(cds,
                      top_specific_marker_ids,
                      group_cells_by="partition",
                      ordering_type="cluster_row_col",
                      max.size=3)
  
  
  assigned_type_marker_test_res <- top_markers(cds,
                                               group_cells_by="assigned_cell_type",
                                               reference_cells=1000,
                                               cores=8)
  
  garnett_markers <- assigned_type_marker_test_res %>%
    filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
    group_by(cell_group) %>%
    top_n(5, marker_score)
  # Exclude genes that are good markers for more than one cell type:
  garnett_markers <- garnett_markers %>% 
    group_by(gene_short_name) %>%
    filter(n() == 1)
  
  generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")
  
  
  cds <- cluster_cells(cds, reduction_method = "UMAP")  
  
  #cds <- reduce_dimension(cds)
  plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "partition")
  
  # Learn the trajectory graph
  cds <- learn_graph(cds)
  
  plot_cells(cds,
             color_cells_by = "partition",
             label_groups_by_cluster=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE)
  
  # a helper function to identify the root principal points:
  get_earliest_principal_node <- function(cds, time_bin="130-170"){
    cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
  }
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
  
  
  
  
  plot_cells(cds,
             color_cells_by = "pseudotime",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1.5) 
  
  
  
  cds_3d <- reduce_dimension(cds, max_components = 3)
  cds_3d <- cluster_cells(cds_3d)
  cds_3d <- learn_graph(cds_3d)
  cds_3d <- order_cells(cds_3d)#, root_pr_nodes=get_earliest_principal_node(cds))
  
  cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
  
  
  
  # Order the cells along the trajectory
  cds <- order_cells(cds)
  
  # Plot cells with pseudotime
  plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, cell_size = 1, alpha = 0.7)
  ggsave("pseudotime_trajectory.png")
  
  # Plot clusters
  plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1, alpha = 0.7)
  ggsave("clusters_trajectory.png")
  
  # Save the Monocle3 object
  saveRDS(cds, file = ifelse(is_h5, "paga_monocle3.rds", "paga_seurat_monocle3.rds"))
}

# Main function to parse arguments and run the trajectory analysis
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  # Parse arguments
  parser <- argparse::ArgumentParser(description = "Run basic trajectory given either a Seurat RDS object or an h5 file")
  parser$add_argument("input", metavar = "input_file", type = "character", help = "Input file: either Seurat RDS or h5 file")
  parser$add_argument("--h5", dest = "is_h5", action = "store_true", help = "File provided is h5 file")
  parser$add_argument("-w", "--workdir", dest = "workdir", type = "character", default = "./", help = "Working directory to save files in")
  
  # Get parsed arguments
  args <- parser$parse_args(args)
  
  # Run the trajectory analysis
  runTrajectory(args$input, args$is_h5, args$workdir)
}

# Execute the main function if the script is run directly
if (interactive()) {
  main()
} else {
  main(commandArgs(trailingOnly = TRUE))
}
