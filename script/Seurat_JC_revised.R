#lib_paths <- c("/mnt/nasapps/development/R/r_libs/4.2.2", "/home/ccrsfifx/R/x86_64-redhat-linux-gnu-library/4.3")
# Set the library paths
#.libPaths(lib_paths)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(cluster)
library(optparse)

roundUp <- function(x, to = 10) {
  to * (x %/% to + as.logical(x %% to))
}

# Define command-line options
option_list <- list(
  make_option(c("-w", "--workdir"), type = "character", default = NULL, help = "Working directory", metavar = "CHARACTER"),
  make_option(c("-p", "--path"), type = "character", default = NULL, help = "Path to input data", metavar = "CHARACTER"),
  make_option(c("-n", "--name"), type = "character", default = NULL, help = "Dataset name", metavar = "CHARACTER"),
  make_option(c("-g", "--genome"), type = "character", default = NULL, help = "Genome version", metavar = "CHARACTER")
)

# Parse command-line arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

if (!is.null(args$help)) {
  args <- args[names(args) != "help"]
}

print(args)

# Extract key components from parsed arguments
workdir <- args$workdir
path <- args$path
name <- args$name
genome <- args$genome

# Show the structure of the parsed arguments
str(args)

# Generate Seurat object
#if (file.exists(path)) {
  # If the path is a file, read the matrix from the TSV
#  matrix_data <- read.delim(path, row.names = 1, check.names = FALSE)
#} else 
if (dir.exists(path)) {
  # If the path is a directory, read the 10X Cell Ranger matrix
  matrix_data <- Read10X(data.dir = path)
} else {
  stop("Invalid path. The provided path is neither a file nor a directory.")
}

seur <- CreateSeuratObject(counts = matrix_data, min.cells = 3, min.features = 100, project = name)

# Housekeeping
seur <- NormalizeData(seur)
seur <- FindVariableFeatures(seur)
seur <- ScaleData(seur)

# Set working directory
setwd(workdir)

# Save gene expression matrix
genematrix <- as.matrix(GetAssayData(seur, layer = "data"))
write.csv(genematrix, file = "Gene_Expression_Matrix.csv", quote = FALSE, row.names = TRUE)

# Calculate and save the number of genes detected per cell
gene_counts_per_cell <- data.frame(Barcodes = colnames(genematrix), Genes = colSums(genematrix != 0))
write.csv(gene_counts_per_cell, file = "Number_of_genes_per_barcode.csv", quote = FALSE, row.names = FALSE)

# Plot Genes Per Barcode
plot_genes_per_barcode <- function(x) {
  mean_x <- mean(x)
  median_x <- median(x)
  density_x <- density(x)
  
  hist_breaks <- 20
  hist_color <- "grey"
  hist_x <- hist(x, breaks = hist_breaks, plot = FALSE)
  xmax <- max(roundUp(hist_x$breaks, 500))
  xmin <- 0
  
  png("GenesPerBarcodeEditedPlot.png", height = 7, width = 7, units = 'in', res = 300)
  hist(x, breaks = hist_breaks, col = hist_color, labels = TRUE, yaxt = 'n', 
       xlab = "Number of Genes", ylab = "Number of Cells", main = "Gene Count Per Cell",
       xlim = c(xmin, xmax), xaxp = c(xmin, xmax, xmax / 500))
  
  abline(v = mean_x, lwd = 2, col = "blue")
  abline(v = median_x, lwd = 2, col = "orange")
  
  par(new = TRUE)
  hist(x, breaks = hist_breaks, col = hist_color, freq = FALSE, ann = FALSE, 
       xlim = c(xmin, xmax), xaxp = c(xmin, xmax, xmax / 500))
  lines(density_x, col = "red", lwd = 2)
  
  mtext(paste("Mean = ", round(mean_x, 2)), side = 3, adj = 0.24, col = "blue", lwd = 2)
  mtext(paste("Median = ", round(median_x, 2)), side = 3, adj = 0.65, col = "orange")
  
  dev.off()
}

# Plot Cells Per Gene
plot_cells_per_gene <- function(x, num_genes) {
  mean_x <- mean(x)
  median_x <- median(x)
  density_x <- density(x)
  
  hist_breaks <- 20
  hist_color <- "grey"
  hist_x <- hist(x, breaks = hist_breaks, plot = FALSE)
  xmax <- max(roundUp(hist_x$breaks, 100))
  xmin <- 0
  
  png("CellsPerGenePlot.png", height = 7, width = 7, units = 'in', res = 300)
  hist(x, breaks = hist_breaks, col = hist_color, labels = TRUE, yaxt = 'n',
       xlab = "Number of Cells", ylab = "Number of Genes", main = "Gene Expression",
       xlim = c(xmin, xmax), xaxp = c(xmin, xmax, xmax / 100))
  
  abline(v = mean_x, lwd = 2, col = "blue")
  abline(v = median_x, lwd = 2, col = "orange")
  
  par(new = TRUE)
  hist(x, breaks = hist_breaks, col = hist_color, freq = FALSE, ann = FALSE,
       xlim = c(xmin, xmax), xaxp = c(xmin, xmax, xmax / 100))
  lines(density_x, col = "red", lwd = 2)
  
  mtext(paste("Mean = ", round(mean_x, 2)), side = 3, adj = 0.1, col = "blue", lwd = 2)
  mtext(paste("Median = ", round(median_x, 2)), side = 3, adj = 0.4, col = "orange")
  mtext(paste("Total Genes = ", num_genes), side = 3, adj = 0.80, col = "purple")
  
  par(new = FALSE)
  dev.off()
}

# Run the plots
x <- gene_counts_per_cell$Genes
plot_genes_per_barcode(x)

genes_per_cell <- data.frame(Genes = rownames(genematrix), Cells = rowSums(genematrix != 0))
write.csv(genes_per_cell, file = "Number_of_cells_per_gene.csv", quote = FALSE, row.names = TRUE)
genes_for_plot <- genes_per_cell[genes_per_cell$Cells > 0, ]
x <- genes_for_plot$Cells
num_genes <- nrow(genes_for_plot)
plot_cells_per_gene(x, num_genes)

# Determine the percentage of mitochondrial genes based on the genome type
if (genome == "mm10") {
  seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern = "^mt-")
} else if (genome == "hg19" || genome == "hg38") {
  seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
}

# Create a violin plot to visualize quality control metrics
png("VlnPlot.png", height = 7, width = 7, units = 'in', res = 300)
VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# Generate scatter plots for feature relationships
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Save scatter plots to a PNG file using patchwork for combining plots
png("PreFilter_Gene_Plot.png", height = 7, width = 10, units = 'in', res = 300)
print(plot1 + plot2)
dev.off()

# Run PCA and UMAP in Seurat
seur <- RunPCA(seur, npcs = 50)
seur <- RunUMAP(seur, dims = 1:10)

# Plot the elbow to estimate optimal PCs
png("ElbowPlot.png", height = 7, width = 7, units = 'in', res = 300)
ElbowPlot(seur, ndims = 50)
dev.off()

# Run JackStraw analysis
seur <- JackStraw(seur, num.replicate = 100)
seur <- ScoreJackStraw(seur, dims = 1:20)
numlist <- which(seur@reductions$pca@jackstraw@overall.p.values < 0.05)
numPCs <- max(10, length(numlist))

# Perform clustering for Doublet
seur <- FindNeighbors(seur, dims = 1:numPCs)
seur <- FindClusters(seur, resolution = 0.5)

# Homotypic Doublet Proportion Estimate
annotations <- seur@meta.data$seurat_clusters
homotypic_prop <- modelHomotypic(annotations)  # Estimate homotypic doublet proportion

# Calculate expected number of doublets
nExp_poi <- round(0.075 * nrow(seur@meta.data))  # 7.5% doublet formation rate
nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

# Use paramSweep to determine the optimal pK
sweep.res.list <- paramSweep(seur, PCs = 1:numPCs, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Plot to visualize the optimal pK
pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))

# Run DoubletFinder with the identified pK
seur <- doubletFinder(seur, PCs = 1:numPCs, pN = 0.25, pK = pK, nExp = nExp_poi_adj)

# Add doublet status to metadata
column_name <- paste0("DF.classifications_0.25_", pK, "_", nExp_poi_adj)
seur$doublet_status <- seur@meta.data[[column_name]]

# Visualize with UMAP
png("UMAP_Singlet_Doublet_Plot.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(seur, reduction = "umap", group.by = "doublet_status", label = TRUE)
dev.off()

# Filtering Parameters
p_hi <- 1e-3  # High p-value for doublet filtering
p_lo <- 1e-2  # Low p-value for poor libraries

# Create a data frame with library size and gene detection
cell_stats <- data.frame(libSize = seur$nCount_RNA, geneDetect = seur$nFeature_RNA)

# Filter out zero or NA values
data_vector <- cell_stats$libSize
data_vector <- data_vector[data_vector > 0 & !is.na(data_vector)]

# Fit a Normal Distribution
mean_data <- mean(data_vector)
var_data <- var(data_vector)
sd_data <- sqrt(var_data)

# Define limits for filtering based on UMI counts
umi_upper_limit <- qnorm(1 - p_hi, mean = mean_data, sd = sd_data)
umi_lower_limit <- qnorm(p_lo, mean = mean_data, sd = sd_data)

# Define limits for filtering based on gene counts
gene_upper_limit <- qnorm(1 - p_hi, mean = mean(cell_stats$geneDetect), sd = sqrt(var(cell_stats$geneDetect)))
gene_lower_limit <- qnorm(p_lo, mean = mean(cell_stats$geneDetect), sd = sqrt(var(cell_stats$geneDetect)))

# Identify cells to be filtered based on UMI and gene counts
temp_doublets <- (cell_stats$libSize > umi_upper_limit) | (cell_stats$geneDetect > gene_upper_limit)
temp_crapLibs <- (cell_stats$libSize < umi_lower_limit) | (cell_stats$geneDetect < gene_lower_limit)

# Mitochondrial Content Analysis
cur_mad <- mad(seur$percent.mito)
cur_med <- median(seur$percent.mito)
mito_upper_limit <- cur_med + 4 * cur_mad

# Store Filter Information
filter_info <- data.frame(
  Doublets = sum(temp_doublets), 
  Poor_Quality = sum(temp_crapLibs), 
  Mitochondrial = sum(seur$percent.mito > mito_upper_limit),
  Total_Filtered = sum(temp_doublets | temp_crapLibs | seur$percent.mito > mito_upper_limit)
)
write.table(t(filter_info), 'FilterNumbers.csv', sep = ',', quote = FALSE, row.names = TRUE)

# Subset Seurat Object for Filtering
seur <- subset(seur, subset = nFeature_RNA > gene_lower_limit & nFeature_RNA < gene_upper_limit)
seur <- subset(seur, subset = nCount_RNA > umi_lower_limit & nCount_RNA < umi_upper_limit)
seur <- subset(seur, subset = percent.mito < mito_upper_limit)

# Create and save scatter plots for quality control
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Save plots to PNG file
png("PostFilter_Gene_Plot.png", height = 7, width = 10, units = 'in', res = 300)
print(plot1 + plot2)
dev.off()

# Normalize with SCTransform
seur <- SCTransform(seur, vars.to.regress = "percent.mito", return.only.var.genes = FALSE)

# Save expression matrices
write.csv(as.matrix(GetAssayData(seur, layer = "counts")), file = "Filtered_Gene_Expression_Matrix.csv", quote = FALSE, row.names = TRUE)
write.csv(as.matrix(GetAssayData(seur)), file = "Filtered_Normalized_Gene_Expression_Matrix.csv", quote = FALSE, row.names = TRUE)

saveRDS(seur, file = "seur_preprocessed_object.rds")

# Identify the top 10 variable features
top10 <- head(VariableFeatures(seur), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seur)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Save the combined plot
png("VariableFeatures.png", height = 7, width = 10, units = 'in', res = 300)
print(plot1 + plot2)
dev.off()

# Perform PCA with selected variable features
seur <- RunPCA(seur, features = VariableFeatures(seur), npcs = max(20, numPCs), verbose = FALSE)

# Plot PCA dimension loadings
pdf("VizPCAPlot.pdf")
for (i in seq(1, 20, 2)) {
  print(VizDimLoadings(seur, dims = i:(i + 1)))
}
dev.off()

# Plot PCA results in pairs
pdf("AllPCAPlot.pdf")
for (i in 1:10) {
  print(DimPlot(seur, dims = c(i, i + 1), reduction = "pca"))
}
dev.off()

# Generate PC heatmaps
pdf("PC_HeatmapPlot.pdf")
for (i in 1:10) {
  DimHeatmap(seur, dims = i, cells = 500, balanced = TRUE)
}
dev.off()

# Find neighbors and clustering for different resolutions
seur <- FindNeighbors(seur, dims = 1:numPCs)
seur <- RunTSNE(seur, dims = 1:numPCs)
seur <- RunUMAP(seur, dims = 1:numPCs)

write.csv(Embeddings(seur, reduction = "tsne"), file = "tSNECoordinates.csv")
write.csv(Embeddings(seur, reduction = "umap"), file = "UMAPCoordinates.csv")

# Set up the resolutions to explore
resolutions <- c(0.1, 0.3, 0.6, 0.8)

# Initialize lists for TSNE and UMAP plots
tsnePlots <- list()
umapPlots <- list()
runRes <- c()

# Iterate through resolutions to find clusters and generate plots
for (res in resolutions) {
  seur <- FindClusters(seur, resolution = res)
  
  # Create t-SNE plot with the current resolution
  tsne <- DimPlot(seur, reduction = "tsne") +
    ggtitle(paste(numPCs, "PCs_res", res, sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Store and save the t-SNE plot
  tsnePlots[[as.character(res)]] <- tsne
  png(paste("TSNEPlot_", numPCs, "PCs_", res, ".png", sep = ""), height = 7, width = 7, units = 'in', res = 300)
  print(tsne)
  dev.off()
  
  # Create UMAP plot with the current resolution
  umap <- DimPlot(seur, reduction = "umap") +
    ggtitle(paste(numPCs, "PCs_res", res, sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Store and save the UMAP plot
  umapPlots[[as.character(res)]] <- umap
  png(paste("UMAPPlot_", numPCs, "PCs_", res, ".png", sep = ""), height = 7, width = 7, units = 'in', res = 300)
  print(umap)
  dev.off()
  
  # Find markers for each cluster
  seur.markers <- FindAllMarkers(seur, logfc.threshold = 0.25, only.pos = TRUE)

    # Check if seur.markers is not empty and contains the cluster column
  if (!is.null(seur.markers) && "cluster" %in% colnames(seur.markers)) {
    write.csv(seur.markers %>% group_by(cluster) %>% top_n(100, p_val),
              paste("top100markers_", numPCs, "_res", res, ".csv", sep = ""))
  } else {
    warning(paste("No clusters found for resolution", res))
  } 

  runRes <- append(runRes, res)
  saveRDS(seur.markers, paste("markers_res", res, ".rds", sep = ""))
}

# Save the Seurat object
saveRDS(seur, file = "seur_cluster_object.rds")

# Create PDF with all TSNE plots
pdf("TSNEPlots.pdf")
for (tsne in tsnePlots) {
  print(tsne)
}
dev.off()

# Create PDF with all UMAP plots
pdf("UMAPPlots.pdf")
for (umap in umapPlots) {
  print(umap)
}
dev.off()

for (res in runRes) {
  # Extract PCA coordinates and cluster identities
  coord <- Embeddings(seur, reduction = 'pca')[, 1:numPCs]
  Idents(seur) <- seur@meta.data[[paste0('SCT_snn_res.', res)]]
  clusters <- Idents(seur)
  
  # Calculate Euclidean distances and silhouette values
  d <- dist(coord, method = "euclidean")
  sil <- silhouette(as.numeric(clusters), d)
  
  # Create and save the silhouette plot
  pdf(paste0("SilhouettePlot_res", res, ".pdf"))
  plot(sil, col = as.factor(clusters[order(clusters, decreasing = FALSE)]), 
       main = paste("Silhouette Plot - Resolution", res), lty = 2)
  abline(v = mean(sil[, 3]), col = "red4", lty = 2)  # Mean silhouette value
  dev.off()
}

# Remove metadata for resolutions that didn't generate markers
unsuccessful_res <- setdiff(resolutions, runRes)
for (res in unsuccessful_res) {
  seur@meta.data[[paste0('SCT_snn_res.', res)]] <- NULL  # Clean up unused resolutions
}

# Write minimum successful resolution to a text file
write(min(runRes), "minRes.txt")
