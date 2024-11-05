library(Seurat)
library(SingleR)
library(ggplot2)
library(celldex)
library(SingleCellExperiment)
library(optparse)
library(scater)

# Define command-line options
option_list <- list(
  make_option(c("-w", "--workdir"), type = "character", default = NULL, help = "Working directory", metavar = "CHARACTER"),
  make_option(c("-r", "--rdsPath"), type = "character", default = NULL, help = "Path to Seurat RDS file", metavar = "CHARACTER"),
  make_option(c("-g", "--genome"), type = "character", default = NULL, help = "Genome version", metavar = "CHARACTER"),
  make_option(c("-m", "--markers"), type = "character", default = NULL, help = "Path to marker genes CSV file", metavar = "CHARACTER")
)

# Parse command-line arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Extract arguments
workdir <- args$workdir
rdsPath <- args$rdsPath
genome <- args$genome
markers <- args$markers
markers_path <- "/mnt/ccrsf-ifx/Software/scripts/bin/currentsnake/single_cell/gene_lists/"

# Determine species and reference dataset based on genome
species <- ifelse(genome == "mm10", "Mouse", "Human")
if (genome == "mm10") {
  marker <- paste0(markers_path, 'mouse_gene_list.csv')
} else if (genome == "hg38") {  
  marker <- paste0(markers_path, 'human_gene_list.csv')
} else {
  stop('Please specify genome as mm10 or hg38')
}

# Read input data
seur <- readRDS(rdsPath)
geneList <- read.csv(marker, header = FALSE, row.names = 1, stringsAsFactors = FALSE)

sce <- as.SingleCellExperiment(seur)
sce <- logNormCounts(sce)

# Set working directory
setwd(workdir)

# Function to generate plots
generatePlots <- function(name, pred, seur, markers) {
  # Save prediction results
  saveRDS(pred, paste0('pred_', name, '.rds'))
  
  # Plot score heatmap
  pdf(paste0('heatmap_', name, '.pdf'))
  print(plotScoreHeatmap(pred, clusters = seur$seurat_clusters))
  dev.off()
  
  # Plot score distribution
  pdf(paste0('ScoreDistribution_', name, '.pdf'), height = 15)
  print(plotDeltaDistribution(pred, show = "delta.med", ncol = 3))
  dev.off()
  
  # Add predictions to Seurat metadata
  seur@meta.data[name] <- ''
  seur@meta.data[rownames(pred),name] <- pred$pruned.labels
  Idents(seur) <- seur@meta.data[,name]
  
  # Plot UMAP
  pdf(paste0("UMAP_", name, ".pdf"))
  print(DimPlot(seur, reduction = "umap"))
  dev.off()
  
  # Plot tSNE
  pdf(paste0("TSNE_", name, ".pdf"))
  print(DimPlot(seur, reduction = "tsne"))
  dev.off()
  
  # Create directory for gene list plots
  dir.create("gene_list_plots", showWarnings = FALSE)
  
  # Plot Violin plots for marker genes
  for (i in 1:nrow(markers)) {
    genes.use <- as.character(markers[i, ])
    genes.use <- genes.use[genes.use %in% rownames(seur)]
    numGenes <- length(genes.use)
    if (numGenes > 1) {
      png(paste0('gene_list_plots/ViolinPlot_', name, '_', rownames(markers)[i], '.png'), height = 5, width = (numGenes + 1) * 5, units = 'in', res = 300)
      print(VlnPlot(seur, features = genes.use, ncol = numGenes))
      dev.off()
    }
  }
  
  return(seur)
}

# Function to generate feature plots
generateMarkerPlots <- function(seur, markers) {
  dir.create("gene_list_plots", showWarnings = FALSE)
  for (i in 1:nrow(markers)) {
    genes.use <- as.character(markers[i, ])
    genes.use <- genes.use[genes.use %in% rownames(seur)]
    numGenes <- length(genes.use)
    if (numGenes > 1) {
      png(paste0('gene_list_plots/FeaturePlot_', rownames(markers)[i], '.png'), height = 5, width = (numGenes + 1) * 5, units = 'in', res = 300)
      print(FeaturePlot(seur, features = genes.use, ncol = numGenes))
      dev.off()
    }
  }
}

# Processing based on genome type
if (genome == "mm10") {
  immgen <- celldex::ImmGenData()
  mouserna <- celldex::MouseRNAseqData()

  pred.multi <- SingleR(test = sce, ref = list(IG = immgen, MRNA = mouserna), labels = list(immgen$label.main, mouserna$label.main))
  seur <- generatePlots('immgen_mouserna', pred.multi, seur, geneList)
  write.csv(sort(table(pred.multi$pruned.labels), decreasing = TRUE), 'annotations_Immgen-MouseRNA_general.csv', row.names = FALSE)

  pred.immgen <- SingleR(test = sce, ref = immgen, labels = immgen$label.main)
  seur <- generatePlots('immgen', pred.immgen, seur, geneList)
  write.csv(sort(table(pred.immgen$pruned.labels), decreasing = TRUE), 'annotations_Immgen_general.csv', row.names = FALSE)

  pred.mouserna <- SingleR(test = sce, ref = mouserna, labels = mouserna$label.main)
  seur <- generatePlots('mouserna', pred.mouserna, seur, geneList)
  write.csv(sort(table(pred.mouserna$pruned.labels), decreasing = TRUE), 'annotations_MouseRNA_general.csv', row.names = FALSE)

  write.csv(seur@meta.data[, c('immgen_mouserna', 'immgen', 'mouserna')], 'annotations_compiled_barcode.csv', quote = FALSE)
} else if (genome == "hg38") {
  hpca <- celldex::HumanPrimaryCellAtlasData()
  blueprint <- celldex::BlueprintEncodeData()
  
  pred.multi <- SingleR(test = sce, ref = list(BP = blueprint, HPCA = hpca), labels = list(blueprint$label.main, hpca$label.main))
  seur <- generatePlots('blueprintencode_hpca', pred.multi, seur, geneList)
  write.csv(sort(table(pred.multi$pruned.labels), decreasing = TRUE), 'annotations_BlueprintENCODE-HPCA_general.csv', row.names = FALSE)

  pred.blueprint <- SingleR(test = sce, ref = blueprint, labels = blueprint$label.main)
  seur <- generatePlots('blueprintencode', pred.blueprint, seur, geneList)
  write.csv(sort(table(pred.blueprint$pruned.labels), decreasing = TRUE), 'annotations_BlueprintENCODE_general.csv', row.names = FALSE)

  pred.hpca <- SingleR(test = sce, ref = hpca, labels = hpca$label.main)
  seur <- generatePlots('hpca', pred.hpca, seur, geneList)
  write.csv(sort(table(pred.hpca$pruned.labels), decreasing = TRUE), 'annotations_HPCA_general.csv', row.names = FALSE)

  write.csv(seur@meta.data[, c('blueprintencode_hpca', 'blueprintencode', 'hpca')], 'annotations_compiled_barcode.csv', quote = FALSE)
}

# Generate marker plots
suppressWarnings(generateMarkerPlots(seur, geneList))

# Save the Seurat object
saveRDS(seur, 'seur_cluster_singler.rds')

# Display session information
sessionInfo()
