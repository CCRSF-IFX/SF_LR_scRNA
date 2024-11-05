library("enrichR")
library("dplyr")

sessionInfo()

args<-commandArgs(TRUE);
show(args);
workdir<-args[1];
genome<-args[2];
seuratPath<-args[3];
seuratClusterObject<-args[4];

paths <- Sys.glob(paste(seuratPath, '/markers_res*.rds', sep=""))

resolutions <- sub(".rds", '', sub(paste(seuratPath,'/markers_res', sep=""),'',paths))

markers <- vector("list")
for (res in resolutions){
  markers[[res]] <- readRDS(paste(seuratPath, "/markers_res", res, ".rds", sep=""))
}

if (genome == "mm10") {
  dbs <-c("Mouse_Gene_Atlas")
}else if(genome == "hg38"){
  dbs <-c("Human_Gene_Atlas")
}


setwd(workdir)

for (res in resolutions) {
  degs <- markers[[res]]

  clusters <- unique(degs$cluster)

  temp <- degs %>% group_by(cluster)

  result <- vector("list")
  for(i in clusters){
    result[i] <- list(temp$gene[temp$cluster==i])
  }

  dir.create(paste("res", res, sep=""))

  for (i in clusters) {
    interest <- result[[i]]
    enriched <- enrichr(interest[1:100], dbs)
    for (j in names(enriched)) {
      while('"expired": true' %in% enriched[[j]][1,]){
      	enriched <- enrichr(interest[1:100], dbs)
      }
      write.csv(enriched[[j]], file = paste("res", res, "/resolution", res, "_cluster", i, "_", j, "_top100.csv", sep=""))
    }

    enriched <- enrichr(interest, dbs)
    for (j in names(enriched)) {
      while('"expired": true' %in% enriched[[j]][1,]){
      	enriched <- enrichr(interest[1:100], dbs)
      }
      write.csv(enriched[[j]], file = paste("res", res, "/resolution", res, "_cluster", i, "_", j, ".csv", sep=""))
    }
  }
}
