library(scran)
library(scater)
library(destiny)
Celltypes <- read.table("lineage_names.txt")$V1

chimeraT <- readRDS("../results/chimeraT_2020_extended_mapping_rev.rds")
chimeraWT <- readRDS("../results/chimeraWT_2020_extended_mapping_rev.rds")
chimeraMixl1 <- readRDS("../results/chimeraMixl1_2020_extended_mapping_rev.rds")

correlation_pseudotime_T <- list()
mnn_pseudotime_T <- list()

correlation_pseudotime_WT <- list()
mnn_pseudotime_WT <- list()

correlation_pseudotime_Mixl1 <- list()
mnn_pseudotime_Mixl1 <- list()

for (j in 1:length(Celltypes)){
  Celltype <- Celltypes[j]
  temp_file <- list.files(paste0("results_no_split/",Celltype),pattern="cor_cells_Mixl1")
  if (length(temp_file) > 0){
    cor_T <- readRDS(paste0("results_no_split/",Celltype,"/",Celltype,"_cor_cells_T.rds"))
    cor_Mixl1 <- readRDS(paste0("results_no_split/",Celltype,"/",Celltype,"_cor_cells_Mixl1.rds"))
    cor_WT <- readRDS(paste0("results_no_split/",Celltype,"/",Celltype,"_cor_cells_WT.rds"))
    pseudotime_files <- list.files(paste0("results_no_split/",Celltype),pattern="pseudotime_sublineage")
    pseudotime_files <- pseudotime_files[grepl(".rds",pseudotime_files)]
    pseudotimes <- readRDS(paste0("results_no_split/",Celltype,"/",pseudotime_files[1]))
    correlation_pseudotime_T[[j]] <- apply(cor_T, 2, function(x) mean(pseudotimes[rownames(cor_T)[order(x,decreasing=T) %in% 1:10]]))
    mnn_pseudotime_T[[j]] <- pseudotimes[chimeraT$closest.cell[match(colnames(cor_T),chimeraT$cell)]]
    correlation_pseudotime_WT[[j]] <- apply(cor_WT, 2, function(x) mean(pseudotimes[rownames(cor_WT)[order(x,decreasing=T) %in% 1:10]]))
    mnn_pseudotime_WT[[j]] <- pseudotimes[chimeraWT$closest.cell[match(colnames(cor_WT),chimeraWT$cell)]]
    correlation_pseudotime_Mixl1[[j]] <- apply(cor_Mixl1, 2, function(x) mean(pseudotimes[rownames(cor_Mixl1)[order(x,decreasing=T) %in% 1:10]]))
    mnn_pseudotime_Mixl1[[j]] <- pseudotimes[chimeraMixl1$closest.cell[match(colnames(cor_Mixl1),chimeraMixl1$cell)]]
  }
 
}

saveRDS(correlation_pseudotime_Mixl1,file="results_no_split/correlation_pseudotime_Mixl1.rds")
saveRDS(mnn_pseudotime_Mixl1,file="results_no_split/mnn_pseudotime_Mixl1.rds")
saveRDS(correlation_pseudotime_T,file="results_no_split/correlation_pseudotime_T.rds")
saveRDS(mnn_pseudotime_T,file="results_no_split/mnn_pseudotime_T.rds")
saveRDS(correlation_pseudotime_WT,file="results_no_split/correlation_pseudotime_WT.rds")
saveRDS(mnn_pseudotime_WT,file="results_no_split/mnn_pseudotime_WT.rds")

# Compare correlation and mnn based pseudotimes

correlation_pseudotime_T <- readRDS("results_no_split/correlation_pseudotime_T.rds")
mnn_pseudotime_T <- readRDS("results_no_split/mnn_pseudotime_T.rds")
names(correlation_pseudotime_T) <- Celltypes
names(mnn_pseudotime_T) <- Celltypes
xx <- unlist(lapply(correlation_pseudotime_T,function(x)!(is.null(x))))
correlation_pseudotime_T <- correlation_pseudotime_T[xx]
mnn_pseudotime_T <- mnn_pseudotime_T[xx]
correlation_pseudotimes <- sapply(1:length(correlation_pseudotime_T),function(j) cor(correlation_pseudotime_T[[j]],mnn_pseudotime_T[[j]]))
names(correlation_pseudotimes) <- names(correlation_pseudotime_T)
sort(correlation_pseudotimes)
write.table(sort(correlation_pseudotimes),file="correlation_between_pseudotimes_T.txt",col.names=F,quote=F)
