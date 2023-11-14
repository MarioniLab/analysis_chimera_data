library(scran)
library(scater)
library(destiny)
Celltypes <- read.table("lineage_names.txt")$V1

chimeraT <- readRDS("../results/chimeraT_2020_extended_mapping_rev.rds")
chimeraWT <- readRDS("../results/chimeraWT_2020_extended_mapping_rev.rds")
chimeraMixl1 <- readRDS("../results/chimeraMixl1_2020_extended_mapping_rev.rds")

correlation_pseudotime_T <- list()

correlation_pseudotime_WT <- list()

correlation_pseudotime_Mixl1 <- list()

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
    correlation_pseudotime_WT[[j]] <- apply(cor_WT, 2, function(x) mean(pseudotimes[rownames(cor_WT)[order(x,decreasing=T) %in% 1:10]]))
    correlation_pseudotime_Mixl1[[j]] <- apply(cor_Mixl1, 2, function(x) mean(pseudotimes[rownames(cor_Mixl1)[order(x,decreasing=T) %in% 1:10]]))
  }
 
}

saveRDS(correlation_pseudotime_Mixl1,file="results_no_split/correlation_pseudotime_Mixl1.rds")
saveRDS(correlation_pseudotime_T,file="results_no_split/correlation_pseudotime_T.rds")
saveRDS(correlation_pseudotime_WT,file="results_no_split/correlation_pseudotime_WT.rds")

