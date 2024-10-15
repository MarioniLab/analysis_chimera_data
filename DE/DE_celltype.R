# This script computes differential expression results for T and Mixl1 chimeras,
# for all cell types with at least 100 cells for both the WT and the respective target
# chimeras. 

## ----echo=FALSE,output=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(nebula)
chimeraWT <- readRDS("../results/chimeraWT_2020_extended_mapping_rev.rds")
chimeraMixl1 <- readRDS("../results/chimeraMixl1_2020_extended_mapping_rev.rds")
chimeraT <- readRDS("../results/chimeraT_2020_extended_mapping_rev.rds")

# First we read in the data and perform batch normalisation across the different chimeras.
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(scran)
library(scater)
library(batchelor)
celltypes_WT <- names(table(chimeraWT$celltype.mapped))[table(chimeraWT$celltype.mapped)>=100]
# DE only for those cell types with at least 100 cells for both target and WT chimeras
celltypes_Mixl1 <- names(table(chimeraMixl1$celltype.mapped))[table(chimeraMixl1$celltype.mapped)>=100]
celltypes_T <- names(table(chimeraT$celltype.mapped))[table(chimeraT$celltype.mapped)>=100]
genes_Mixl1 <- intersect(rownames(chimeraWT),rownames(chimeraMixl1))
genes_T <- intersect(rownames(chimeraWT),rownames(chimeraT))
coldata_intersect <- intersect(names(colData(chimeraT)),intersect(names(colData(chimeraWT)),
                                names(colData(chimeraMixl1))))
colData(chimeraMixl1) <- colData(chimeraMixl1)[,coldata_intersect]
colData(chimeraWT) <- colData(chimeraWT)[,coldata_intersect]
colData(chimeraT) <- colData(chimeraT)[,coldata_intersect]
reducedDims(chimeraMixl1) <- NULL
reducedDims(chimeraWT) <- NULL
reducedDims(chimeraT) <- NULL
logcounts(chimeraMixl1) <- NULL

# same of the barcodes are identical across the different data sets
# so we rename them before merging the data sets
colnames(chimeraWT) <- paste0(colnames(chimeraWT),"_WT")
colnames(chimeraT) <- paste0(colnames(chimeraT),"_T")
colnames(chimeraMixl1) <- paste0(colnames(chimeraMixl1),"_Mixl1")

sce_batch_normalised_Mixl1 <- multiBatchNorm(chimeraWT[genes_Mixl1,],chimeraMixl1[genes_Mixl1,])
sce_Mixl1 <- cbind(sce_batch_normalised_Mixl1[[1]], sce_batch_normalised_Mixl1[[2]])
sce_Mixl1$target <- c(rep("WT",ncol(sce_batch_normalised_Mixl1[[1]])),
                    rep("Mixl1",ncol(sce_batch_normalised_Mixl1[[2]])))

sce_Mixl1$sample[1:ncol(sce_batch_normalised_Mixl1[[1]])] <- 
  paste0(sce_Mixl1$sample[1:ncol(sce_batch_normalised_Mixl1[[1]])],"_WT")
sce_Mixl1$sample[-(1:ncol(sce_batch_normalised_Mixl1[[1]]))] <- 
  paste0(sce_Mixl1$sample[-(1:ncol(sce_batch_normalised_Mixl1[[1]]))],"_Mixl1")

sce_batch_normalised_T <- multiBatchNorm(chimeraWT[genes_T,],chimeraT[genes_T,])
sce_T <- cbind(sce_batch_normalised_T[[1]], sce_batch_normalised_T[[2]])
sce_T$target <- c(rep("WT",ncol(sce_batch_normalised_T[[1]])),
                      rep("T",ncol(sce_batch_normalised_T[[2]])))

sce_T$sample[1:ncol(sce_batch_normalised_T[[1]])] <- 
  paste0(sce_T$sample[1:ncol(sce_batch_normalised_T[[1]])],"_WT")
sce_T$sample[-(1:ncol(sce_batch_normalised_T[[1]]))] <- 
  paste0(sce_T$sample[-(1:ncol(sce_batch_normalised_T[[1]]))],"_T")

gene_conv <-  read.table("../data/genes.tsv")

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
celltypes_WT_Mixl1 <- intersect(celltypes_WT,celltypes_Mixl1)
celltypes_WT_T <- intersect(celltypes_WT,celltypes_T)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ----DE_Mixl1------------------------------------------------------------------------------------------------------------------------------------------------------------

markers <- list()
sce_Mixl1$target_tomato <- FALSE
sce_Mixl1$target_tomato[sce_Mixl1$tomato & sce_Mixl1$target=="Mixl1"] <- TRUE
celltypes_WT_Mixl1 <- intersect(celltypes_WT,celltypes_Mixl1)
sce_markers <- sce_Mixl1[,sce_Mixl1$celltype.mapped %in% celltypes_WT_Mixl1]
for (j in 1:length(celltypes_WT_Mixl1)){
  sce_temp <- sce_markers[,sce_markers$celltype.mapped == celltypes_WT_Mixl1[j]]
  sce_data <- scToNeb(obj = sce_temp, assay = "counts", id = "sample", pred = c("target","tomato","target_tomato"), offset="sizeFactor")
  df <- model.matrix(~target+tomato+target_tomato, data=sce_data$pred)
  
  markers_temp <- nebula(sce_data$count,sce_data$id,pred=df[,c("(Intercept)","targetWT","tomatoTRUE","target_tomatoTRUE")],offset=sce_data$offset)
  temp <- markers_temp$summary
  temp$DE_gene_name <- gene_conv$V2[match(temp$gene,gene_conv$V1)]
  temp$FDR <- p.adjust(temp$p_target_tomatoTRUE,method="BH")
  markers[[j]] <- temp
  save_file_name <- paste0(celltypes_WT_Mixl1[j],"_Mixl1_DE.csv")
  save_file_name <- gsub(" ","_",save_file_name)
  save_file_name <- gsub("/","_",save_file_name)
  temp <- temp[!(is.na(temp$FDR)),]
  temp <- temp[!(is.na(temp$FDR)),]
  write.table(temp,file=paste0("markers_DE_celltype/",save_file_name),sep=",",col.names = TRUE,row.names=FALSE)
}
markers <- lapply(markers,function(x) x[,c("DE_gene_name","gene","p_target_tomatoTRUE","FDR","logFC_target_tomatoTRUE")])
markers <- do.call(rbind,markers)
colnames(markers) <- c("DE_gene_name","DE_gene","p_value","FDR","logFC")

write.table(markers[markers$FDR<0.1,],file="markers_DE_celltype/Mixl1_DE.csv",sep=",",col.names = TRUE,row.names=FALSE)


file_name_from_CT <- function(CT){
  save_file_name <- paste0(CT,"_Mixl1_DE.csv")
  save_file_name <- gsub(" ","_",save_file_name)
  save_file_name <- gsub("/","_",save_file_name)
  return(save_file_name)
}
file_names_WT_Mixl1 <- intersect(sapply(celltypes_WT_Mixl1,file_name_from_CT),list.files("markers_DE_celltype"))
markers_Mixl1 <- lapply(file_names_WT_Mixl1,function(x) read.table(paste0("markers_DE_celltype/",x),sep=",",header=TRUE))
for (j in 1:length(markers_Mixl1)){
  markers_Mixl1[[j]]$cell_type <- strsplit(file_names_WT_Mixl1[j],"_Mixl1")[[1]][1]
}

markers_Mixl1 <- lapply(markers_Mixl1,function(x) x[,c("cell_type","DE_gene_name","gene","p_target_tomatoTRUE","FDR","logFC_target_tomatoTRUE")])
markers_Mixl1 <- do.call(rbind,markers_Mixl1)
colnames(markers_Mixl1) <- c("cell_type","DE_gene_name","DE_gene","p_value","FDR","logFC")

write.table(markers_Mixl1[markers_Mixl1$FDR<0.1,],file="markers_DE_celltype/Mixl1_DE.csv",sep=",",col.names = TRUE,row.names=FALSE)


## ----DE_T------------------------------------------------------------------------------------------------------------------------------------------------------------

markers <- list()
sce_T$target_tomato <- FALSE
sce_T$target_tomato[sce_T$tomato & sce_T$target=="T"] <- TRUE
celltypes_WT_T <- intersect(celltypes_WT,celltypes_T)
sce_markers <- sce_T[,sce_T$celltype.mapped %in% celltypes_WT_T]
for (j in 1:length(celltypes_WT_T)){
  sce_temp <- sce_markers[,sce_markers$celltype.mapped == celltypes_WT_T[j]]
  sce_data <- scToNeb(obj = sce_temp, assay = "counts", id = "sample", pred = c("target","tomato","target_tomato"), offset="sizeFactor")
  df <- model.matrix(~target+tomato+target_tomato, data=sce_data$pred)
 
  markers_temp <- nebula(sce_data$count,sce_data$id,pred=df[,c("(Intercept)","targetWT","tomatoTRUE","target_tomatoTRUE")],offset=sce_data$offset)
  temp <- markers_temp$summary
  temp$DE_gene_name <- gene_conv$V2[match(temp$gene,gene_conv$V1)]
  temp$FDR <- p.adjust(temp$p_target_tomatoTRUE,method="BH")
  markers[[j]] <- temp
  save_file_name <- paste0(celltypes_WT_T[j],"_T_DE.csv")
  save_file_name <- gsub(" ","_",save_file_name)
  save_file_name <- gsub("/","_",save_file_name)
  temp <- temp[!(is.na(temp$FDR)),]
  temp <- temp[!(is.na(temp$FDR)),]
  write.table(temp,file=paste0("markers_DE_celltype/",save_file_name),sep=",",col.names = TRUE,row.names=FALSE)
}
markers <- lapply(markers,function(x) x[,c("DE_gene_name","gene","p_target_tomatoTRUE","FDR","logFC_target_tomatoTRUE")])
markers <- do.call(rbind,markers)
colnames(markers) <- c("DE_gene_name","DE_gene","p_value","FDR","logFC")

write.table(markers[markers$FDR<0.1,],file="markers_DE_celltype/T_DE.csv",sep=",",col.names = TRUE,row.names=FALSE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_name_from_CT <- function(CT){
  save_file_name <- paste0(CT,"_T_DE.csv")
  save_file_name <- gsub(" ","_",save_file_name)
  save_file_name <- gsub("/","_",save_file_name)
  return(save_file_name)
}
file_names_WT_T <- intersect(sapply(celltypes_WT_T,file_name_from_CT),list.files("markers_DE_celltype"))
markers_T <- lapply(file_names_WT_T,function(x) read.table(paste0("markers_DE_celltype/",x),sep=",",header=TRUE))
for (j in 1:length(markers_T)){
  markers_T[[j]]$cell_type <- strsplit(file_names_WT_T[j],"_T")[[1]][1]
}

markers_T <- lapply(markers_T,function(x) x[,c("cell_type","DE_gene_name","gene","p_target_tomatoTRUE","FDR","logFC_target_tomatoTRUE")])
markers_T <- do.call(rbind,markers_T)
colnames(markers_T) <- c("cell_type","DE_gene_name","DE_gene","p_value","FDR","logFC")

write.table(markers_T[markers_T$FDR<0.1,],file="markers_DE_celltype/T_DE.csv",sep=",",col.names = TRUE,row.names=FALSE)

