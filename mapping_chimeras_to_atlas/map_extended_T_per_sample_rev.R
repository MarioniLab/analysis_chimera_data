#This scripts maps each sample from the T chimera to the atlas separately, using the mutual nearest neighbours method (Haghverdi , with principal components based on the reference atlas only, and not the chimera data set, in order to map the chimera cells based on the underlying structure of variation of the reference data set.

sample <-  commandArgs(trailingOnly=TRUE)
library(Matrix)
library(scran)
library(scater)
library(knitr)
library(pheatmap)
library(reshape2)
library(edgeR)
library(BiocNeighbors)
library(BiocParallel)
ncores = 20
mcparam = SnowParam(workers = ncores)
register(mcparam)

source("chimera_core_functions_mapping_2020_revised.R")


atlas_sce <- readRDS("../atlas_sce_sub.rds")
atlas_meta <- readRDS("../atlas_meta_sub.rds")
chimeraT <- readRDS("chimeraT_original_mapping.rds")
names(colData(chimeraT))[names(colData(chimeraT)) == "Sample"] <- "sample"
sce_chimera <- SingleCellExperiment(assays = list(counts=counts(chimeraT[,complete.cases(reducedDims(chimeraT)$pca.corrected.E8.5)])))
sce_chimera <- logNormCounts(sce_chimera,size.factors=chimeraT$sizeFactor[complete.cases(reducedDims(chimeraT)$pca.corrected.E8.5)])

meta_chimera <- colData(chimeraT[,complete.cases(reducedDims(chimeraT)$pca.corrected.E8.5)])

wrapper_map_function(atlas_sce, atlas_meta, sce_chimera[-nrow(sce_chimera), meta_chimera$sample == strtoi(sample)], meta_chimera[meta_chimera$sample == strtoi(sample),],
           target_name="T",saveFileName=paste0("mapping_T_rev_sample_",sample,".rds"))


