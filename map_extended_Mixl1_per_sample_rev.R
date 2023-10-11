#like 2020, but PCA space based on atlas only
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
chimeraMixl1 <- readRDS("../results/Mixl1/chimeraMixl1.rds")
names(colData(chimeraMixl1))[names(colData(chimeraMixl1)) == "Sample"] <- "sample"
sce_chimera <- SingleCellExperiment(assays = list(counts=counts(chimeraMixl1),logcounts=logcounts(chimeraMixl1)))

meta_chimera <- colData(chimeraMixl1)

wrapper_map_function(atlas_sce, atlas_meta, sce_chimera[, meta_chimera$sample == strtoi(sample)], meta_chimera[meta_chimera$sample == strtoi(sample),],
           target_name="Mixl1",saveFileName=paste0("mapping_Mixl1_rev_sample_",sample,".rds"))


