library(Matrix)
library(scran)
library(scater)

atlas_sce <- readRDS("../data/big_atlas/big_atlas.rds")
atlas_meta = readRDS("../fromIvan_new_atlas/integrated_meta_celltype_clus.rds")


#Downsample to 10k cells per sample (or maximum number of cells in sample)
set.seed(42)
keep = lapply(unique(atlas_meta$stage), function(x){
  if(x %in% c("E9.5","mixed_gastrulation")){
    return(c())
  } else if(sum(atlas_meta$stage == x) < 10000) {
    return(which(atlas_meta$stage == x))
  } else {
    hits = which(atlas_meta$stage == x)
    return(sample(hits, 10000))
  }
})
keep = do.call(c, keep)

atlas_sce = atlas_sce[,keep]
atlas_meta = atlas_meta[keep,]
saveRDS(atlas_sce,file="../atlas_sce_sub.rds")
saveRDS(atlas_meta,file="../atlas_meta_sub.rds")