library(scran)
library(irlba)
library(scater)
library(SingleCellExperiment)
library(Matrix)
library(scales)
library(BiocParallel)
library(BiocNeighbors)
library(batchelor)



getHVGs <- function(log_counts,min_mean = 0.001,FDR = 0.01)
{
  #identifying highly variable genes from a matrix of log counts
  # min_mean: 
  exclude <- c("ENSMUSG00000086503",read.table("../data/ygenes.tab", stringsAsFactors = FALSE)[,1],"tomato-td")
  stats <- modelGeneVar(log_counts,subset.row = setdiff(rownames(log_counts),exclude))
  stats <- stats[stats$mean > min_mean,]
  top_genes <- rownames(stats[stats$FDR < FDR,])
  return(top_genes)
} 



getmode <- function(v, dist) {
  # taken from
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}


mnnMap = function(atlas_pca, atlas_meta, map_pca, meta_chimera, k_map = 10, return.pca = FALSE,prop=NULL,BPPARAM=SerialParam()){
  # taken from (with minor adaptations)
  # https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
  
  correct = reducedMNN(atlas_pca, map_pca,prop.k=prop,BPPARAM=BPPARAM)$corrected
  atlas = 1:nrow(atlas_pca)
  correct_atlas = correct[atlas,]
  correct_map = correct[-atlas,]
  
  pre <- buildIndex(correct_atlas, BNPARAM=KmknnParam())
  
  knns = BiocNeighbors::queryKNN(correct_atlas, correct_map, k = 10, 
                                 get.index = TRUE, get.distance = TRUE,
                                 BPPARAM=BPPARAM,precomputed=pre)
  
  k.mapped <- matrix(0,nrow=nrow(knns$index),ncol=ncol(knns$index))
  for (j in 1:nrow(knns$index))
  {
    k.mapped[j,] <- as.vector(atlas_meta$cell[knns$index[j,]])
  }
  celltypes <-  matrix(0,nrow=nrow(knns$index),ncol=ncol(knns$index))
  stages <-  matrix(0,nrow=nrow(knns$index),ncol=ncol(knns$index))
  for (j in 1:nrow(knns$index))
  {
    xx <- match(k.mapped[j,], atlas_meta$cell)
    celltypes[j,] <- atlas_meta$celltype.clustering[xx]
    stages[j,] <- as.vector(atlas_meta$stage)[xx]
  }
  
  celltype.mapped <- rep("",nrow(celltypes))
  for (j in 1:nrow(celltypes)){
    celltype.mapped[j] <- getmode(celltypes[j,],1:ncol(celltypes))
  }
  stage.mapped <- rep("",nrow(stages))
  for (j in 1:nrow(stages)){
    stage.mapped[j] <- getmode(stages[j,],1:ncol(stages))
  }
  
  out = lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped = k.mapped[x,],
         celltype.mapped = celltype.mapped[x],
         stage.mapped = stage.mapped[x])
  })
  
  names(out) = meta_chimera$cell
  
  if(return.pca){
    return(list(out,list("atlas" = correct_atlas,
                         "mapped" = correct_map)))
  }
  
  
  return(out)
  
}




doBatchCorrect2 = function(counts, 
                           timepoints, 
                           samples, 
                           timepoint_order = c("E6.5", "E6.75", "E7.0", "mixed_gastrulation", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5","E8.75","E9.0","E9.25"), 
                           sample_order=NULL, 
                           npc = 50,
                           pc_override = NULL, 
                           BPPARAM = SerialParam()){
  #taken from 
  # https://github.com/MarioniLab/TChimeras2020/blob/master/core_functions/embryos_core_functions.R
  #default: samples largest to smallest
  if(is.null(sample_order)){
    tab = table(samples)
    sample_order = names(tab)[order(tab), decreasing = TRUE]
  }
  
  #remove timepoints that are not present
  timepoint_order = timepoint_order[timepoint_order %in% timepoints]
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list = lapply(unique(timepoints), function(tp){
    sub_pc = pca[timepoints == tp, , drop = FALSE]
    sub_samp = samples[timepoints == tp]
    list = lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) = unique(sub_samp)
    return(list)
  })
  
  names(pc_list) = unique(timepoints)
  
  #arrange to match timepoint order
  pc_list = pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list = lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list = lapply(pc_list, function(x){
    if(length(x) > 1){
      return(reducedMNN(x, BPPARAM = BPPARAM)$corrected)
    } else {
      return(x[[1]])
    }
  })
  #perform correction over list
  if(length(correct_list)>1){
    correct = reducedMNN(correct_list, BPPARAM = BPPARAM)$corrected
  } else {
    correct = correct_list[[1]]
  }
  
  correct = correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}



wrapper_map_function  <-  function(atlas_sce, atlas_meta, sce_chimera, meta_chimera, return.pca = FALSE, nPC = 50,
                      target_name,latent_time=NULL,BPPARAM=SerialParam(),saveFileName="temp",prop=NULL){
  # taken from the mapWrap function 
  # https://github.com/MarioniLab/TChimeras2020/blob/master/core_functions/embryos_core_functions.R
  
  #prevent duplicate rownames
  colnames(sce_chimera) = paste0("map_", colnames(sce_chimera))
  meta_chimera$cell = paste0("map_", meta_chimera$cell)
  sce_chimera$sizeFactor <- NULL
  hvgs = getHVGs(logcounts(atlas_sce),FDR=0.05)
  big_pca = multiBatchPCA(atlas_sce[hvgs,],sce_chimera[hvgs,], weights = c(1,0))
  
  names(big_pca) <- c("atlas","chimera")
  rownames(big_pca$atlas) = colnames(atlas_sce) 
  rownames(big_pca$chimera) = colnames(sce_chimera)
  pca_atlas = big_pca[[1]]
  map_pca = big_pca[[2]]
  
  #correct the atlas first
  order_df = atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells = sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage = factor(order_df$stage, 
                          levels = rev(c("E9.25","E9.0","E8.75","E8.5", 
                                         "E8.25", 
                                         "E8.0", 
                                         "E7.75", 
                                         "E7.5", 
                                         "E7.25", 
                                         "mixed_gastrulation", 
                                         "E7.0", 
                                         "E6.75", 
                                         "E6.5")))
  order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage = as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected = doBatchCorrect2(counts = logcounts(atlas_sce[hvgs,]), 
                                    timepoints = atlas_meta$stage, 
                                    samples = atlas_meta$sample, 
                                    timepoint_order = order_df$stage, 
                                    sample_order = order_df$sample, 
                                    pc_override = pca_atlas,
                                    BPPARAM=BPPARAM)
  
  
  
  mapping = mnnMap(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca = map_pca,
                   meta_chimera = meta_chimera,
                   return.pca = FALSE,prop=prop,BPPARAM=BPPARAM)
  saveRDS(extract_mapping(mapping),file=saveFileName)}

extract_mapping <- function(mapping){
  cn = substr(names(mapping), 5, nchar(names(mapping))) # remove
  ct = sapply(mapping, function(x) x$celltype.mapped)
  #ctc = sapply(mapping, function(x) x$celltype.mapped,clustering)
  st = sapply(mapping, function(x) x$stage.mapped)
  closest = sapply(mapping, function(x) x$cells.mapped[1])
  #tm = sapply(mapping, function(x) x$time.mapped)
  out = data.frame(cell = cn, celltype.mapped = ct, stage.mapped = st, closest.cell = closest )
  return(out)
  
}



