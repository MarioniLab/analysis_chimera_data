# This script assigns cells from the extended mouse gastrulation atlas to lineages, 
# by clustering WOT scores using mixture models of skewed t-distributions. 

Celltype <- commandArgs(trailingOnly=TRUE)[1]#the final cell type of the lineage
split_trajectory <- FALSE
#split_trajectory <- as.logical(strtoi(commandArgs(trailingOnly=TRUE)[2])) #an additional
# parameter used to split trajectories into sub-trajectories, not used here

set.seed(444)
setwd("/nfs/research/marioni/magda/chimera/scripts_for_paper")
knitr::opts_chunk$set(warning = FALSE, message = FALSE,cache=TRUE,cache.lazy = FALSE)
library(SingleCellExperiment)
library(cluster)
library(biomaRt)
library(patchwork)
source("chimera_core_functions_big_atlas.R")

if(split_trajectory == TRUE){
  dir.create(paste0(Celltype,"_split"))
}

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Reading in the meta data for the extended mouse gastrulation atlas (provided by I. Imaz-Rosshandler)
atlas_meta  <- readRDS("../fromIvan_new_atlas/integrated_meta_celltype_clus.rds")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Reading in WOT scores from the extended mouse gastrulation atlas (provided by I. Imaz-Rosshandler)
# Subsetting the matrix of WOT scores to the lineage considered
library(data.table)
wot_data <- fread("wot_output/trajectory_from_9.25.txt_trajectory.txt",header=TRUE)
wot_data <- wot_data[match(atlas_meta$cell[atlas_meta$stage %in% c("E9.25","E9.0","E8.75","E8.5","E8.25","E8.0","E7.75","E7.5")],wot_data$id),]

wot_scores_Celltype <- wot_data[[Celltype]]
names(wot_scores_Celltype) <- wot_data$id


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Creating a list of the WOT scores per stage.
# As a first step, we subset to the cells with the highest 5% of scores.
stages_sub <- c("E9.25","E9.0","E8.75","E8.5","E8.25","E8.0","E7.75","E7.5")
wot_scores_per_stage <- list()
wot_scores_per_stage[[1]] <- NULL
for (j in 2:length(stages_sub)){
  wot_scores_per_stage[[j]] <- wot_scores_Celltype[intersect(names(wot_scores_Celltype),atlas_meta$cell[atlas_meta$stage==stages_sub[j]])]
  index_temp <-  wot_scores_per_stage[[j]] > quantile(wot_scores_per_stage[[j]],0.95)
  wot_scores_per_stage[[j]]  <- wot_scores_per_stage[[j]][index_temp]
}
names(wot_scores_per_stage) <- stages_sub


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Identifying the cells that are part of the lineage by applying mixture modelling to the WOT scores (based on the highest 5% of WOT scores used above)
library(mixsmsn)
for (j in 2:length(stages_sub))
{
  fits <- list()
  for (k in 1:2){
    fits[[k]] <- smsn.mix(log2(wot_scores_per_stage[[j]]), nu = 3, g = k+1, get.init = TRUE, criteria = TRUE,
  group = TRUE, family = "Skew.t", calc.im=FALSE)
  }
  fit <- fits[[which.max(c(fits[[1]]$bic,fits[[2]]$bic))]]
  table(fit$group)
  wot_scores_per_stage[[j]] <- wot_scores_per_stage[[j]][fit$group == which.max(fit$mu)] 
}




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Discarding cell types that do not constitute at least 10% of the cells of the lineage at any stage
for (j in 2:length(stages_sub))
{ 
  temp <- table(atlas_meta$celltype.clustering[match(names( wot_scores_per_stage[[j]]), atlas_meta$cell)])
  celltypes_keep <- names(temp)[temp >= 0.1*sum(temp)]
  wot_scores_per_stage[[j]] <- wot_scores_per_stage[[j]][atlas_meta$celltype.clustering[match(names( wot_scores_per_stage[[j]]), atlas_meta$cell)]%in% celltypes_keep]
}
wot_scores_per_stage_save <- wot_scores_per_stage



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
celltypes_trajectory <- c()
for (j in 1:length(stages_sub))
{ 
  print(stages_sub[j])
  temp <- table(atlas_meta$celltype.clustering[match(names( wot_scores_per_stage[[j]]), atlas_meta$cell)])
  print(temp)
  celltypes_trajectory <- c(celltypes_trajectory,names(temp))}
celltypes_trajectory <- unique(celltypes_trajectory)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  From the time point onward from which the final cell type is the most frequent 
#cell type in that trajectory we exclude other cell types from the trajectory.
is_dominant <- rep(FALSE,length(stages_sub)-1)
for (j in 1:length(stages_sub))
{ 
  temp <- table(atlas_meta$celltype.clustering[match(names( wot_scores_per_stage[[j]]), atlas_meta$cell)])
  print(temp)
  is_dominant[j-1] <- names(temp)[which.max(temp)] == Celltype
}

cutoff <- min((1:length(is_dominant))[!(is_dominant)])-1
cutoff <- min(cutoff,5)

if (cutoff >= 1){
  for (j in 2:(1+cutoff)){
  wot_scores_per_stage[[j]] <- wot_scores_per_stage[[j]][intersect(names(wot_scores_per_stage[[j]]),atlas_meta$cell[atlas_meta$celltype.clustering == Celltype])]
}
}



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for (j in 1:length(stages_sub))
{ 
  print(stages_sub[j])
  temp <- table(atlas_meta$celltype.clustering[match(names( wot_scores_per_stage[[j]]), atlas_meta$cell)])
  print(temp)
}




## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Celltype_atlas <- atlas_meta$cell[atlas_meta$celltype.clustering == gsub("_"," ",Celltype)]



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
atlas_sce <- readRDS("../data/big_atlas/big_atlas.rds")
atlas_sce$celltype.clustering <- atlas_meta$celltype.clustering
atlas_sce$stage <- atlas_meta$stage


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
initial_celltypes <- atlas_meta$celltype.clustering[match(names(wot_scores_per_stage[["E7.5"]]),atlas_meta$cell)]
print(table(initial_celltypes))



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
initial_celltypes_keep <- names(table(initial_celltypes))[table(initial_celltypes) >= 0.2*length(initial_celltypes)]
initial_celltypes_keep <- intersect(initial_celltypes_keep,names(table(initial_celltypes))[table(initial_celltypes) >= 20])

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
starting_cells <- list()
for (j in 1:length(initial_celltypes_keep)){
  starting_cells[[j]] <- names(wot_scores_per_stage[["E7.5"]])[initial_celltypes%in%initial_celltypes_keep[j]]
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#the following is only relevant in case of subsplitting of trajectories
if (length(starting_cells) > 1 && split_trajectory){
  trans_7.5_7.75 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_7.5_7.75.rds")
  trans_7.5_7.75 <- trans_7.5_7.75[names(wot_scores_per_stage[[7]]),names(wot_scores_per_stage[[8]])]
  gc()
  
  trans_7.75_8.0 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_7.75_8.0.rds")
  trans_7.75_8.0 <- trans_7.75_8.0[names(wot_scores_per_stage[[6]]),names(wot_scores_per_stage[[7]])]
  gc()
  
  trans_8.0_8.25 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_8.0_8.25.rds")
  trans_8.0_8.25 <- trans_8.0_8.25[names(wot_scores_per_stage[[5]]),names(wot_scores_per_stage[[6]])]
  gc()
  
  trans_8.25_8.5 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_8.25_8.5.rds")
  trans_8.25_8.5 <- trans_8.25_8.5[names(wot_scores_per_stage[[4]]),names(wot_scores_per_stage[[5]])]
  gc()
  
  trans_8.5_8.75 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_8.5_8.75.rds")
  trans_8.5_8.75 <- trans_8.5_8.75[names(wot_scores_per_stage[[3]]),names(wot_scores_per_stage[[4]])]
  gc()
  
  trans_8.75_9.0 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_8.75_9.0.rds")
  trans_8.75_9.0  <- trans_8.75_9.0[names(wot_scores_per_stage[[2]]),names(wot_scores_per_stage[[3]])]
  gc()
  
  trans_9.0_9.25 <- readRDS("../fromIvan_new_atlas/embryo_mapped_collapsed_9.0_9.25.rds")
  trans_9.0_9.25 <- trans_9.0_9.25[intersect(rownames(trans_9.0_9.25),Celltype_atlas),names(wot_scores_per_stage[[2]])]
  gc()
  
  trans_matrix_list <- list(trans_9.0_9.25,trans_8.75_9.0,trans_8.5_8.75,trans_8.25_8.5,trans_8.0_8.25,trans_7.75_8.0,trans_7.5_7.75)
  names(trans_matrix_list) <- c("trans_9.0_9.25","trans_8.75_9.0","trans_8.5_8.75","trans_8.25_8.5","trans_8.0_8.25","trans_7.75_8.0","trans_7.5_7.75")
  
  
  subtrajectory_probabilities <- vector(mode = "list", length = length(starting_cells))
  for (j in 1:length(starting_cells)){
    subtrajectory_probabilities[[j]] <- vector(mode = "list", length = 6)
    names(subtrajectory_probabilities[[j]]) <- names(wot_scores_per_stage)[seq(7,2,by=-1)]
    subtrajectory_probabilities[[j]][[1]] <- trans_matrix_list[[7]][,starting_cells[[j]]] %*%     rep(1/length(starting_cells[[j]]),length(starting_cells[[j]]))
    subtrajectory_probabilities[[j]][[1]] <- subtrajectory_probabilities[[j]][[1]]/sum(subtrajectory_probabilities[[j]][[1]])
    for (k in 2:6){
      subtrajectory_probabilities[[j]][[k]] <- trans_matrix_list[[7-k+1]]%*% subtrajectory_probabilities[[j]][[k-1]]
      subtrajectory_probabilities[[j]][[k]] <- subtrajectory_probabilities[[j]][[k]]/sum(subtrajectory_probabilities[[j]][[k]])
    }
    names(subtrajectory_probabilities[[j]]) <- c("E7.75","E8.0","E8.25","E8.5","E8.75","E9.0")
  }
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  decision_trajectory <- function(p){
    return(sample(1:length(p),size=1,prob=p))
  }
  ind_trajectory <- c()
  temp <- c()
  for (j in 1:length(starting_cells)){
    temp <- cbind(temp,subtrajectory_probabilities[[j]]$E7.75)}
  colnames(temp) <- initial_celltypes_keep
  ind_trajectory <- c(ind_trajectory,apply(temp,1,decision_trajectory))
  
  
  temp <- c()
  for (j in 1:length(starting_cells)){
    temp <- cbind(temp,subtrajectory_probabilities[[j]]$E8.0)}
  colnames(temp) <- initial_celltypes_keep
  ind_trajectory <- c(ind_trajectory,apply(temp,1,decision_trajectory))
  
  temp <- c()
  for (j in 1:length(starting_cells)){
    temp <- cbind(temp,subtrajectory_probabilities[[j]]$E8.25)}
  colnames(temp) <- initial_celltypes_keep
  ind_trajectory <- c(ind_trajectory,apply(temp,1,decision_trajectory))
  
  
  temp <- c()
  for (j in 1:length(starting_cells)){
    temp <- cbind(temp,subtrajectory_probabilities[[j]]$E8.5)}
  colnames(temp) <- initial_celltypes_keep
  ind_trajectory <- c(ind_trajectory,apply(temp,1,decision_trajectory))
  
  temp <- c()
  for (j in 1:length(starting_cells)){
    temp <- cbind(temp,subtrajectory_probabilities[[j]]$E8.75)}
  colnames(temp) <- initial_celltypes_keep
  ind_trajectory <- c(ind_trajectory,apply(temp,1,decision_trajectory))
  
  
  temp <- c()
  for (j in 1:length(starting_cells)){
    temp <- cbind(temp,subtrajectory_probabilities[[j]]$E9.0)}
  colnames(temp) <- initial_celltypes_keep
  ind_trajectory <- c(ind_trajectory,apply(temp,1,decision_trajectory))
  
  
  ind_trajectory <- unlist(ind_trajectory)
  
}else{
  starting_cells[[1]] <- do.call(c,starting_cells)
  if (length(starting_cells) > 1){
    for (k in 2:length(starting_cells)){
      starting_cells[[k]] <- NULL
    }
  }
  all_cells <- unlist(sapply(wot_scores_per_stage[1:7],names))
  ind_trajectory <- rep(1,length(all_cells))
  names(ind_trajectory) <- all_cells
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
yGenes=read.table("../data/ygenes.tab")
#We exclude Xist and genes on the y chromosome
exclude <- c("ENSMUSG00000086503",yGenes)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
genes.go <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                   filters = 'go', values = list('GO:0007049'), mart = ensembl)
genes.go <- unique(genes.go$ensembl_gene_id)

library(GO.db)
library(org.Mm.eg.db)
results <- AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0007049"), columns = c('ENSEMBL'), keytype = "GOALL")
exclude <- unique(exclude,results$ENSEMBL)

saveRDS(exclude,file="exclude_genes_pseudotime.rds")

exclude <- readRDS("exclude_genes_pseudotime.rds")

#Computing and plotting diffusion maps

atlas_sub_Celltype <- atlas_sce[,match(unlist(sapply(wot_scores_per_stage,names)),colnames(atlas_sce))]
dpts <- list()
library(destiny)

dir.create(Celltype)
#setwd(Celltype)
pseudotime_vecs <- vector(mode = "list", length = length(starting_cells))
# preselect genes for pseudotime
# need genes that are differential across different time points
# We use genes correlated with actual time, then perform a linear regression of gene expression on time. 
# We then discard those genes whose residuals from this regression are correlated with batch, to avoid influence on batch effects on the diffusion map. 

  exclude1 <- read.table("../data/spatial_genes.txt")$V1 #exclude the impact of spatial genes on diffusion map
  initial_celltypes_keep_no_tab <- sapply(initial_celltypes_keep,function(x) gsub(" ","_",x))
for (j in 1:length(starting_cells)){
  saveRDS(c(starting_cells[[j]],names(ind_trajectory[ind_trajectory == j])),file=paste0(Celltype,"/",Celltype,"_cells_sublineage_",initial_celltypes_keep_no_tab[j],".rds"))
  sce_temp <- atlas_sub_Celltype[setdiff(rownames(atlas_sub_Celltype),unique(c(exclude1,exclude))),colnames(atlas_sub_Celltype) %in% c(starting_cells[[j]],names(ind_trajectory[ind_trajectory == j]))]
  time_vec <- as.double(sapply(as.vector(atlas_meta$stage[match(colnames(sce_temp),atlas_meta$cell)]),function(x) strsplit(x,"E")[[1]][2]))
    sce_temp <- sce_temp[,!(is.na(time_vec))]
  time_vec <- time_vec[!(is.na(time_vec))]
  lf_init_end <- rowMeans(logcounts(sce_temp)[,time_vec == max(time_vec)])- rowMeans(logcounts(sce_temp)[,time_vec == min(time_vec)])
  sce_temp <- sce_temp[abs(lf_init_end) > 0.5,]
  p_vals <- p.adjust(apply(logcounts(sce_temp),1,function(v) summary(lm(v ~ time_vec))$coefficients[2,4]))
  residuals <- apply(logcounts(sce_temp),1,function(v) lm(v ~ time_vec)$residuals)
  sce_temp <- sce_temp[p_vals<10^-10,]
  residuals <- residuals[,p_vals<10^-10]
  batch_vec <- atlas_meta$sample[match(colnames(sce_temp),atlas_meta$cell)]
  p_vals_res <- p.adjust(apply(residuals,2,function(v) summary(lm(v ~ batch_vec))$coefficients[2,4]))
  sce_temp <- sce_temp[p_vals_res > 0.2,]
  write.table(rownames(sce_temp),file=paste0(Celltype,"/",Celltype,"genes_sublineage",initial_celltypes_keep_no_tab[j],".txt"),row.names=F,col.names=F,quote=F)
  
  
  dm <- DiffusionMap(t(as.matrix(logcounts(sce_temp))),n_pcs=10)
  saveRDS(dm,file=paste0(Celltype,"/",Celltype,"_diffusion_map_sublineage_",
                          initial_celltypes_keep_no_tab[j],".rds"))
  xx <- sample(1:length(dm$DC1))
  tmp <- data.frame(DC1 = eigenvectors(dm)[xx, 1],
                  DC2 = eigenvectors(dm)[xx, 2],
                  celltype = atlas_meta$celltype.clustering[match(colnames(sce_temp)[xx],atlas_meta$cell)],
                  stage = atlas_meta$stage[match(colnames(sce_temp)[xx],atlas_meta$cell)])

  p1 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = celltype)) +
    geom_point(size=0.2) + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic(base_size=11) + theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
      labs(color="")+ guides(color = guide_legend(nrow=2,override.aes = list(size = 3))) + ggtitle(paste0(Celltype,", ",
                    initial_celltypes_keep_no_tab[j]," sublineage"))
  
 
  
  ggsave(p1,file=paste0(Celltype,"/pseudotime_celltype_sublineage_",initial_celltypes_keep_no_tab[j],".pdf"),height=4,width=7)
  p2 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = stage)) +
    geom_point(size=1) + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic(base_size=11) + theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
      labs(color="")+ guides(color = guide_legend(override.aes = list(size = 3)))+ ggtitle(paste0(Celltype,", ",
                                                                                                  initial_celltypes_keep_no_tab[j]," sublineage"))
  ggsave(p2,file=paste0(Celltype,"/pseudotime_stage_sublineage_",initial_celltypes_keep_no_tab[j],".pdf"),height=4,width=7)

  #saveRDS(dm,file=paste0("diffusion_map_",initial_celltypes_keep_no_tab[j],".rds"))
  cors <- cor(dm@eigenvectors,time_vec)
  dpt <- dm@eigenvectors[,which.max(abs(cors))]
  if (cors[which.max(abs(cors))] < 0 ){
   dpt <- -dpt 
  }
  tmp <- data.frame(DC1=dm$DC1[xx],DC2=dm$DC2[xx],dpt=dpt[xx],celltype = atlas_meta$celltype.clustering[match(colnames(sce_temp)[xx],atlas_meta$cell)])
  p3 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = dpt)) +
    geom_point(alpha=0.5,size=2) + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic(base_size=11) +scale_color_viridis_c()+ theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
      labs(color="")+ guides(color = guide_legend(override.aes = list(size = 3)))+ ggtitle(paste0(Celltype,", ",
                                                                                                  initial_celltypes_keep_no_tab[j]," sublineage"))
  ggsave(p3,file=paste0(Celltype,"/pseudotime_sublineage_",initial_celltypes_keep_no_tab[j],".pdf"),height=4,width=7)
  pseudotime_vecs[[j]] <- dpt
  
  
  names(pseudotime_vecs[[j]]) <- colnames(sce_temp)
  saveRDS(pseudotime_vecs[[j]],file=paste0(Celltype,"/",Celltype,"_pseudotime_sublineage_",
                                initial_celltypes_keep_no_tab[j],".rds"))
  write.table(dm@eigenvectors[,1:2],file=paste0(Celltype,"/diffusion_components.txt"))
     #Plot intervals for different cell types
   p4 <- ggplot(tmp,aes(x = dpt,color=celltype)  ) + geom_density()+ ggtitle(paste0(Celltype,", ",
                                                                                    initial_celltypes_keep_no_tab[j]," sublineage"))
   ggsave(p4,file=paste0(Celltype,"/pseudotime_celltype_densities_sublineage_",initial_celltypes_keep_no_tab[j],".pdf"),height=4,width=7)

}

