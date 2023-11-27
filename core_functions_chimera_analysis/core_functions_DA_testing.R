# This script implements perturbSuite_DA for chimera data sets. 
# It contains the same functions as core_functions_DA_testing_general.R,
# but with input parameters adapted to the chimera data, and an additional 
# function to test for differential abundance per mapped cell type and stage

library(SingleCellExperiment)
library(cluster)
library(biomaRt)
library(pROC)
library(patchwork)
library(scran)
library(scater)
source("chimera_core_functions_big_atlas.R")

#' Performs perturbSuite_DA_celltype
#' Computes and plots differential abundance per mapped cell type for perturbation experiments
#' with control experiment involving several cell types, using the normalisation step 
#' included in perturbSuite_DA. 
#' @title da_per_celltype
#' @param sce_chimera SingleCellExperiment for the target chimera dataset, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker indicating the knockout, and 'celltype.mapped'
#' @param chimera_WT SingleCellExperiment for the control experiment, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker and 'celltype.mapped'
#' @param target string, name of the gene that was the target of the knockout
#' @param alpha significance level, default is 10%
#' @return data.frame with the following columns: celltype, odds_ratio-odds ratio of the percentage of 
#' tomato positive cells among the target chimera cells of the celltype versus the respective ratio for the control experiment, 
#' p_values, sig- whether a celltype is significantly enriched, depleted, or not enriched,
#' @export
#' @examples 


da_per_celltype <- function(sce_chimera,chimera_WT, target,alpha=0.1){
  sce_chimera_tomato <- sce_chimera[,sce_chimera$tomato==TRUE]
  sce_chimera_tomato_neg <- sce_chimera[,sce_chimera$tomato==FALSE]
  chimera_WT_tomato_neg <- chimera_WT[,!(as.logical(chimera_WT$tomato))]
  chimera_WT_tomato <- chimera_WT[,as.logical(chimera_WT$tomato)]
  celltypes <- names(table(chimera_WT_tomato$celltype.mapped))[table(chimera_WT_tomato$celltype.mapped) >= 30]
  ratio_WT <- rep(0,length(celltypes))
  ratio_target <- rep(0, length(celltypes))
  for (j in 1:length(celltypes)){
    ratio_WT[j] <- sum(chimera_WT_tomato$celltype.mapped == celltypes[j])/sum(chimera_WT_tomato_neg$celltype.mapped == celltypes[j])
    ratio_target[j] <- sum(sce_chimera_tomato$celltype.mapped == celltypes[j])/sum(sce_chimera_tomato_neg$celltype.mapped == celltypes[j])
  } 
  names(ratio_target) <- celltypes
  names(ratio_WT) <- celltypes
  norm_factor_WT <- median(ratio_WT[(!(is.infinite(ratio_WT)))&(!(is.na(ratio_WT)))])
  norm_factor_target <- median(ratio_target[(!(is.infinite(ratio_target)))&(!(is.na(ratio_target)))]) 
  sample_number_target <- floor(min(ncol(sce_chimera_tomato),ncol(sce_chimera_tomato_neg)/norm_factor_target))
  sample_number_WT <- floor(min(ncol(chimera_WT_tomato),ncol(chimera_WT_tomato_neg)/norm_factor_WT))
  p_values <- list()
  odds_ratio <- list()
  for (k in 1:100){
    sce_chimera_tomato <- sce_chimera_tomato[,sample(1:ncol(sce_chimera_tomato),sample_number_target)]
    sce_chimera_tomato_neg <- sce_chimera_tomato_neg[,sample(1:ncol(sce_chimera_tomato_neg),floor(sample_number_target*norm_factor_target))]
    
    chimera_WT_tomato <- chimera_WT_tomato[,sample(1:ncol(chimera_WT_tomato),sample_number_WT)]
    chimera_WT_tomato_neg <- chimera_WT_tomato_neg[,sample(1:ncol(chimera_WT_tomato_neg),floor(sample_number_WT*norm_factor_WT))]
    
    odds_ratio[[k]] <- rep(0,length(celltypes))
    p_values[[k]] <- rep(0,length(celltypes))
    for (j in 1:length(celltypes)){
      table_fisher <- c(sum(sce_chimera_tomato$celltype.mapped == celltypes[j]),sum(sce_chimera_tomato_neg$celltype.mapped == celltypes[j]),
                        sum(chimera_WT_tomato$celltype.mapped == celltypes[j]),sum(chimera_WT_tomato_neg$celltype.mapped == celltypes[j]))
      dim(table_fisher) <- c(2,2)
      aa <- fisher.test(table_fisher)
      odds_ratio[[k]][j] <- aa$estimate
      p_values[[k]][j] <- aa$p.value
    } 
  }
  p_values <- do.call(cbind,p_values)
  odds_ratio <- do.call(cbind,odds_ratio)
  p_values <- apply(p_values,1,median)
  odds_ratio <- apply(odds_ratio,1,median)
  p_values <- p.adjust(p_values)
  fisher_test_celltypes <- data.frame(celltype=celltypes,p_values=p_values,odds_ratio=odds_ratio)
  fisher_test_celltypes <- fisher_test_celltypes[order(fisher_test_celltypes$p_values),]
  fisher_test_celltypes$sig <- "enriched"
  fisher_test_celltypes$sig[ fisher_test_celltypes$odds_ratio < 1] <- "depleted"
  fisher_test_celltypes$sig[fisher_test_celltypes$p_values > alpha] <- "not significant"

  p <- ggplot( fisher_test_celltypes , aes(x=log10(odds_ratio), y=-log10(p_values+1e-100),text=celltype)) +
    labs(x=expression(paste("log"[10]," of odds ratio")),y=expression(paste("-log"[10]," of FDR adjusted p-values"))) +
    geom_point(aes(color=sig),size=3) +
    scale_color_manual(values=c("not significant"="grey","enriched"="darkblue","depleted" = "darkred"),name="") +
    ggrepel::geom_text_repel(data=fisher_test_celltypes[(apply(cbind(fisher_test_celltypes$odds_ratio,1/fisher_test_celltypes$odds_ratio),1,max)>1.5) | (fisher_test_celltypes$p_values < 0.1),], 
                             aes(x=log10(odds_ratio), y=-log10(p_values+1e-100), label=celltype), max.overlaps=Inf, size=4) +
    theme_classic(base_size=14) + scale_y_continuous(trans=pseudo_log_trans(sigma=1,base=10))+
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black')) + ggtitle(paste0(target,"- per celltype enrichment"))+
    annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "blue",alpha=0.1) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "red",alpha=0.1) 
  print(p)
  ggsave(p,file=paste0("da_celltypes_",target,".pdf"))
  return(fisher_test_celltypes)
}


#' Computes and plots differential abundance per mapped cell type for perturbation experiments
#' with control experiment involving several cell types, WITHOUT perturbSuite_DA normalisation.
#' @title da_per_celltype_without_normalising
#' @param sce_chimera SingleCellExperiment for the target chimera dataset, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker indicating the knockout, and 'celltype.mapped'
#' @param chimera_WT SingleCellExperiment for the control experiment, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker and 'celltype.mapped'
#' @param target string, name of the gene that was the target of the knockout
#' @param alpha significance level, default is 10%
#' @return data.frame with the following columns: celltype, odds_ratio-odds ratio of the percentage of 
#' tomato positive cells among the target chimera cells of the celltype versus the respective ratio for the control experiment, 
#' p_values, sig- whether a celltype is significantly enriched, depleted, or not enriched,
#' @export
#' @examples 

da_per_celltype_without_normalising <- function(sce_chimera,chimera_WT, target,alpha=0.1){
  sce_chimera_tomato <- sce_chimera[,sce_chimera$tomato==TRUE]
  sce_chimera_tomato_neg <- sce_chimera[,sce_chimera$tomato==FALSE]
  chimera_WT_tomato_neg <- chimera_WT[,!(as.logical(chimera_WT$tomato))]
  chimera_WT_tomato <- chimera_WT[,as.logical(chimera_WT$tomato)]
  celltypes <- names(table(chimera_WT_tomato$celltype.mapped))[table(chimera_WT_tomato$celltype.mapped) >= 10]
  #celltypes <- unique(chimera_WT_tomato$celltype.mapped)
  
  odds_ratio <- rep(0,length(celltypes))
  p_values <- rep(0,length(celltypes))
  for (j in 1:length(celltypes)){
    table_fisher <- c(sum(sce_chimera_tomato$celltype.mapped == celltypes[j]),sum(sce_chimera_tomato_neg$celltype.mapped == celltypes[j]),
                      sum(chimera_WT_tomato$celltype.mapped == celltypes[j]),sum(chimera_WT_tomato_neg$celltype.mapped == celltypes[j]))
    dim(table_fisher) <- c(2,2)
    aa <- fisher.test(table_fisher)
    p_values[j] <-aa$p.value
    odds_ratio[j] <- aa$estimate
  } 
  p_values <- p.adjust(p_values)
  fisher_test_celltypes <- data.frame(celltype=celltypes,p_values=p_values,odds_ratio=odds_ratio)
  fisher_test_celltypes <- fisher_test_celltypes[order(fisher_test_celltypes$p_values),]
  fisher_test_celltypes$sig <- "enriched"
  fisher_test_celltypes$sig[ fisher_test_celltypes$odds_ratio < 1] <- "depleted"
  fisher_test_celltypes$sig[fisher_test_celltypes$p_values > alpha] <- "not significant"
  
  p <- ggplot( fisher_test_celltypes , aes(x=log10(odds_ratio), y=-log10(p_values+1e-100),text=celltype)) +
    labs(x=expression(paste("log"[10]," of odds ratio")),y=expression(paste("-log"[10]," of FDR adjusted p-values"))) +
    geom_point(aes(color=sig),size=3) +
    scale_color_manual(values=c("not significant"="grey","enriched"="darkblue","depleted" = "darkred"),name="") +
    ggrepel::geom_text_repel(data=fisher_test_celltypes[fisher_test_celltypes$sig!="not significant",], 
                             aes(x=log10(odds_ratio), y=-log10(p_values+1e-100), label=celltype), max.overlaps=Inf, size=4) +
    theme_classic(base_size=14) + scale_y_continuous(trans=pseudo_log_trans(sigma=1,base=10))+
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black')) + ggtitle(paste0(target,"- per celltype enrichment"))
  print(p)
  ggsave(p,file=paste0("da_celltypes_",target,".pdf"))
  return(fisher_test_celltypes)
}

#' Computes and plots differential abundance per mapped cell type and stage for perturbation experiments
#' with control experiment involving several cell types.  
#' @title da_per_celltype_per_stage
#' @param sce_chimera SingleCellExperiment for the target chimera dataset, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker indicating the knockout, and 'celltype.mapped'
#' @param chimera_WT SingleCellExperiment for the control experiment, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker and 'celltype.mapped'
#' @param target string, name of the gene that was the target of the knockout
#' @param alpha significance level, default is 10%
#' @param stages the mapped stages for which to do differential abundance testing
#' @param subset_plot by default NULL, if not NULL then it is the vector of names of the cell types to which to subset the plot
#' @return matrix of size number of celltypes x number of stages of odd ratios of the percentage of 
#' tomato positive cells among the target chimera cells of the celltype versus the respective ratio for the control experiment
#' @export
#' @examples 

da_per_celltype_per_stage <- function(sce_chimera,chimera_WT, target,alpha=0.1,stages = c("E7.75","E8.0","E8.25","E8.5","E8.75"),
                  subset_plot=NULL){
  sce_chimera_tomato <- sce_chimera[,as.logical(sce_chimera$tomato)]
  sce_chimera_tomato_neg <- sce_chimera[,!(as.logical(sce_chimera$tomato))]
  chimera_WT_tomato_neg <- chimera_WT[,!(as.logical(chimera_WT$tomato))]
  chimera_WT_tomato <- chimera_WT[,as.logical(chimera_WT$tomato)]
  celltypes <- names(table(chimera_WT_tomato$celltype.mapped))[table(chimera_WT_tomato$celltype.mapped) >= 30]
  ratio_WT <- matrix(NA,nrow=length(celltypes),ncol=length(stages))
  ratio_target <- matrix(NA, nrow=length(celltypes),ncol=length(stages))
  for (k in 1:length(stages)){
    for (j in 1:length(celltypes)){
      nr_WT_stage_ct <- sum(chimera_WT_tomato$celltype.mapped == celltypes[j]&chimera_WT_tomato$stage.mapped==stages[k])
      if (nr_WT_stage_ct >= 30){
        ratio_target[j,k] <- sum(sce_chimera_tomato$celltype.mapped == celltypes[j]&sce_chimera_tomato$stage.mapped==stages[k])/sum(sce_chimera_tomato_neg$celltype.mapped == celltypes[j]&sce_chimera_tomato_neg$stage.mapped==stages[k])
        ratio_WT[j,k] <- sum(chimera_WT_tomato$celltype.mapped == celltypes[j]&chimera_WT_tomato$stage.mapped==stages[k])/sum(chimera_WT_tomato_neg$celltype.mapped == celltypes[j]&chimera_WT_tomato_neg$stage.mapped==stages[k])
      }
    } 
  }
  rownames(ratio_target) <- celltypes
  rownames(ratio_WT) <- celltypes
  colnames(ratio_target) <- stages
  colnames(ratio_WT) <- stages
  
  norm_factor_target <- apply(ratio_target,2,function(x) median(x[(!(is.na(x))) & (!is.infinite(x))]))
  norm_factor_WT <- apply(ratio_WT,2,function(x) median(x[(!(is.na(x))) & (!is.infinite(x))]))
  
  p_values <- list()
  odds_ratio <- list()
  
  for (ij in 1:30){
    odds_ratio[[ij]] <- matrix(NA,nrow=length(celltypes),ncol=length(stages))
    p_values[[ij]] <- matrix(1,nrow=length(celltypes),ncol=length(stages))
    for (k in 1:length(stages)){
      sce_chimera_tomato_stage <- sce_chimera_tomato[,sce_chimera_tomato$stage.mapped == stages[k]]
      sce_chimera_tomato_neg_stage <- sce_chimera_tomato_neg[,sce_chimera_tomato_neg$stage.mapped == stages[k]]
      sample_number_target <- floor(min(sum(sce_chimera_tomato$stage.mapped==stages[k]),sum(sce_chimera_tomato_neg$stage.mapped==stages[k])/norm_factor_target[k]))
      sce_chimera_tomato_stage <- sce_chimera_tomato_stage[,sample(1:ncol(sce_chimera_tomato_stage),sample_number_target)]
      sce_chimera_tomato_neg_stage <- sce_chimera_tomato_neg_stage[,sample(1:ncol(sce_chimera_tomato_neg_stage),floor(sample_number_target*norm_factor_target[k]))]
      
      chimera_WT_tomato_stage <- chimera_WT_tomato[,chimera_WT_tomato$stage.mapped == stages[k]]
      chimera_WT_tomato_neg_stage <- chimera_WT_tomato_neg[,chimera_WT_tomato_neg$stage.mapped == stages[k]]
      sample_number_WT <- floor(min(sum(chimera_WT_tomato$stage.mapped==stages[k]),sum(chimera_WT_tomato_neg$stage.mapped==stages[k])/norm_factor_WT[k]))
      chimera_WT_tomato_stage <- chimera_WT_tomato_stage[,sample(1:ncol(chimera_WT_tomato_stage),sample_number_WT)]
      chimera_WT_tomato_neg_stage <- chimera_WT_tomato_neg_stage[,sample(1:ncol(chimera_WT_tomato_neg_stage),floor(sample_number_WT*norm_factor_WT[k]))]
      
      for (j in 1:length(celltypes)){
        table_fisher <- c(sum(sce_chimera_tomato_stage$celltype.mapped == celltypes[j]),sum(sce_chimera_tomato_neg_stage$celltype.mapped == celltypes[j]),
                          sum(chimera_WT_tomato_stage$celltype.mapped == celltypes[j]),sum(chimera_WT_tomato_neg_stage$celltype.mapped == celltypes[j]))
        dim(table_fisher) <- c(2,2)
        #if less than 30 cells are assigned to the target chimera (tomato+ and tomato-), or less than 30 cells are assigned to the WT chimera (tomato+ or tomato-), then we do not perform the test, 
        #but set the odds_ratio to NA, to avoid giving too much impact to potential misassignment of celltype
        rowSumsMatrix <- rowSums(table_fisher)
        #if (all(rowSumsMatrix >= 30)){
        aa <- fisher.test(table_fisher)
        p_values[[ij]][j,k] <- aa$p.value
        odds_ratio[[ij]][j,k] <- aa$estimate
        # }
      } 
    }
  }
  odds_ratio <- apply(simplify2array(odds_ratio), 1:2, median)
  p_values <- apply(simplify2array(p_values), 1:2, median)
  p_values <- apply(p_values,1,p.adjust)
  p_values <- apply(p_values,2,p.adjust)
  colnames(odds_ratio) <- stages
  rownames(odds_ratio) <- celltypes
  p_values <- t(p_values)
  colnames(p_values) <- stages
  rownames(p_values) <- celltypes
  odds_ratio <- odds_ratio[apply(odds_ratio,1,function(x) any(!(is.na(x)))),]
  odds_ratio_temp <- odds_ratio
  odds_ratio_temp[odds_ratio_temp==Inf] <- max(odds_ratio_temp[odds_ratio_temp!=Inf])
  a <- heatmap(odds_ratio_temp)
  odds_ratio <- odds_ratio[a$rowInd,]
  p_values <- p_values[a$rowInd,]
  df <- reshape2::melt(odds_ratio)
  df2 <- reshape2::melt(p_values)
  df$p_value <- df2$value
  colnames(df) <- c("celltype_mapped","stage_mapped","odds_ratio","p_value")
  df$sig <- df$p_value < 0.1
  df <- df[!(is.na(df$odds_ratio)),]
  df$log_odds_ratio <- log2(df$odds_ratio+0.0000001)
  df$log_odds_ratio<- cut(df$log_odds_ratio, c(Inf, c(3,1.5,0,-1.5,-3), -Inf))
  order <- rownames(odds_ratio)[order(apply(odds_ratio,1,function(x) max(x[!(is.na(x))])))]
  df$celltype_mapped <- factor(df$celltype_mapped,levels=order)
  if (!is.null(subset_plot)){
    df <- df[df$celltype_mapped%in%subset_plot,]
  }
  pp <- ggplot(df,mapping = aes(y = celltype_mapped, x = stage_mapped)) + 
    geom_raster(aes(fill=log_odds_ratio)) + 
    labs(y="mapped celltype", x="mapped stage", title=paste0(target,"- mapped celltypes")) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11)) + labs(fill="log2 of\n odds ratio")+
                       theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11)) + labs(fill="log2 of\n odds ratio")+
                       scale_fill_brewer(palette="RdYlBu") + geom_point(aes(size=ifelse(sig, "significant", "not_significant"))) +
    scale_size_manual(values=c(significant=2, not_significant=NA),guide="none")
    
  print(pp)
  
  return(df)
}



#' Computes and plots differential fate probability for perturbation experiments
#' with control experiment involving several trajectories.  
#'  
#' @title differential_fate_probability
#' @param sce_chimera SingleCellExperiment for the target chimera dataset, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker indicating the knockout, and 'celltype.mapped'
#' @param chimera_WT SingleCellExperiment for the control experiment, 
#' colData needs to include the slots 'tomato' (TRUE/FALSE) for the presence
#' of the fluorescent marker and 'celltype.mapped'
#' @param target string, name of the gene that was the target of the knockout
#' @param trajectory_scores dataframe whose first column is the cell-id, and whose remaining columns
#' are scores for the respective trajectory for each cell. 
#' @param alpha significance level, default is 10%
#' @return data.frame with the following columns: trajectory, odds_ratio-odds ratio of the percentage of 
#' tomato positive cells among the target chimera cells of the celltype versus the respective ratio for the control experiment, 
#' p_values, sig- whether a celltype is significantly enriched, depleted, or not enriched,
#' @export
#' @examples 


differential_fate_probability <- function(sce_chimera,chimera_WT,target,trajectory_scores,alpha=0.1){
  trajectories <- colnames(trajectory_scores)
  trajectory_scores$id <- rownames(trajectory_scores)
  sce_chimera_tomato <- sce_chimera[,as.logical(sce_chimera$tomato)]
  sce_chimera_tomato_neg <- sce_chimera[,!(as.logical(sce_chimera$tomato))]
  chimera_WT_tomato_neg <- chimera_WT[,!(as.logical(chimera_WT$tomato))]
  chimera_WT_tomato <- chimera_WT[,as.logical(chimera_WT$tomato)]
  
  quantiles <- rep(NA,length(trajectories))
  ratio_WT <- rep(NA,length(trajectories))
  ratio_target <-rep(NA,length(trajectories))
  
  names(ratio_target) <- trajectories
  names(ratio_WT) <- trajectories
  trajectories <- as.matrix(trajectories)
  
  sce_chimera_score_tomato <- trajectory_scores[trajectory_scores$id %in% sce_chimera_tomato$closest.cell,-ncol(trajectory_scores)]
  sce_chimera_score_tomato <- sce_chimera_score_tomato[apply(sce_chimera_score_tomato,1,function(x) max(x) > 10^-4),]
  
  
  sce_chimera_score_tomato_neg <- trajectory_scores[trajectory_scores$id %in% sce_chimera_tomato_neg$closest.cell,-ncol(trajectory_scores)]
  sce_chimera_score_tomato_neg <- sce_chimera_score_tomato_neg[apply(sce_chimera_score_tomato_neg,1,function(x) max(x) > 10^-4),]
  
  
  chimera_WT_score <- trajectory_scores[trajectory_scores$id %in% chimera_WT_tomato$closest.cell,-ncol(trajectory_scores)]
  chimera_WT_score  <- chimera_WT_score[apply(chimera_WT_score ,1,function(x) max(x) > 10^-4),] 
  
  chimera_WT_score_neg <- trajectory_scores[trajectory_scores$id %in% chimera_WT_tomato_neg$closest.cell,-ncol(trajectory_scores)] 
  chimera_WT_score_neg  <- chimera_WT_score_neg[apply(chimera_WT_score_neg ,1,function(x) max(x) > 10^-4),] 
  
  odds_ratio_all <- list()
  p_values_all <- list()
  
  for (jk in 1:30){
    sce_chimera_score_max <- apply(sce_chimera_score_tomato,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #sce_chimera_score_max <- apply(sce_chimera_score_tomato,1,function(x) return(names(x)[which.max(x)]))
    sce_chimera_score_max_table <- table(sce_chimera_score_max)
    sce_chimera_score_neg_max <- apply(sce_chimera_score_tomato_neg,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #sce_chimera_score_neg_max <- apply(sce_chimera_score_tomato_neg,1,function(x) return(names(x)[which.max(x)]))
    sce_chimera_score_neg_max_table <- table(sce_chimera_score_neg_max)
    chimera_WT_score_max <- apply(chimera_WT_score,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #chimera_WT_score_max <- apply(chimera_WT_score,1,function(x) return(names(x)[which.max(x)]))
    chimera_WT_score_max_table <- table(chimera_WT_score_max)
    chimera_WT_score_neg_max <- apply(chimera_WT_score_neg,1,function(x) {y <- x>10^-4; temp <- x[y];
    return(sample(names(temp),1,prob=temp/sum(temp)))})
    #chimera_WT_score_neg_max <- apply(chimera_WT_score_neg,1,function(x) return(names(x)[which.max(x)]))
    chimera_WT_score_neg_max_table <- table(chimera_WT_score_neg_max)
    
    
    nr_target_tomato <- rep(0, length(trajectories))
    names(nr_target_tomato) <- trajectories
    nr_target_tomato[names(sce_chimera_score_max_table)] <- sce_chimera_score_max_table
    
    nr_target_tomato_neg <- rep(0, length(trajectories))
    names(nr_target_tomato_neg) <- trajectories
    nr_target_tomato_neg[names(sce_chimera_score_neg_max_table)] <- sce_chimera_score_neg_max_table
    
    nr_WT_tomato <- rep(0, length(trajectories))
    names(nr_WT_tomato) <- trajectories
    nr_WT_tomato[names(chimera_WT_score_max_table)] <- chimera_WT_score_max_table
    
    nr_WT_tomato_neg <- rep(0, length(trajectories))
    names(nr_WT_tomato_neg) <- trajectories
    nr_WT_tomato_neg[names(chimera_WT_score_neg_max_table)] <- chimera_WT_score_neg_max_table
    
    xy <- nr_WT_tomato >= 10
    trajectories_keep <- names(nr_target_tomato)[xy]
    
    nr_target_tomato <- nr_target_tomato[trajectories_keep]
    nr_target_tomato_neg <- nr_target_tomato_neg[trajectories_keep]
    nr_WT_tomato <-  nr_WT_tomato[trajectories_keep]
    nr_WT_tomato_neg <- nr_WT_tomato_neg[trajectories_keep]
    
    ratio_target <- nr_target_tomato/nr_target_tomato_neg
    
    ratio_WT <- nr_WT_tomato/nr_WT_tomato_neg
    
    norm_factor_target <- median(ratio_target)
    norm_factor_WT <- median(ratio_WT)
    
    sample_number_WT <- floor(min(length(chimera_WT_score_max ),length(chimera_WT_score_neg_max)/norm_factor_WT))
    sample_number_target <- floor(min(length(sce_chimera_score_max ),length(sce_chimera_score_neg_max)/norm_factor_target))
    odds_ratio <- list()
    p_values <- list()
    for (k in 1:30){
      chimera_WT_score_max <- chimera_WT_score_max[sample(1:length(chimera_WT_score_max),sample_number_WT)]
      chimera_WT_score_neg_max <- chimera_WT_score_neg_max[sample(1:length(chimera_WT_score_neg_max),sample_number_WT*norm_factor_WT)]
      
      sce_chimera_score_max <- sce_chimera_score_max[sample(1:length(sce_chimera_score_max),sample_number_target)]
      sce_chimera_score_neg_max <- sce_chimera_score_neg_max[sample(1:length(sce_chimera_score_neg_max),sample_number_target*norm_factor_target)]
      
      odds_ratio[[k]] <- rep(NA,length(trajectories_keep))
      p_values[[k]] <- rep(NA,length(trajectories_keep))
      
      for (j in 1:length(trajectories_keep)){
        table_fisher <- c(sum(sce_chimera_score_max== trajectories_keep[j]),sum(sce_chimera_score_neg_max == trajectories_keep[j]),sum(chimera_WT_score_max == trajectories_keep[j]),sum(chimera_WT_score_neg_max== trajectories_keep[j]))
        dim(table_fisher) <- c(2,2)
        aa <- fisher.test(table_fisher)
        p_values[[k]][j] <- aa$p.value
        odds_ratio[[k]][j] <- aa$estimate
      } 
    }
    odds_ratio <- do.call(cbind,odds_ratio)
    odds_ratio <- apply(odds_ratio,1,median)
    p_values <- do.call(cbind,p_values)
    p_values <- apply(p_values,1,median)
    p_values_all[[jk]] <- p_values
    odds_ratio_all[[jk]] <- odds_ratio
  }
  
  odds_ratio <- do.call(cbind,odds_ratio_all)
  odds_ratio <- apply(odds_ratio,1,median)
  p_values <- do.call(cbind,p_values_all)
  p_values <- apply(p_values,1,median)
  p_values <- p.adjust(p_values)
  
  fisher_test_trajectories <- data.frame(trajectory=trajectories_keep,p_values=p_values,odds_ratio=odds_ratio)
  fisher_test_trajectories <- fisher_test_trajectories[order(fisher_test_trajectories$p_values),]
  fisher_test_trajectories$sig <- "enriched"
  fisher_test_trajectories$sig[fisher_test_trajectories$odds_ratio < 1] <- "depleted"
  fisher_test_trajectories$sig[fisher_test_trajectories$p_values> alpha] <- "not significant"
  p <- ggplot( fisher_test_trajectories , aes(x=log10(odds_ratio), y=-log10(p_values+1e-100),color=trajectory)) +
    labs(x=expression(paste("log"[10],"(odds ratio)")),y=expression(paste("-log"[10],"(p_values)"))) +
    scale_color_manual(values=c("not significant"="black","enriched"="darkblue","depleted" = "darkred"),name="") + geom_point()+
    ggrepel::geom_text_repel(data=fisher_test_trajectories[apply(cbind(fisher_test_trajectories$odds_ratio,1/fisher_test_trajectories$odds_ratio),1,max)>1.5,], 
                             aes(x=log10(odds_ratio), y=-log10(p_values+1e-100), label=trajectory,color=sig), max.overlaps=Inf, size=5,segment.color = NA) +
    theme_classic() + scale_y_continuous(trans=pseudo_log_trans(sigma=1,base=10))+
    theme(
      axis.text = element_text(size=rel(0.75), color='black'),
      axis.title = element_text(size=rel(1.0), color='black')) + ggtitle(paste0(target,"- trajectory enrichment"))+annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "blue",alpha=0.1) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "red",alpha=0.1) 
  print(p)
  ggsave(p,file=paste0("da_trajectories_",target,".pdf"))
  return(fisher_test_trajectories)
}