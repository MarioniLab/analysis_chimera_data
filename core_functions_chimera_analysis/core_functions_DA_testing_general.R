# Core functions for perturbSuite_DA

library(SingleCellExperiment)
library(patchwork)
library(scran)
library(scater)
library(scales)

#' @title da_per_celltype
#' @param sce_case SingleCellExperiment for the case data set, 
#' colData needs to include the slots 'perturbed' (TRUE/FALSE) for the presence
#' of a perturbation (CRISPR, drug, disease etc.)
#' @param sce_control SingleCellExperiment for the control data set, 
#' colData needs to include the slots 'perturbed' (TRUE/FALSE) for the presence
#' of a perturbation (CRISPR, drug, disease etc.)
#' @param alpha significance level, default is 10% FDR
#' @return data.frame with the following columns: celltype, odds_ratio-odds ratio of the percentage of 
#' perturbed cells among the cells in the case data versus the respective ratio for the control experiment, 
#' p_values, sig- whether a celltype is significantly enriched, depleted, or not enriched
#' @param thresh_perturbed_control minimum number of perturbed cells per cell type in the control data set
#' @param thresh_unperturbed_case minimum number of unperturbed cells per cell type in the case data set
#' @export
#' @examples 


da_per_celltype <- function(sce_case,sce_control,alpha=0.1,thresh_perturbed_control=5,thresh_unperturbed_case=5){
  sce_case_perturbed <- sce_case[,sce_case$perturbed]
  sce_case_unperturbed <- sce_case[,!(sce_case$perturbed)]
  sce_control_perturbed <- sce_control[,sce_control$perturbed]
  sce_control_unperturbed <- sce_control[,!(sce_control$perturbed)]
  celltypes <- names(table(sce_control_perturbed$celltype))[table(sce_control_perturbed$celltype) >= thresh_perturbed_control]
  celltypes_1 <- names(table(sce_case_unperturbed$celltype))[table(sce_case_unperturbed$celltype) >= thresh_unperturbed_case]
  celltypes <- intersect(celltypes,celltypes_1)
  ratio_control <- rep(0,length(celltypes))
  ratio_target <- rep(0, length(celltypes))
  for (j in 1:length(celltypes)){
    ratio_control[j] <- sum(sce_control_perturbed$celltype == celltypes[j])/sum(sce_control_unperturbed$celltype == celltypes[j])
    ratio_target[j] <- sum(sce_case_perturbed$celltype == celltypes[j])/sum(sce_case_unperturbed$celltype == celltypes[j])
  } 
  names(ratio_target) <- celltypes
  names(ratio_control) <- celltypes
  
  norm_factor_control <- median(ratio_control[(!(is.infinite(ratio_control)))&(!(is.na(ratio_control)))])
  norm_factor_target <- median(ratio_target[(!(is.infinite(ratio_target)))&(!(is.na(ratio_target)))]) 
  sample_number_target <- floor(min(ncol(sce_case_perturbed),ncol(sce_case_unperturbed)/norm_factor_target))
  sample_number_control <- floor(min(ncol(sce_control_perturbed),ncol(sce_control_unperturbed)/norm_factor_control))
  p_values <- list()
  odds_ratio <- list()
  for (k in 1:100){
    sce_case_perturbed <- sce_case_perturbed[,sample(1:ncol(sce_case_perturbed),sample_number_target)]
    sce_case_unperturbed <- sce_case_unperturbed[,sample(1:ncol(sce_case_unperturbed),floor(sample_number_target*norm_factor_target))]
    
    sce_control_perturbed <- sce_control_perturbed[,sample(1:ncol(sce_control_perturbed),sample_number_control)]
    sce_control_unperturbed <- sce_control_unperturbed[,sample(1:ncol(sce_control_unperturbed),floor(sample_number_control*norm_factor_control))]
    
    odds_ratio[[k]] <- rep(0,length(celltypes))
    p_values[[k]] <- rep(0,length(celltypes))
    for (j in 1:length(celltypes)){
      table_fisher <- c(sum(sce_case_perturbed$celltype == celltypes[j]),sum(sce_case_unperturbed$celltype == celltypes[j]),
                        sum(sce_control_perturbed$celltype == celltypes[j]),sum(sce_control_unperturbed$celltype == celltypes[j]))
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
      axis.title = element_text(size=rel(1.0), color='black')) +
    annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "blue",alpha=0.1) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = -Inf, fill= "red",alpha=0.1) 
  print(p)
  return(fisher_test_celltypes)
}



