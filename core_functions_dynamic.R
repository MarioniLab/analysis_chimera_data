library(destiny)
library(SingleCellExperiment)
library(ggthemes)
library(glmnet)
library(scran)
library(scater)
library(pROC)
library(ggrepel)

find_dynamic_genes <- function(sce,excluded_genes=NULL)
{
  sce <- sce[,!(is.na(sce$stage))]
  stage_vec <- sce$stage
  lf_init_end <- rowMeans(logcounts(sce)[,stage_vec == max(stage_vec)])- rowMeans(logcounts(sce)[,stage_vec == min(stage_vec)])
  sce <- sce[abs(lf_init_end) > 0.5,]
  stage_vec <- sce$stage
  p_vals <- p.adjust(apply(logcounts(sce),1,function(v) summary(lm(v ~ stage_vec))$coefficients[2,4]))
  residuals <- apply(logcounts(sce),1,function(v) lm(v ~ stage_vec)$residuals)
  sce <- sce[p_vals<10^-10,]
  residuals <- residuals[,p_vals<10^-10]
  p_vals_res <- p.adjust(apply(residuals,2,function(v) summary(lm(v ~ sce$batch))$coefficients[2,4]))
  sce <- sce[p_vals_res > 0.2,]
  return(sce)
  }

compute_and_plot_pseudotime <- function(sce){
  dm <- DiffusionMap(t(as.matrix(logcounts(sce))),n_pcs=10)
  cors <- cor(dm@eigenvectors,sce$stage)
  dm@eigenvectors[,cors <0 ] <- -dm@eigenvectors[,cors <0 ]
  sce$dpt <- dm@eigenvectors[,which.max(abs(cors))]
    xx <- sample(1:length(dm$DC1))
  tmp <- data.frame(DC1 = eigenvectors(dm)[xx, 1],
                    DC2 = eigenvectors(dm)[xx, 2],
                    celltype = sce$celltype[xx],
                    dpt = sce$dpt[xx])
  p1 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = celltype)) +
    geom_point() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic(base_size=11) + theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
    labs(color="")+ guides(color = guide_legend(nrow=2,override.aes = list(size = 3)))+ scale_color_colorblind()
  print(p1)
  
  p2 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = dpt)) +
    geom_point(alpha=0.5,size=2) + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic(base_size=11) +scale_color_viridis_c()+ theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin())+
    labs(color="")+ guides(color = guide_legend(override.aes = list(size = 3)))
  print(p2)
  return(sce)
}


correlation_pseudotime <- function(reference_sce,perturbed_sce){
  correlation_matrix <- matrix(ncol=ncol(reference_sce),nrow=ncol(perturbed_sce))
  rownames(correlation_matrix) <- colnames(perturbed_sce)
  colnames(correlation_matrix) <- colnames(reference_sce)
  logcounts_reference <- logcounts(reference_sce)
  logcounts_perturbed <- logcounts(perturbed_sce)
  for (j in 1:ncol(correlation_matrix)){
    for (k in 1:nrow(correlation_matrix)){
      correlation_matrix[k,j] <- cor(logcounts_reference[,j],logcounts_perturbed[,k])
    }
  }
  perturbed_sce$dpt <- apply(correlation_matrix, 1, function(x) 
    mean(reference_sce$dpt[colnames(correlation_matrix)[order(x,decreasing=T) %in% 1:10]]))
    output <- c()
    output$perturbed_sce <- perturbed_sce
    output$correlation_matrix <- correlation_matrix
  return(output)
}

Wilcoxon_test_perturbed_vs_normal <- function(perturbed_sce){
  wilcox_test <- wilcox.test(perturbed_sce$dpt[perturbed_sce$perturbation_status=="perturbed"],
                             perturbed_sce$dpt[perturbed_sce$perturbation_status=="normal"],conf.int = TRUE)
  return(data.frame(lower_bound = wilcox_test$conf.int[1],
                    upper_bound = wilcox_test$conf.int[2],
                    p_value = wilcox_test$p.value,
                    estimate = wilcox_test$estimate ))
}



DE_expression_static <- function(sce_case,sce_control){
  y_case <- sce_case$perturbed
  X_case_static <- logcounts(sce_case)
  y_control <- sce_control$perturbed
  X_control_static <- logcounts(sce_control)

  
  X_static <- cbind(X_control_static,X_case_static)
  perturbation_all <- c(sce_control$perturbed,sce_case$perturbed)
  target <- c(rep(FALSE,length(sce_control$perturbed)),rep(TRUE,length(sce_case$perturbed)))
  genes <- rownames(sce_case)
  estimate_interaction_static <- matrix(ncol=4,nrow=length(genes))
  colnames( estimate_interaction_static) <- c("estimate","std_error","z_value","p_value")
  estimate_static <- matrix(ncol=4,nrow=length(genes))
  colnames( estimate_static) <- c("estimate","std_error","z_value","p_value")
  for (jk in 1:nrow(X_static)){
    temp <- summary(glm(perturbation_all
                        ~X_static[jk,]+target+target*X_static[jk,],family="binomial"))$coefficient
    if (nrow(temp) ==4){
      estimate_interaction_static[jk,] <- temp[4,]}
    temp2 <- summary(glm(sce_case$perturbed
                         ~X_case_static[jk,],family="binomial"))$coefficient
    if (nrow(temp2) ==2){
      estimate_static[jk,] <- temp2[2,]}
  }
  estimate_interaction_static <- as.data.frame(estimate_interaction_static)
  estimate_interaction_static$gene <- rownames(X_static)
  estimate_interaction_static$type <- "static"
  estimate_interaction_static$FDR <- p.adjust(estimate_interaction_static$p_value)
 
  estimate_static <- as.data.frame(estimate_static)
  estimate_static$gene <- rownames(X_static)
  estimate_static$type <- "static"
  estimate_static$FDR <- p.adjust(estimate_static$p_value)
  
 
  
  return(data.frame(static_effect_contrasted_with_control= estimate_interaction_static,
                    static_effect=estimate_static))
}


DE_expression_dynamic <- function(sce_case,sce_control,pseudotime_case,pseudotime_control){
  y_case <- sce_case$perturbed
  X_case_dynamic <- apply(logcounts(sce_case),1,function(x) x*pseudotime_case[colnames(sce_case)])
  colnames(X_case_dynamic) <- paste0(colnames(X_case_dynamic),"- dynamic")
  
  y_control <- sce_control$perturbed
  X_control_dynamic <- apply(logcounts(sce_control),1,function(x) x*pseudotime_control[colnames(sce_control)])
  colnames(X_control_dynamic) <- paste0(colnames(X_control_dynamic),"- dynamic")
  
  perturbation_all <- c(sce_control$perturbed,sce_case$perturbed)
  target <- c(rep(FALSE,length(sce_control$perturbed)),rep(TRUE,length(sce_case$perturbed)))
  genes <- rownames(sce_case)
  
  X_dynamic <- rbind(X_control_dynamic,X_case_dynamic)
  estimate_interaction_dynamic <- matrix(ncol=4,nrow=length(genes))
  colnames( estimate_interaction_dynamic) <- c("estimate","std_error","z_value","p_value")
  estimate_dynamic <- matrix(ncol=4,nrow=length(genes))
  colnames( estimate_dynamic) <- c("estimate","std_error","z_value","p_value")
  for (jk in 1:ncol(X_dynamic)){
    temp <- summary(glm(perturbation_all
                        ~X_dynamic[,jk]+target+target*X_dynamic[,jk],family="binomial"))$coefficient
    if (nrow(temp) ==4){
      estimate_interaction_dynamic[jk,] <- temp[4,]
    }
    temp2 <- summary(glm(sce_case$perturbed
                         ~X_case_dynamic[,jk],family="binomial"))$coefficient
    if (nrow(temp2) ==2){
      estimate_dynamic[jk,] <- temp2[2,]
    }
  }
  
  estimate_interaction_dynamic <- as.data.frame(estimate_interaction_dynamic)
  estimate_interaction_dynamic$gene <- colnames(X_dynamic)
  estimate_interaction_dynamic$type <- "dynamic"
  estimate_interaction_dynamic$FDR <- p.adjust(estimate_interaction_dynamic$p_value)
  
  estimate_dynamic <- as.data.frame(estimate_dynamic)
  estimate_dynamic$gene <- colnames(X_dynamic)
  estimate_dynamic$type <- "dynamic"
  estimate_dynamic$FDR <- p.adjust(estimate_dynamic$p_value)
  
  return(data.frame(dynamic_effect_contrasted_with_control= estimate_interaction_dynamic,
                    dynamic_effect=estimate_dynamic))
}

volcano_plot_static <- function(df_effect, FDR_static=1,max_highlight=20,highlight_extra_genes=NULL){
  df_static_interaction <- df_effect[,c("static_effect_contrasted_with_control.estimate",
                                        "static_effect_contrasted_with_control.FDR","static_effect_contrasted_with_control.gene" )]
  colnames(df_static_interaction) <- c("logFC","FDR","gene")
  
  
  df_static_simple_interaction <- cbind(df_static_interaction,
                                        df_effect$static_effect.estimate
                                      [match(df_effect$static_effect_contrasted_with_control.gene,
                                             df_effect$static_effect.gene)],
                                      df_effect$static_effect.FDR
                                      [match(df_effect$static_effect_contrasted_with_control.gene,
                                             df_effect$static_effect.gene)])
  colnames(df_static_simple_interaction) <- c("logFC_interaction","FDR_interaction","gene","logFC_simple","FDR_simple")
  df_static_simple_interaction$sig <- "not significant"
  df_static_simple_interaction$sig[df_static_simple_interaction$FDR_interaction < FDR_static] <- "interaction"
  df_static_simple_interaction$sig[df_static_simple_interaction$FDR_simple < FDR_static] <- "simple"
  df_static_simple_interaction$sig[df_static_simple_interaction$FDR_simple < FDR_static &
                                     df_static_simple_interaction$FDR_interaction < FDR_static] <- "both"
  max_highlight <- min(max_highlight,nrow(df_static_simple_interaction))
  highlight_static <- df_static_simple_interaction$FDR_interaction <FDR_static & df_static_simple_interaction$FDR_simple <FDR_static & 
    (df_static_simple_interaction$logFC_interaction <= sort(df_static_simple_interaction$logFC_interaction)[max_highlight]|
       df_static_simple_interaction$logFC_interaction >= sort(df_static_simple_interaction$logFC_interaction,decreasing=TRUE)[max_highlight])
  highlight_static[highlight_static & (df_static_simple_interaction$logFC_interaction > 0)] <- "up"
  highlight_static[highlight_static==TRUE & df_static_simple_interaction$logFC_interaction<0] <- "down"
  highlight_static[df_static_simple_interaction$gene %in% highlight_extra_genes] <- "gene_of_interest"
  df_static_simple_interaction$highlight <- highlight_static
  p_static_interaction <- ggplot(df_static_simple_interaction, aes(x=logFC_interaction, y=-log10(FDR_interaction),color=highlight)) +
    geom_point(size=3, alpha=0.4) +
    ggtitle('static') + scale_color_manual(values=c("FALSE"="grey","up"="darkblue","down"="darkred","gene_of_interest"="purple"))+
    labs(y=expression('-Log'[10]*' FDR'), x=expression('logFC_interaction')) +
    theme_classic(base_size=20) +
    theme(legend.position="none", plot.title = element_text(size = rel(1), hjust = 0.5))+
    geom_text_repel(data=df_static_simple_interaction[!(highlight_static==FALSE),],
                    aes(x = logFC_interaction, y = -log10(FDR_interaction),label=gene),max.overlaps=100)+ 
    scale_y_continuous(trans=scales::pseudo_log_trans(sigma=5,base = 1.5))+
    geom_hline(yintercept = 1)
    
  return(p_static_interaction)
}

volcano_plot_dynamic <- function(df_effect, FDR_dynamic=1,max_highlight=20,highlight_extra_genes=NULL){
  df_dynamic_interaction <- df_effect[,c("dynamic_effect_contrasted_with_control.z_value",
                                        "dynamic_effect_contrasted_with_control.FDR","dynamic_effect_contrasted_with_control.gene" )]
  colnames(df_dynamic_interaction) <- c("z_score","FDR","gene")
  
  df_dynamic_simple_interaction <- cbind(df_dynamic_interaction,df_effect$dynamic_effect.z_value
                                        [match(df_effect$dynamic_effect_contrasted_with_control.gene,
                                               df_effect$dynamic_effect.gene)],
                                        df_effect$dynamic_effect.FDR
                                        [match(df_effect$dynamic_effect_contrasted_with_control.gene,
                                               df_effect$dynamic_effect.gene)])
  colnames(df_dynamic_simple_interaction) <- c("z_score_interaction","FDR_interaction","gene","z_score_simple","FDR_simple")
  df_dynamic_simple_interaction$gene <- sapply(df_dynamic_simple_interaction$gene,function(x) strsplit(x,"- dynamic")[[1]][1])
  df_dynamic_simple_interaction$sig <- "not significant"
  df_dynamic_simple_interaction$sig[df_dynamic_simple_interaction$FDR_interaction < FDR_dynamic] <- "interaction"
  df_dynamic_simple_interaction$sig[df_dynamic_simple_interaction$FDR_simple < FDR_dynamic] <- "simple"
  df_dynamic_simple_interaction$sig[df_dynamic_simple_interaction$FDR_simple < FDR_dynamic &
                                     df_dynamic_simple_interaction$FDR_interaction < FDR_dynamic] <- "both"
  
  highlight_dynamic <- df_dynamic_simple_interaction$FDR_interaction <FDR_dynamic & df_dynamic_simple_interaction$FDR_simple <FDR_dynamic & 
    (df_dynamic_simple_interaction$z_score_interaction <= sort(df_dynamic_simple_interaction$z_score_interaction)[max_highlight]|
       df_dynamic_simple_interaction$z_score_interaction >= sort(df_dynamic_simple_interaction$z_score_interaction,decreasing=TRUE)[max_highlight])
  highlight_dynamic[highlight_dynamic & (df_dynamic_simple_interaction$z_score_interaction > 0)] <- "up"
  highlight_dynamic[highlight_dynamic==TRUE & df_dynamic_simple_interaction$z_score_interaction<0] <- "down"
  highlight_dynamic[df_dynamic_simple_interaction$gene %in% highlight_extra_genes] <- "gene_of_interest"
  df_dynamic_simple_interaction$highlight <- highlight_dynamic
  p_dynamic_interaction <- ggplot(df_dynamic_simple_interaction, aes(x=z_score_interaction, y=-log10(FDR_interaction),color=highlight)) +
    geom_point(size=3, alpha=0.4) +
    ggtitle('dynamic') + scale_color_manual(values=c("FALSE"="grey","up"="darkblue","down"="darkred","gene_of_interest"="purple"))+
    labs(y=expression('-Log'[10]*' FDR'), x=expression('z_score_interaction')) +
    theme_classic(base_size=20) +
    theme(legend.position="none", plot.title = element_text(size = rel(1), hjust = 0.5))+
    geom_text_repel(data=df_dynamic_simple_interaction[!(highlight_dynamic==FALSE),],
                    aes(x = z_score_interaction, y = -log10(FDR_interaction),label=gene),max.overlaps=100)+ 
    scale_y_continuous(trans=scales::pseudo_log_trans(sigma=5,base = 1.5))+
    geom_hline(yintercept = 1)
 return(p_dynamic_interaction)
  
}

refine_celltype_stage <- function(sce_chimera,cor_matrix){
  k.mapped <- matrix(0,nrow=ncol(cor_matrix),10)
  for (j in 1:ncol(cor_matrix))
  {
    k.mapped[j,] <- rownames(cor_matrix)[order(cor_matrix[,j],decreasing=TRUE)[1:10]]
  }
  celltypes <-  matrix(0,nrow=ncol(cor_matrix),ncol=10)
  stages <-   matrix(0,nrow=ncol(cor_matrix),ncol=10)
  for (j in 1:ncol(cor_matrix))
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
  sce_chimera$stage.refined <- NA
  sce_chimera$stage.refined[match(colnames(cor_matrix),colnames(sce_chimera))] <- stage.mapped
  sce_chimera$celltype.refined <- NA
  sce_chimera$celltype.refined[match(colnames(cor_matrix),colnames(sce_chimera))] <- celltype.mapped
  return(sce_chimera)
}

