#A set of core functions for the analysis of chimera data sets with the 
# extended mouse gastrulation atlas as the reference set


library(scran)
library(scater)
library(SingleCellExperiment)
library(viridis)


#' Identifies highly variable genes in chimera data sets. 
#' @title getHVGs
#' @param sce chimera data set as SingleCellExperiment
#' @param min_mean minimum mean expression level above which a gene is included
#' in the variance modelling to find the highly variable genes. The default is 0.001.
#' @param FDR FDR threshold for the variance modelling, cutoff for including genes 
#' in the highly variable gene list
#' @return vector of highly variable genes
#' @export
getHVGs <- function(sce,min_mean = 0.001,FDR = 0.01,yGenes=read.table("data/ygenes.tab", stringsAsFactors = FALSE)[,1])
{
  
  exclude <- c("ENSMUSG00000086503",yGenes,"tomato-td")#
  stats <- modelGeneVar(sce,subset.row = setdiff(rownames(sce),exclude),
                        block=sce$sample)
  stats <- stats[stats$mean > min_mean,]
  top_genes <- rownames(stats[stats$FDR < FDR,])
  return(top_genes)
} 


#' Identifies the most frequent label (e.g. cell type) among the
#' k closest cells of a chimera cell in the atlas
#' taken from
# https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
#' @title getmode
#' @param v vector of labels, such as cell type
#' @param dist distance of the cells from the chimera cell
#' @return most frequent label (in case of ties the one of the neighbour cell
#' closest to the chimera cell)
#' @export
getmode <- function(v, dist) {
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


###

###

plot_target_gene_expression <- function(sce, target_gene)
{
  target_count <- counts(sce)[rowData(sce)$SYMBOL == target_gene,]
  sce$target_count <- target_count > 0
  umap_pl <-plotReducedDim(sce,dimred="umap",colour_by="target_count",text_by="celltype",
                           point_size=0.5,text_size=3,point_alpha=0.3,text_colour="black") +
    scale_color_manual(values=c("lightgrey","blue"))+
    labs(color=paste0(target_gene," expressed"))
  
  print(umap_pl)
}



plot_celltypes <- function(chimeraTarget,saveFileName="temp.pdf",target=""){
  Target_tomato_neg <- chimeraTarget[,chimeraTarget$tomato==F]
  Target_tomato_pos <- chimeraTarget[,chimeraTarget$tomato==T]
  Target_tomato_pos <- Target_tomato_pos[,Target_tomato_pos$celltype.mapped %in% unique(Target_tomato_neg$celltype.mapped)]
  prop_celltypes_Target_8.5 <- 100*c(as.vector(table(factor(Target_tomato_neg$celltype.mapped,levels=unique(Target_tomato_neg$celltype.mapped))))/
                                       length(Target_tomato_neg$celltype.mapped),
  as.vector(table(factor(Target_tomato_pos$celltype.mapped,levels=unique(Target_tomato_neg$celltype.mapped))))/
    length(Target_tomato_pos$celltype.mapped))
  xx <- order(prop_celltypes_Target_8.5[1:(length(prop_celltypes_Target_8.5)*0.5)])
  prop_celltypes_Target_8.5 <- c(prop_celltypes_Target_8.5[xx],prop_celltypes_Target_8.5[xx+(length(prop_celltypes_Target_8.5)*0.5)])
  celltypes_Target_8.5 <- factor(rep(names(table(factor(Target_tomato_neg$celltype.mapped,
                          levels=unique(Target_tomato_neg$celltype.mapped))))[xx],2),
                      levels=names(table(factor(Target_tomato_neg$celltype.mapped,
                                                levels=unique(Target_tomato_neg$celltype.mapped))))[xx])
  dataset <- c(rep("tomato-neg",length(unique(Target_tomato_neg$celltype.mapped))),
               rep("tomato-pos",length(unique(Target_tomato_neg$celltype.mapped))))
  p1 = ggplot(mapping = aes(x = celltypes_Target_8.5,y=prop_celltypes_Target_8.5, fill = dataset)) +
    geom_bar(stat="identity",position="dodge") +
    scale_fill_manual(values=c( "tomato-neg" = "darkblue","tomato-pos" = "red"))+
    labs(x = "", y = "% of cells") +
    theme(text = element_text(size=15)) +
    scale_x_discrete(breaks = as.factor(celltypes_Target_8.5[1:(0.5*length(celltypes_Target_8.5))]), 
                     drop = FALSE) + coord_flip() + ggtitle(target)
  print(p1)
  ggsave(p1,file=saveFileName,width=10,height=20)
}

plot_celltypes_DA <- function(chimeraTarget,saveFileName="temp.pdf",target=""){
  # plot cellytpes for which chimera knockout cells are over- or underrepresented by 
  # a factor of at least 1.2
  # not using normalisation 
  Target_tomato_neg <- chimeraTarget[,chimeraTarget$tomato==F]
  Target_tomato_pos <- chimeraTarget[,chimeraTarget$tomato==T]
  celltypes <- unique(chimeraTarget$celltype.mapped)
  celltypes <- celltypes[!(grepl("ExE",celltypes))]
  prop_celltypes_Target_8.5 <- 100*c(as.vector(table(factor(Target_tomato_neg$celltype.mapped,levels=celltypes)))/length(Target_tomato_neg$celltype.mapped),
                                     as.vector(table(factor(Target_tomato_pos$celltype.mapped,levels=celltypes)))/length(Target_tomato_pos$celltype.mapped))
  prop_celltypes_Target_8.5_rel <- prop_celltypes_Target_8.5[(1+0.5*length(prop_celltypes_Target_8.5)):length(prop_celltypes_Target_8.5)]/prop_celltypes_Target_8.5[1:(0.5*length(prop_celltypes_Target_8.5))]
  xx_1 <- which((prop_celltypes_Target_8.5_rel > 1.2) | (prop_celltypes_Target_8.5_rel < 1/1.2))
  prop_celltypes_Target_8.5_sum <- prop_celltypes_Target_8.5[(1+0.5*length(prop_celltypes_Target_8.5)):length(prop_celltypes_Target_8.5)]+prop_celltypes_Target_8.5[1:(0.5*length(prop_celltypes_Target_8.5))]
  xx_1 <- intersect(xx_1,which(prop_celltypes_Target_8.5_sum > 1))
  xx <- c(xx_1,xx_1+length(prop_celltypes_Target_8.5_rel))
  celltypes_Target_8.5 <- rep(celltypes,2)
  dataset <- c(rep("tomato-neg",length(unique(Target_tomato_neg$celltype.mapped))),
               rep("tomato-pos",length(unique(Target_tomato_pos$celltype.mapped))))
  p1 = ggplot(mapping = aes(x = celltypes_Target_8.5[xx],y=prop_celltypes_Target_8.5[xx], fill = dataset[xx])) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c( "tomato-neg" = "darkblue","tomato-pos" = "red"), name = paste0(target," chimera\n data sets") )+
    labs(x = "", y = "% of cells") +
    theme(text = element_text(size=15)) +
    scale_x_discrete(breaks = unique(celltypes_Target_8.5)[order(unique(celltypes_Target_8.5))], 
                     drop = FALSE) + coord_flip()
  print(p1)
  ggsave(p1,file=saveFileName,width=10,height=20)
}

plot_stage <- function(chimeraTarget,saveFileName="temp.pdf",target="target",scale="relative"){
  Target_tomato_neg <- chimeraTarget[,chimeraTarget$tomato==F]
  Target_tomato_pos <- chimeraTarget[,chimeraTarget$tomato==T]
  levels <- unique(c(Target_tomato_neg$stage.mapped,Target_tomato_pos$stage.mapped))
  if (scale =="absolute"){
    prop_stages_Target_8.5 <- 100*c(as.vector(table(factor(Target_tomato_neg$stage.mapped, levels=levels)))/
                                      length(chimeraTarget$stage.mapped),
                                    as.vector(table(factor(Target_tomato_pos$stage.mapped, levels=levels)))/
                                      length(chimeraTarget$stage.mapped))
  }
  if (scale =="relative"){
    prop_stages_Target_8.5 <- 100*c(as.vector(table(factor(Target_tomato_neg$stage.mapped, levels=levels)))/
                                      length(Target_tomato_neg$stage.mapped),
                                    as.vector(table(factor(Target_tomato_pos$stage.mapped, levels=levels)))/
                                      length(Target_tomato_pos$stage.mapped))
  }
  stages_Target_8.5 <- c(names(table(factor(Target_tomato_neg$stage.mapped, levels=levels))),
                         names(table(factor(Target_tomato_pos$stage.mapped, levels=levels))))
  dataset <- c(rep("tomato-neg",length(stages_Target_8.5)*0.5),
               rep("tomato-pos",length(stages_Target_8.5)*0.5))
  p1 = ggplot(mapping = aes(x = stages_Target_8.5,y=prop_stages_Target_8.5, fill = dataset)) +
    geom_bar(stat="identity",position="dodge") +
    scale_fill_manual(values=c( "tomato-neg" = "darkblue","tomato-pos" = "red"), name = target)+
    labs(x = "", y = "% of cells") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.1),text = element_text(size=12)) +
    scale_x_discrete(breaks = unique(stages_Target_8.5)[order(unique(stages_Target_8.5))], 
                     drop = FALSE) +ylim(0,35)+ggtitle(target)
  print(p1)
  ggsave(p1,file=saveFileName,width=20)
}





cell_cycle_scoring <- function(sce)
{
  sce$S.Score <- rep(NA,ncol(sce))
  sce$G2M.Score <- rep(NA,ncol(sce))
  sce$phase <- rep(NA,ncol(sce))
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  s.genes = getLDS(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", 
                   values = cc.genes$s.genes, mart = human, attributesL = c("ensembl_gene_id"), 
                   martL = mouse, uniqueRows=T)[,2]
  g2m.genes = getLDS(attributes = c("ensembl_gene_id"), filters = "hgnc_symbol", 
                     values = cc.genes$g2m.genes, mart = human, attributesL = c("ensembl_gene_id"), 
                     martL = mouse, uniqueRows=T)[,2]
  for (j in 1:length(unique(sce$sample))){
    print(j)
    jj <- unique(sce$sample)[j]
    SeuratObj <- Seurat::CreateSeuratObject(as.matrix(counts(sce[,sce$sample ==jj])))
    SeuratObj <- Seurat::NormalizeData(SeuratObj)
    SeuratObj <- Seurat::FindVariableFeatures(SeuratObj, selection.method = "vst")
    SeuratObj <- Seurat::ScaleData(SeuratObj, features = rownames(SeuratObj))
    SeuratObj <- Seurat::CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    sce$S.Score[sce$sample == jj] <- SeuratObj$S.Score
    sce$G2M.Score[sce$sample == jj] <- SeuratObj$G2M.Score
    sce$phase[sce$sample == jj]  <- SeuratObj$Phase
  }
  return(sce)
}



#The following colour-scheme was provided by I. Imaz-Rosshandler (Imaz-Rosshandler et al. (2023). bioRxiv. )
celltype_colours_final = c(
  "Epiblast" = "#635547",
  "Primitive Streak" = "#DABE99",
  "Caudal epiblast" = "#9E6762",
  "PGC" = "#FACB12",
  "Anterior Primitive Streak" = "#C19F70",
  "Node" = "#153B3D",
  "Notochord" = "#0F4A9C",
  "Gut tube" = "#EF5A9D",
  "Hindgut" = "#F397C0",
  "Midgut" = "#FF00B2",
  "Foregut" = "#FFB7FF",
  "Pharyngeal endoderm" = "#95E1FF",
  "Thyroid primordium" = "#97BAD3",
  "Nascent mesoderm" = "#C594BF",
  "Intermediate mesoderm" = "#139992",
  "Caudal mesoderm" = "#3F84AA",
  "Lateral plate mesoderm" = "#F9DFE6",
  "Limb mesoderm" = "#E35F82",
  "Forelimb" = "#D02D75",
  "Kidney primordium" = "#E85639",
  "Presomitic mesoderm" = "#5581CA", #"0000ff", #"blue",
  "Somitic mesoderm" = "#005579",
  "Posterior somitic tissues" = "#5ADBE4", #"40e0d0",#"turquoise",
  "Paraxial mesoderm" = "#8DB5CE",
  "Cranial mesoderm" = "#456722",#"#006400",#darkgreen",
  "Anterior somitic tissues" = "#D5E839",
  "Sclerotome" = "#E3CB3A", #"ffff00", #"yellow",
  "Dermomyotome" = "#00BFC4", #"a52a2a", #"brown",
  "Pharyngeal mesoderm" = "#C9EBFB",
  "Cardiopharyngeal progenitors" = "#556789",
  "Anterior cardiopharyngeal progenitors" = "#683ED8",
  "Allantois" = "#532C8A",
  "Mesenchyme" = "#CC7818",
  "YS mesothelium" = "#FF7F9C",
  "Epicardium" = "#F79083",
  "Embryo proper mesothelium" = "#FF487D",
  "Cardiopharyngeal progenitors FHF" = "#D780B0",
  "Cardiomyocytes FHF 1" = "#A64D7E",
  "Cardiomyocytes FHF 2" = "#B51D8D",
  "Cardiopharyngeal progenitors SHF" = "#4B7193",
  "Cardiomyocytes SHF 1" = "#5D70DC",
  "Cardiomyocytes SHF 2" = "#332C6C",
  "Haematoendothelial progenitors" = "#FBBE92",
  "Blood progenitors" = "#6C4B4C",
  "Erythroid" = "#C72228",
  "Chorioallantoic-derived erythroid progenitors" = "#E50000",
  "Megakaryocyte progenitors" = "#E3CB3A",
  "MEP" = "#EF4E22",
  "EMP" = "#7C2A47",
  "YS endothelium" = "#FF891C",
  "YS mesothelium-derived endothelial progenitors" = "#AE3F3F",
  "Allantois endothelium" = "#2F4A60",
  "Embryo proper endothelium" = "#90E3BF",
  "Venous endothelium" = "#BD3400",
  "Endocardium" = "#9D0049",
  "NMPs/Mesoderm-biased" = "#89C1F5",
  "NMPs" = "#8EC792",
  "Ectoderm" = "#FF675C",
  "Optic vesicle" = "#BD7300",
  "Ventral forebrain progenitors" = "#A0B689",
  "Early dorsal forebrain progenitors" = "#0F8073",
  "Late dorsal forebrain progenitors" = "#7A9941",
  "Midbrain/Hindbrain boundary" = "#8AB3B5",
  "Midbrain progenitors" = "#9BF981",
  "Dorsal midbrain neurons" = "#12ED4C",
  "Ventral hindbrain progenitors" = "#7E907A",
  "Dorsal hindbrain progenitors" = "#2C6521",
  "Hindbrain floor plate" = "#BF9DA8",
  "Hindbrain neural progenitors" = "#59B545",
  "Neural tube" = "#233629",
  "Migratory neural crest" = "#4A6798",
  "Branchial arch neural crest" = "#BD84B0",
  "Frontonasal mesenchyme" = "#D3B1B1",
  "Spinal cord progenitors" = "#6B2035",
  "Dorsal spinal cord progenitors" = "#E273D6",
  "Non-neural ectoderm" = "#F7F79E",
  "Surface ectoderm" = "#FCFF00",
  "Epidermis" = "#FFF335",
  "Limb ectoderm" = "#FFD731",
  "Amniotic ectoderm" = "#DBB400",
  "Placodal ectoderm" = "#FF5C00",
  "Otic placode" = "#F1A262",
  "Otic neural progenitors" = "#00B000",
  "Visceral endoderm" = "#F6BFCB",
  "ExE endoderm" = "#7F6874",
  "ExE ectoderm" = "#989898",
  "Parietal endoderm" = "#1A1A1A"
)

colours_up_down_regulation = c("up in T only" = "#F392E5", "up in Mixl only"="#FEED07", "up in T and Mixl1" = "#FE0707",
                              "down in T only" = "#32CBF6", "down in Mixl1 only" ="#303B30", "down in T and Mixl1"="#1B1FFC",
                              "down in T, up in Mixl1" = "#27AC23", "up in T, down in Mixl1" ="#A24BCE")


                              
