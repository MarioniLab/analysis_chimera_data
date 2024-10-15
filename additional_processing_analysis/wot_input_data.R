# This script extracts cells for each final cell type at E9.25 and stores them in gmt files for 
#input to WOT.

atlas_meta  <- readRDS("../fromIvan_new_atlas/integrated_meta_celltype_clus.rds")
atlas_meta <- atlas_meta[!(atlas_meta$stage%in%c("E9.5","mixed_gastrulation")),]
celltypes_9.25 <- unique(atlas_meta$celltype.clustering[atlas_meta$stage=="E9.25"])
library("ActivePathways")
gene_list <- list()
for (j in 1:length(celltypes_9.25)){
  gene_list[[j]] <- atlas_meta$cell[atlas_meta$celltype.clustering==celltypes_9.25[j]]
}
names(gene_list) <- celltypes_9.25
lapply(gene_list, function(x) write.table( data.frame(x), "gene_list_wot.gmt"  , append= T ))
