# This script combines the per-sample mapping from map_extended_T_per_sample_rev.R
chimeraT <- readRDS("chimeraT_original_mapping.rds")
chimeraT <- chimeraT[,complete.cases(reducedDims(chimeraT)$pca.corrected.E8.5)]
mappings = lapply(c(1,2,5,6,7,8,9,10),function(x) readRDS(paste0("mapping_T_rev_sample_",toString(x),".rds")))
mappings <- do.call(rbind,mappings)

chimeraT$stage.mapped = as.character(mappings$stage.mapped[match(chimeraT$cell, mappings$cell)])
chimeraT$celltype.mapped = as.character(mappings$celltype.mapped[match(chimeraT$cell, mappings$cell)])
chimeraT$closest.cell = mappings$closest.cell[match(chimeraT$cell, mappings$cell)]
saveRDS(chimeraT,file="../results/chimeraT_2020_extended_mapping_rev.rds")
