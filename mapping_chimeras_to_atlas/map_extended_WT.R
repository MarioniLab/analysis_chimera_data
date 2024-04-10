# This script combines the per-sample mapping from map_extended_WT_per_sample_rev.R
chimeraWT <- readRDS("chimeraWT_original_mapping.rds")
names(colData(chimeraWT))[names(colData(chimeraWT)) == "Sample"] <- "sample"
chimeraWT <- chimeraWT[,complete.cases(reducedDims(chimeraWT)$pca.corrected.E8.5)]
mappings = lapply(c(5,6,7,8,9,10),function(x) readRDS(paste0("mapping_WT_rev_sample_",toString(x),".rds")))
mappings <- do.call(rbind,mappings)

chimeraWT$stage.mapped = as.character(mappings$stage.mapped[match(chimeraWT$cell, mappings$cell)])
chimeraWT$celltype.mapped = as.character(mappings$celltype.mapped[match(chimeraWT$cell, mappings$cell)])
chimeraWT$closest.cell = mappings$closest.cell[match(chimeraWT$cell, mappings$cell)]
saveRDS(chimeraWT,file="../results/chimeraWT_2020_extended_mapping_rev.rds")
