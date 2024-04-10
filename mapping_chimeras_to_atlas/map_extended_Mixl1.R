# This script combines the per-sample mapping from map_extended_Mixl1_per_sample_rev.R
mappings = lapply(c(1,2,3,4,5,6),function(x) readRDS(paste0("mapping_Mixl1_rev_sample_",toString(x),".rds")))
mappings <- do.call(rbind,mappings)

chimeraMixl1$stage.mapped = as.character(mappings$stage.mapped[match(chimeraMixl1$cell, mappings$cell)])
chimeraMixl1$celltype.mapped = as.character(mappings$celltype.mapped[match(chimeraMixl1$cell, mappings$cell)])
chimeraMixl1$closest.cell = mappings$closest.cell[match(chimeraMixl1$cell, mappings$cell)]
saveRDS(chimeraMixl1,file="../results/chimeraMixl1_2020_extended_mapping_rev.rds")
