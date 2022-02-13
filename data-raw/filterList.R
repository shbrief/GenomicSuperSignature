## The list of 'low' quality RAVs
## Output containing RAVs belong to this list will have a message.
## 1) single-element clusters (`s_ind`)
## 2) RAVs with no or too-many enriched gene sets (`gsea_filtered_ind`).
## 3) PCs from only one study (`single_study_ind`)


##### MSigDB C2 ################################################################
dir <- "~/data2/GenomicSuperSignatureLibrary/refinebioRseq/RAVmodel_536"
RAVmodel <- readRDS(file.path(dir, "refinebioRseq_RAVmodel_C2_20220115.rds"))

## Single-element clusters
s_ind <- which(metadata(RAVmodel)$size == 1)

## No or too many (5% of input) enriched gene sets
num_enriched_gs <- sapply(gsea(RAVmodel), nrow)
no_gs_ind <- which(num_enriched_gs == 0)
cut <- 276   ##### 5% of 5,529 MSigDB C2 gene sets
too_many_gs_ind <- which(num_enriched_gs > cut)
gsea_filtered_ind <- unique(c(no_gs_ind, too_many_gs_ind))

## Single-study clusters
ns_ind <- which(metadata(RAVmodel)$size != 1)
s_study <- which(sapply(studies(RAVmodel), length) == 1)
single_study_ind <- intersect(ns_ind, s_study)

##### PLIERpriors ##############################################################
dir <- "~/data2/GenomicSuperSignatureLibrary/refinebioRseq/RAVmodel_536"
RAVmodel <- readRDS(file.path(dir, "refinebioRseq_RAVmodel_PLIERpriors_20220117.rds"))

## No or too many (5% of input) enriched gene sets
num_enriched_gs <- sapply(gsea(RAVmodel), nrow)
no_gs_ind <- which(num_enriched_gs == 0)
cut <- 31   ##### 5% of all 628 PLIERpriors gene sets
too_many_gs_ind <- which(num_enriched_gs > cut)
gsea_filtered_ind_plier <- unique(c(no_gs_ind, too_many_gs_ind))

##### Save #####################################################################
filterList <- list("Cluster_Size_filter" = s_ind,
                   "GSEA_C2_filter" = gsea_filtered_ind,
                   "GSEA_PLIERpriors_filter" = gsea_filtered_ind_plier,
                   "Redundancy_filter" = single_study_ind)
save(filterList, file = "data/filterList.RData")

