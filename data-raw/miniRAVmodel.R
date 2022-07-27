##### Script to build miniRAVmodel #############################################

## Select RAVs to keep in the miniRAVmodel
dir <- "~/data2/GenomicSuperSignaturePaper/inst/extdata"
RAVmodel <- readRDS(file.path(dir, "refinebioRseq_RAVmodel_PLIERpriors_20220117.rds"))

# Validated index
library(bcellViper)
data(bcellViper)
val_all <- validate(dset, RAVmodel)
validated_ind <- validatedSignatures(val_all, RAVmodel, num.out = 3,
                                     swCutoff = 0, indexOnly = TRUE,
                                     filterMessage = FALSE)  # 2538 1139 884

# Keyword-containing index
keyword_ind <- findSignature(RAVmodel, "Bcell", k = 5)  # 695 1994

# RAVs to keep (haven't found the script I got others)
keep_ind <- c(1076, 338, 1467, 1614, 294, 3071, 1694, 438, 725, 1497, 501, 941,
              validated_ind, keyword_ind,
              312, 468) # for TCGA validation


## Subset RAVmodel with 17 RAVs
miniRAVmodel <- RAVmodel[, keep_ind]
cl_ind <- which(metadata(miniRAVmodel)$cluster %in% keep_ind)
metadata(miniRAVmodel)$cluster <- metadata(miniRAVmodel)$cluster[cl_ind]
metadata(miniRAVmodel)$size <- metadata(miniRAVmodel)$size[keep_ind]

## Restructure RAVmodel <========================================= not done yet
# colData(miniRAVmodel)$size <- metadata(miniRAVmodel)$size   # move 'size'
# cl_df <- list()
# cl <- metadata(miniRAVmodel)$cluster
# for (i in colnames(miniRAVmodel)) {   # convert vector into table
#     rav_num <- gsub("RAV", "", i)
#     ind <- which(cl == rav_num)
#     cl_df[[i]] <- names(cl[ind])
# }
# colData(miniRAVmodel)$cluster <- cl_df   # move cluster
# metadata(miniRAVmodel) <- metadata(miniRAVmodel)[3:7]   # remove old slots

## Save miniRAVmodel
# Run tools::checkRdaFiles() to determine the best compression for each file
# save(miniRAVmodel, file = "~/data2/GenomicSuperSignature/data/miniRAVmodel.RData")
# tools::checkRdaFiles("~/data2/GenomicSuperSignature/data/miniRAVmodel.RData")

save(miniRAVmodel, file = "~/data2/GenomicSuperSignature/data/miniRAVmodel.RData")
tools::resaveRdaFiles("~/data2/GenomicSuperSignature/data/miniRAVmodel.RData",
                      version = 3)
