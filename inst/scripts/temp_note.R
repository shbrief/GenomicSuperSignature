## Update the package name from 'PCAGenomicSignatures' to 'GenomicSuperSignature'
fname <- "~/data2/GenomicSuperSignature/vignettes/PCAmodel_C2.rds"
PCAmodel <- readRDS(fname)
RVmodel <- GenomicSuperSignature::PCAGenomicSignatures(assays = list(model = model(PCAmodel)))
metadata(RVmodel) <- GenomicSuperSignature::metadata(PCAmodel)
metadata(RVmodel)$version <- ">= 0.0.7"
studies(RVmodel) <- GenomicSuperSignature::studies(PCAmodel)
silhouetteWidth(RVmodel) <- GenomicSuperSignature::silhouetteWidth(PCAmodel)
trainingData(RVmodel) <- GenomicSuperSignature::trainingData(PCAmodel)
gsea(RVmodel) <- GenomicSuperSignature::gsea(PCAmodel)
saveRDS(RVmodel, "~/data2/GenomicSuperSignature/vignettes/RAVmodel_C2.rds")


## Update the avgLoading name 'RAV'
fname <- "~/data2/GenomicSuperSignature/vignettes/RAVmodel_C2.rds"
RVmodel <- readRDS(fname)
names(metadata(RVmodel)$size) <- gsub("PCcluster", "RAV", names(metadata(RVmodel)$size))   # before
names(metadata(RVmodel)$size) <- paste0("RAV", 1:length(metadata(RVmodel)$size))
colnames(RVmodel) <- gsub("PCcluster", "RAV", colnames(RVmodel))
names(gsea(RVmodel)) <- gsub("PCcluster", "RAV", names(gsea(RVmodel)))
names(colData(RVmodel))[1] <- "RAV"
saveRDS(RVmodel, "~/data2/GenomicSuperSignature/vignettes/RAVmodel_C2.rds")




## README.md badges
[![Travis-CI Build Status](https://travis-ci.com/shbrief/GenomicSuperSignature.svg?branch=master)](https://travis-ci.org/shbrief/GenomicSuperSignature)
[![Coverage Status](https://codecov.io/github/shbrief/GenomicSuperSignature/coverage.svg?branch=master)](https://codecov.io/github/shbrief/GenomicSuperSignature?branch=master)
[![R build status](https://github.com/shbrief/GenomicSuperSignature/workflows/R-CMD-check/badge.svg)](https://github.com/shbrief/GenomicSuperSignature/actions)
