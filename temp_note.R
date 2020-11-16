## Update 'PCcluster' to 'RAV'
data.dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq/PCAmodel_536_clNum4"
PCAmodel <- readRDS(file.path(data.dir, "refinebioRseq_PCAmodel_C2.rds"))
names(metadata(PCAmodel)$size) <- gsub("PCcluster", "RAV", names(metadata(PCAmodel)$size))
names(metadata(PCAmodel)$size) <- paste0("RAV", 1:length(metadata(PCAmodel)$size))
colnames(PCAmodel) <- gsub("PCcluster", "RAV", colnames(PCAmodel))
names(gsea(PCAmodel)) <- gsub("PCcluster", "RAV", names(gsea(PCAmodel)))
names(colData(PCAmodel))[1] <- "RAV"
saveRDS(PCAmodel, file.path(data.dir, "refinebioRseq_PCAmodel_C2.rds"))



## Update package name from 'PCAGenomicSignatures' to "GenomicSuperSignature'
fname <- "~/data2/PCAGenomicSignatures/vignettes/PCAmodel_C2.rds"
PCAmodel <- readRDS(fname)
RVmodel <- GenomicSuperSignature::PCAGenomicSignatures(assays = list(model = model(PCAmodel)))
metadata(RVmodel) <- GenomicSuperSignature::metadata(PCAmodel)
metadata(RVmodel)$version <- ">= 0.0.6"
studies(RVmodel) <- GenomicSuperSignature::studies(PCAmodel)
silhouetteWidth(RVmodel) <- GenomicSuperSignature::silhouetteWidth(PCAmodel)
trainingData(RVmodel) <- GenomicSuperSignature::trainingData(PCAmodel)
gsea(RVmodel) <- GenomicSuperSignature::gsea(PCAmodel)
saveRDS(RVmodel, "~/data2/PCAGenomicSignatures/vignettes/RVmodel_C2.rds")


## Update 'RAV' to 'RV'
# fname <- "~/data2/PCAGenomicSignatures/vignettes/RVmodel_C2.rds"
# RVmodel <- readRDS(fname)
# names(metadata(RVmodel)$size) <- paste0("RV", 1:length(metadata(RVmodel)$size))
# colnames(RVmodel) <- gsub("RAV", "RV", colnames(RVmodel))
# names(gsea(RVmodel)) <- gsub("RAV", "RV", names(gsea(RVmodel)))
# names(colData(RVmodel))[1] <- "RV"
# saveRDS(RVmodel, "~/data2/PCAGenomicSignatures/vignettes/RVmodel_C2.rds")
fname <- "~/data2/PCAGenomicSignatures/vignettes/RVmodel_C2.rds"
RVmodel <- readRDS(fname)
names(metadata(RVmodel)$size) <- paste0("RAV", 1:length(metadata(RVmodel)$size))
colnames(RVmodel) <- gsub("RV", "RAV", colnames(RVmodel))
names(gsea(RVmodel)) <- gsub("RV", "RAV", names(gsea(RVmodel)))
names(colData(RVmodel))[1] <- "RAV"
saveRDS(RVmodel, "~/data2/PCAGenomicSignatures/vignettes/RVmodel_C2.rds")
