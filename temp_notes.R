## Update 11.01.20 : PCclsuter --> RAV
names(metadata(PCAmodel)$size) <- gsub("PCcluster", "RAV", names(metadata(PCAmodel)$size))
colnames(PCAmodel) <- gsub("PCcluster", "RAV", colnames(PCAmodel))
names(colData(PCAmodel))[1] <- "RAV"
names(colData(PCAmodel)$gsea) <- gsub("PCcluster", "RAV", names(colData(PCAmodel)$gsea))
save(PCAmodel, file = "~/data2/PCAGenomicSignatures/data/miniPCAmodel.RData")


crc_dir <- "~/data2/PCAGenomicSignaturesPaper/Results/CRC"
load(file.path(crc_dir, "data/eSets/setNames.RData"))
set <- setNames[18]    # manually changed
load(paste0(crc_dir, "/data/eSets_new/", set, '.RData'))
eset.tmp <- get(set)
colnames(pData(eset.tmp)) <- gsub("PCcluster", "RAV", colnames(pData(eset.tmp)))
eset.tmp -> NHS.HPFS_eset
save(NHS.HPFS_eset, file = paste0(crc_dir, "/data/eSets_new/", set, '.RData'))
