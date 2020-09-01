##### PCAmodel from 1,399 studies ##############################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_1399/refinebioRseq_PCAmodel_hclust.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_1399/gsea")
#
# library(PCAGenomicSignatures)
#
# for (i in seq_len(ncol(PCAmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   res <- msigdb_gsea(i, PCAmodel, category = "C2")
#   saveRDS(res, fpath)
# }



##### PCAmodel from 536 studies ################################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536/gsea")
#
# library(PCAGenomicSignatures)
# library(PLIER)
# library(clusterProfiler)
#
# data(canonicalPathways)
# data(bloodCellMarkersIRISDMAP)
# data(svmMarkers)
# allPaths <- combinePaths(canonicalPathways, bloodCellMarkersIRISDMAP, svmMarkers)
#
# source('~/data2/GenomicSuperSignature/R/gmtToMatrix.R')
# term2gene <- matrixToTERM2GENE(allPaths)
#
# for (i in seq_len(ncol(PCAmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   geneList <- model(PCAmodel)[,i]
#   geneList <- sort(geneList, decreasing = TRUE)
#   res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
#   saveRDS(res, fpath)
# }



# ##### PCAmodel from 536 studies with allGenes ##################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536_allGenes/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536_allGenes/gsea")
#
# library(PCAGenomicSignatures)
# library(PLIER)
# library(clusterProfiler)
#
# data(canonicalPathways)
# data(bloodCellMarkersIRISDMAP)
# data(svmMarkers)
# allPaths <- combinePaths(canonicalPathways, bloodCellMarkersIRISDMAP, svmMarkers)
#
# source('~/data2/GenomicSuperSignature/R/gmtToMatrix.R')
# term2gene <- matrixToTERM2GENE(allPaths)
#
# for (i in seq_len(ncol(PCAmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   geneList <- model(PCAmodel)[,i]
#   geneList <- sort(geneList, decreasing = TRUE)
#   res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
#   saveRDS(res, fpath)
# }



# ##### PCAmodel from 536 studies C2.CP ##########################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536/gsea_c2cp")
#
# if (!dir.exists(out_dir)) {dir.create(out_dir)}
#
# library(PCAGenomicSignatures)
# library(clusterProfiler)
#
# # MSigDB C2 CP
# term2gene <- clusterProfiler::read.gmt("~/data2/Genomic_Super_Signature/GSEA/data/c2.cp.v7.1.symbols.gmt")
# colnames(term2gene) <- c("gs_name", "entrez_gene")
#
# for (i in seq_len(ncol(PCAmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   geneList <- model(PCAmodel)[,i]
#   geneList <- sort(geneList, decreasing = TRUE)
#   res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
#   saveRDS(res, fpath)
# }



##### PCAmodel from 536 studies C2 #############################################
dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536/refinebioRseq_PCAmodel.rds"))
out_dir <- file.path(dat_dir, "PCAmodel_536/gsea_c2")

if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(PCAGenomicSignatures)
library(clusterProfiler)

# MSigDB C2
term2gene <- clusterProfiler::read.gmt("~/data2/Genomic_Super_Signature/GSEA/data/c2.all.v7.1.symbols.gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(PCAmodel))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path(out_dir, fname)

  geneList <- model(PCAmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}
