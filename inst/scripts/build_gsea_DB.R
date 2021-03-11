#' Choose the correct GSEA script based on the PCAmodel and the annotation database
#' you are interested in. Unless specified, PCAmodels are built with top 90% varying
#' common genes. Here is the current list:
#'
#' 1. 1399 studies + MSigDB C2
#' 2. 536 studies + PLIER priors
#' 3. 536 studies + all genes + PLIER priors
#' 4. 536 studies + MSigDB C2 CP
#' 5. 536 studies + MSigDB C2 all






# ##### 1. PCAmodel from 1,399 studies #########################################
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



# ##### 2. PCAmodel from 536 studies ###########################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536/gsea_PLIERpriors")
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



# ##### 3. PCAmodel from 536 studies PLIERpriors with allGenes #################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536_allGenes/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536_allGenes/gsea_PLIERpriors")
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



# ##### 4. PCAmodel from 536 studies C2.CP #####################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536/gsea_C2CP")
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



#### 5. PCAmodel from 536 studies C2 ###########################################
dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536/refinebioRseq_PCAmodel.rds"))
# out_dir <- file.path(dat_dir, "PCAmodel_536/gsea_c2")
PCAmodel <- readRDS(file.path(dat_dir, "PCAmodel_536_clNum4/refinebioRseq_PCAmodel.rds"))
out_dir <- file.path(dat_dir, "PCAmodel_536_clNum4/gsea_c2")

if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(PCAGenomicSignatures)
library(clusterProfiler)

# MSigDB C2
term2gene <- clusterProfiler::read.gmt("~/data2/[archive]Genomic_Super_Signature/GSEA/data/c2.all.v7.1.symbols.gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(PCAmodel))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path(out_dir, fname)

  geneList <- model(PCAmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  res <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}
