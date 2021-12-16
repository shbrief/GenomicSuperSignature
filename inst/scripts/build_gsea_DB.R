#' Choose the correct GSEA script based on the RAVmodel and the annotation database
#' you are interested in. Unless specified, RAVmodels are built with top 90% varying
#' common genes. Here is the current list:
#'
#' 1. 1399 studies + MSigDB C2
#' 2. 536 studies + PLIER priors
#' 3. 536 studies + all genes + PLIER priors
#' 4. 536 studies + MSigDB C2 CP
#' 5. 536 studies + MSigDB C2 all






# ##### 1. RAVmodel from 1,399 studies #########################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# RAVmodel <- readRDS(file.path(dat_dir, "RAVmodel_1399/refinebioRseq_RAVmodel_hclust.rds"))
# out_dir <- file.path(dat_dir, "RAVmodel_1399/gsea")
#
# library(GenomicSuperSignatures)
#
# for (i in seq_len(ncol(RAVmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   res <- msigdb_gsea(i, RAVmodel, category = "C2")
#   saveRDS(res, fpath)
# }



# ##### 2. RAVmodel from 536 studies ###########################################
# dat_dir <- "~/data2/GenomicSuperSignatureLibrary/refinebioRseq"
# RAVmodel <- readRDS(file.path(dat_dir, "RAVmodel_536/refinebioRseq_RAVmodel.rds"))
# out_dir <- file.path(dat_dir, "RAVmodel_536/gsea_PLIERpriors")
#
# library(GenomicSuperSignatures)
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
# for (i in seq_len(ncol(RAVmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   geneList <- RAVindex(RAVmodel)[,i]
#   geneList <- sort(geneList, decreasing = TRUE)
#   res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
#   saveRDS(res, fpath)
# }



# ##### 3. RAVmodel from 536 studies PLIERpriors with allGenes #################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# RAVmodel <- readRDS(file.path(dat_dir, "RAVmodel_536_allGenes/refinebioRseq_RAVmodel.rds"))
# out_dir <- file.path(dat_dir, "RAVmodel_536_allGenes/gsea_PLIERpriors")
#
# library(GenomicSuperSignatures)
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
# for (i in seq_len(ncol(RAVmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   geneList <- RAVindx(RAVmodel)[,i]
#   geneList <- sort(geneList, decreasing = TRUE)
#   res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
#   saveRDS(res, fpath)
# }



# ##### 4. RAVmodel from 536 studies C2.CP #####################################
# dat_dir <- "~/data2/PCAGenomicSignatureLibrary/refinebioRseq"
# RAVmodel <- readRDS(file.path(dat_dir, "RAVmodel_536/refinebioRseq_RAVmodel.rds"))
# out_dir <- file.path(dat_dir, "RAVmodel_536/gsea_C2CP")
#
# if (!dir.exists(out_dir)) {dir.create(out_dir)}
#
# library(GenomicSuperSignatures)
# library(clusterProfiler)
#
# # MSigDB C2 CP
# term2gene <- clusterProfiler::read.gmt("~/data2/Genomic_Super_Signature/GSEA/data/c2.cp.v7.1.symbols.gmt")
# colnames(term2gene) <- c("gs_name", "entrez_gene")
#
# for (i in seq_len(ncol(RAVmodel))) {
#   fname <- paste0("gsea_", i, ".rds")
#   fpath <- file.path(out_dir, fname)
#
#   geneList <- RAVindx(RAVmodel)[,i]
#   geneList <- sort(geneList, decreasing = TRUE)
#   res <- GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
#   saveRDS(res, fpath)
# }



#### 5. RAVmodel from 536 studies C2 ###########################################
dat_dir <- "~/data2/GenomicSuperSignatureLibrary/refinebioRseq"
# RAVmodel <- readRDS(file.path(dat_dir, "RAVmodel_536/refinebioRseq_RAVmodel.rds"))
# out_dir <- file.path(dat_dir, "RAVmodel_536/gsea_c2")
RAVmodel <- readRDS(file.path(dat_dir, "RAVmodel_536/refinebioRseq_RAVmodel_20211207.rds"))
out_dir <- file.path(dat_dir, "RAVmodel_536_noLINC/gsea_c2")

if (!dir.exists(out_dir)) {dir.create(out_dir)}

library(GenomicSuperSignatures)
library(clusterProfiler)

# MSigDB C2
term2gene <- clusterProfiler::read.gmt("~/data2/[archive]Genomic_Super_Signature/GSEA/data/c2.all.v7.1.symbols.gmt")
colnames(term2gene) <- c("gs_name", "entrez_gene")

for (i in seq_len(ncol(RAVmodel))) {
  fname <- paste0("gsea_", i, ".rds")
  fpath <- file.path(out_dir, fname)

  geneList <- RAVindex(RAVmodel)[,i]
  geneList <- sort(geneList, decreasing = TRUE)
  res <- clusterProfiler::GSEA(geneList, TERM2GENE = term2gene, pvalueCutoff = 0.05)
  saveRDS(res, fpath)
}
