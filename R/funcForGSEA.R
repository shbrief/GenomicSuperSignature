#' MSigDB GSEA results
#'
#' @param ind An interger. Index of avgLoading/PCcluster to apply GSEA.
#' @param PCAmodel PCAGenomicSignature object.
#' @param category A character vector representing MSigDB category. Options are
#' "H", "C1", "C2"(default), "C3", "C4", "C5", "C6", and "C7"
#' @param n An interger. The number of top and bottom enriched pathways to plot.
#' Default is \code{NULL} and print out all the pathways enriched under \code{pvalueCutoff}.
#' @param pvalueCutoff Cutoff for both pvalue and p.adjust. Default is 0.5.
#' @param minGSSize A mininum size of gene set to be analyzed
#' @param maxGSSize A maximum size of gene set to be analyzed
#' @param pAdjustMethod p-value adjustment methods, which will be used as an input
#' for \code{method} argument of \code{\link[stats]{p.adjust}} function. Available
#' options are "holm", "hochberg", "hommel", "bonferroni", "BH"(default), "BY", "fdr", "none".
#' @param verbose Logical. Default is FALSE.
#' @param seed Logical. Default is FALSE.
#' @param by Available options are \code{c("fgsea", "DOSE")}. Default is "fgsea".
#' @param geneSets Custom genesets to use with MSigDB genesets. It should be in
#' a named list format.
#'
#' @return Barplot of GSEA output. Top and bottom \code{n} genesets based on NES
#' are plotted and qvalues are denoted by color.
#'
#' @export
msigdb_gsea <- function(ind, PCAmodel, category = "C2", n = NULL, pvalueCutoff = 0.5,
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        verbose = FALSE, seed = FALSE, by = "fgsea", geneSets = NULL) {

    ## Target geneList
    al <- model(PCAmodel)[, ind]
    obj <- list()
    obj[[1]] <- names(al)
    names(al) <- EnrichmentBrowser::idMap(obj, "hsa", from = "SYMBOL", to = "ENTREZID")[[1]]
    al <- sort(al, decreasing = TRUE)

    ## Formating
    geneList <- al
    gene <- names(geneList)[abs(geneList) > mean(abs(geneList))]

    ## MSigDB
    m <- msigdbr::msigdbr(species = "Homo sapiens", category = category) %>%
        clusterProfiler.dplyr::select(gs_name, entrez_gene)

    ## Custom genesets
    if (!is.null(geneSets)) {
        if (is.list(geneSets)) {
            geneSets <- EnrichmentBrowser::idMap(geneSets, "hsa", from = "SYMBOL", to = "ENTREZID")
            custom_m <- reshape2::melt(geneSets)
            colnames(custom_m) <- c("entrez_gene", "gs_name")
            custom_m$gs_name <- gsub(" ", "_", custom_m$gs_name)
        }
        m <- rbind(m, custom_m)
    }

    ## GSEA
    ms <- clusterProfiler::GSEA(geneList, TERM2GENE = m, pvalueCutoff = pvalueCutoff,
                                minGSSize = minGSSize, maxGSSize = maxGSSize,
                                pAdjustMethod = pAdjustMethod, verbose = verbose,
                                seed = seed, by = by)

    ## Subset
    y <- clusterProfiler.dplyr::mutate(ms, ordering = abs(NES)) %>%
        clusterProfiler.dplyr::arrange(dplyr::desc(ordering))

    if (is.null(n)) {
        y <- clusterProfiler.dplyr::group_by(y, sign(NES))
        return(y)
    } else {
        y <- clusterProfiler.dplyr::group_by(y, sign(NES)) %>%
            clusterProfiler.dplyr::slice(1:n)
        return(y)
    }
}


#' Subset GSEA output
#'
#'
#' @export
subsetGSEA <- function(gseaRes, n = 20) {

    topPathways <- clusterProfiler.dplyr::mutate(gseaRes, ordering = abs(NES)) %>%
        clusterProfiler.dplyr::arrange(dplyr::desc(ordering)) %>%
        clusterProfiler.dplyr::group_by(sign(NES)) %>%
        clusterProfiler.dplyr::slice(1:n)
    return(topPathways)
}



#' @title Order genes in loading vectors
#' @description This function takes Z matrix (= average loadings) and orders the
#' genes in each loading vector (= PCcluster) in a descending manner.
#'
#' @param LoadingMatrix An avgloading matrix. Rows represent genes and columns
#' represent clusters of principle components
#' @param LoadingVector A list of column names or indexes of \code{LoadingMatrix}
#' you want to check. Default is \code{NULL}, under which the function takes the
#' all column names of \code{LoadingMatrix}
#' @param abs Under the defaul condition (\code{TRUE}), this function will
#' create a gene list based on the absolute value.
#'
#' @return A list of loadings selected by \code{LoadingVector}, where all the
#' genes in each loading are listed in descending order.
#'
#' @export
makeGeneList <- function(LoadingMatrix, LoadingVector = NULL, abs = TRUE) {

    if (is.null(LoadingVector)) {
        LoadingVector <- colnames(LoadingMatrix)
    }

    if (abs == TRUE) {
        LoadingMatrix <- abs(LoadingMatrix)
    }

    geneLists <- list()
    for (x in LoadingVector) {
        if (x %in% colnames(LoadingMatrix)) {
            geneList <- LoadingMatrix[, x]    # feature 1: numeric vector
            names(geneList) <- rownames(LoadingMatrix)    # feature 2: named vector
        } else {
            geneList <- LoadingMatrix[, c(x)]   # If index, not the name of column, is provided.
            names(geneList) <- rownames(LoadingMatrix)
        }

        geneList <- sort(geneList, decreasing = TRUE)    # feature 3: decreasing order
        geneLists[[as.character(x)]] <- geneList
    }
    return(geneLists)
}


#' @title GSEA on pre-ordered gene lists
#' @description This function is a wrapper of \code{\link[clusterProfiler]{GSEA}} function,
#' making it applicable to a list of gene lists. Set seed for reproducible result.
#'
#' @importFrom clusterProfiler GSEA
#'
#' @param geneList A list of genes ordered by rank
#' @param TERM2GENE User input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
#' @param TERM2NAME User input of TERM TO NAME mapping, a data.frame of 2 column with term and name. Optional.
#' @param minGSSize A mininum size of gene set to be analyzed
#' @param maxGSSize A maximum size of gene set to be analyzed
#' @param pvalueCutoff p-value cutoff
#' @param verbose Logical. Default is \code{FALSE}
#' @param ... Any additional argument inherited from \code{\link[clusterProfiler]{GSEA}}.
#'
#' @return A list of \code{gseaResult} objects. \code{NA} if there is no enrichment result.
#'
#' @export
run_gsea <- function(geneList, TERM2GENE, TERM2NAME,
                     minGSSize = 10, maxGSSize = 500,
                     pvalueCutoff = 0.05, verbose = FALSE, ...) {
    gsea <- list()
    for (x in names(geneList)) {
        res <- clusterProfiler::GSEA(geneList[[x]],
                                    TERM2GENE = TERM2GENE,
                                    TERM2NAME = TERM2NAME,
                                    minGSSize = minGSSize,
                                    maxGSSize = maxGSSize,
                                    verbose = verbose,
                                    pvalueCutoff = pvalueCutoff,
                                    ...)

        # collect GSEA outputs with enrichment result under the assigned pvalueCutoff
        if (nrow(res) != 0) {
            gsea[[x]] <- res
        }
    }
    return(gsea)
}


#' Subset enriched pathways of loading vectors
#'
#' This function is renamed from \code{topPathways} to \code{subsetPathways}.
#'
#' @importFrom enrichplot cnetplot
#' @import methods
#'
#' @param PCAmodel PCAGenomicSignatures object. Also an output from \code{\link[clusterProfiler]{GSEA}} can be used.
#' @param ind A numeric vector containing the PCcluster number you want to check
#' enriched pathways. If not specified, this function returns results from all the PCclusters.
#' @param n The number of top and bottom pathways to be selected based on normalized
#' enrichment score (NES).
#' @param both Default is \code{FALSE}, where only the top \code{n} pathways will
#' be printed. If it is set to \code{TRUE}, the ouput will contain both top and
#' bottom \code{n} pathways.
#'
#'
#' @return A data frame with top and bottom \code{n} pathways from the enrichment results.
#'
#' @export
subsetPathways <- function(PCAmodel, ind = NULL, n = 10, both = FALSE) {

    if (is(PCAmodel, "PCAGenomicSignatures")) {
        gsea_loading <- unlist(gsea(PCAmodel))
    } else if (class(PCAmodel) %in% c("data.frame", "matrix", "list")) {
        gsea_loading <- PCAmodel
    }

    res <- list()
    for (name in names(gsea_loading)) {
        x <- gsea_loading[[name]]
        up <- x$Description[order(x$NES, decreasing=TRUE)][1:n]
        down <- x$Description[order(x$NES, decreasing=FALSE)][1:n]
        res[[name]] <- c(up, down)

        ### Adding NES and qvalues?
        # y <- x[,c("Description", "NES", "qvalues")]
        # up <- y[order(y$NES, decreasing=TRUE),][1:n,]
        # down <- y[order(y$NES, decreasing=FALSE),][1:n,]
        # res[[name]] <- rbind(up, down)
    }

    res <- data.frame(res, check.names = FALSE)
    rnames <- c(paste0("Up_", 1:n), paste0("Down_", 1:n))
    rownames(res) <- rnames

    if (both == FALSE) {res <- res[1:n,]}

    if (is.null(ind)) {
        res <- S4Vectors::DataFrame(res)
        return(res)
    } else {
        col_num <- which(colnames(res) %in% c(paste0("PCcluster", ind)))
        res <- res[, col_num, drop = FALSE]
        res <- S4Vectors::DataFrame(res)
        return(res)
    }
}
