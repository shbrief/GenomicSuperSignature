#' Barplot GSEA output
#'
#' @param ind An interger. Index of avgLoading/PCcluster to apply GSEA.
#' @param PCAmodel PCAGenomicSignature object.
#' @param category A character vector representing MSigDB category. Options are
#' "H","C1","C2"(default), "C3","C4","C5","C6", and "C7"
#' @param n An interger. The number of top and bottom enriched pathways to plot. Default is 10.
#'
#' @return Barplot of GSEA output. Top and bottom \code{n} genesets are plotted.
#'
#' @export
msigdb_gsea <- function(ind, PCAmodel, category = "C2", n = 10) {

    ## Target geneList
    al <- model(PCAmodel)[, ind]
    names(al) <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(al), column='ENTREZID', keytype='SYMBOL')
    al <- sort(al, decreasing = TRUE)

    ## Formating
    geneList <- al
    gene <- names(geneList)[abs(geneList) > mean(abs(geneList))]

    ## MSigDB
    m <- msigdbr::msigdbr(species = "Homo sapiens", category = category) %>%
        clusterProfiler.dplyr::select(gs_name, entrez_gene)
    ms <- clusterProfiler::GSEA(geneList, TERM2GENE = m)

    ## Barplot
    y <- clusterProfiler.dplyr::mutate(ms, ordering = abs(NES)) %>% clusterProfiler.dplyr::arrange(desc(ordering))
    y_bar <- clusterProfiler.dplyr::group_by(y, sign(NES)) %>% clusterProfiler.dplyr::slice(1:n)
    ggplot(y_bar, aes(NES, forcats::fct_reorder(Description, NES), fill = qvalues), showCategory=(n*2)) +
        geom_bar(stat='identity') +
        scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
        theme_minimal() + ylab(NULL)
}




#' Dotplot of Gene Ontology Enrichment result
#'
#' @param x Available options are "GeneRatio" (default) and "Count"
#' @param color "Available options are "pvalue", "p.adjust" (default), and "qvalue"
#' @export
do_dotplot <- function(ind, PCAmodel, x = "GeneRatio", color = "p.adjust") {

    ## Target geneList
    al <- model(PCAmodel)[, ind]
    names(al) <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(al), column='ENTREZID', keytype='SYMBOL')
    al <- sort(al, decreasing = TRUE)

    ## Formating
    geneList <- al
    gene <- names(geneList)[abs(geneList) > mean(abs(geneList))]

    edo2 <- DOSE::gseDO(geneList)
    dotplot(edo2, showCategory = 10, font.size = 8,
            x = x, color = color)
}
