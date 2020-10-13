#' Barplot GSEA output
#'
#' @param ind An interger. Index of avgLoading/PCcluster to apply GSEA.
#' @param PCAmodel PCAGenomicSignature object.
#' @param category A character vector representing MSigDB category. Options are
#' "H", "C1", "C2"(default), "C3", "C4", "C5", "C6", and "C7"
#' @param n An interger. The number of top and bottom enriched pathways to plot. Default is 10.
#' @param pvalueCutoff Cutoff for both pvalue and p.adjust. Default is 0.5.
#' @param gseaRes An output from \link{msigdb_gsea} function. If this argument
#' is provided, you don't need to provide the above inputs: ind, PCAmodel, category, n, pvalueCutoff.
#'
#' @return Barplot of GSEA output. Top and bottom \code{n} genesets based on NES
#' are plotted and qvalues are denoted by color.
#'
#' @export
gseaBarplot <- function(ind, PCAmodel, category = "C2", n = 10, pvalueCutoff = 0.5, gseaRes = NULL) {

    ## Binding the variables from res locally to the function
    NES <- Description <- qvalues <- NULL

    if (is.null(gseaRes)) {
        gseaRes <- msigdb_gsea(ind, PCAmodel, category = category, n = n, pvalueCutoff = pvalueCutoff)
    } else {gseaRes <- gseaRes}

    if (nrow(gseaRes) == 0) return("No pathway is enriched")   # Handle empty dataframes

    ## Barplot
    ggplot(gseaRes, aes(NES, forcats::fct_reorder(Description, NES), fill = qvalues), showCategory=(n*2)) +
        geom_bar(stat='identity') +
        scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) +
        theme_minimal() + ylab(NULL)
}


#' Calculate jaccard similarity of two sets
#'
#' @param a A vector
#' @param b A vector
#' @return A numerical value
#'
#' @keywords internal
.jaccard_similarity <- function(a, b) {
    length(intersect(a, b)) / length(union(a, b))
}

#' Calculate overlap similarity of two sets
#'
#' @param a A vector
#' @param b A vector
#' @return A numerical value
#'
#' @keywords internal
.overlap_similarity <- function(a, b) {
    length(intersect(a, b)) / min(length(a), length(b))
}

#' Format a string using placeholders
#'
#' @param string A an unformatted string with placeholders
#' @param ... Variables to format placeholders with
#' @return A formatted string
#'
#' @examples
#' \dontrun{
#' format_str("Format with {1} and {2}", "x", "y")
#' }
#'
#' @keywords internal
.format_str <- function(string, ...) {
    args <- list(...)
    for (i in 1:length(args)) {
        pattern <- paste("\\{", i, "}", sep="")
        replacement <- args[[i]]
        string <- gsub(pattern, replacement, string)
    }
    return(string)
}


#' Plot the network of enriched pathways
#'
#' @param ind An interger. Index of avgLoading/PCcluster to apply GSEA.
#' @param PCAmodel PCAGenomicSignature object.
#' @param category A character vector representing MSigDB category. Options are
#' "H", "C1", "C2"(default), "C3", "C4", "C5", "C6", and "C7"
#' @param n An interger. The number of top and bottom enriched pathways to plot. Default is 10.
#' @param pvalueCutoff Cutoff for both pvalue and p.adjust. Default is 0.5.
#' @param gseaRes An output from \link{msigdb_gsea} function. If this argument
#' is provided, you don't need to provide the above inputs: ind, PCAmodel, category, n, pvalueCutoff.
#' @param similarity_metric Metric to calculate geneset similarity. Available values
#' are \code{c("jaccard_similarity", "overlap_similarity")}.
#' @param similarity_cutoff Geneset similarity cutoff. Default is 0.3.
#' @param title Plot title
#' @return A visNetwork object
#'
#' @note Modified from \code{\link[hypeR]{hyp_emap}}
#'
#' @importFrom purrr when
#' @importFrom dplyr filter
#' @importFrom igraph graph.adjacency V
#' @importFrom visNetwork visNetwork visNodes visEdges visOptions visInteraction toVisNetworkData visIgraphLayout
#'
#' @export
gseaNetwork <- function(ind, PCAmodel, category = category, n = n, pvalueCutoff = pvalueCutoff,
                        gseaRes = NULL,
                        similarity_metric = c("jaccard_similarity", "overlap_similarity"),
                        similarity_cutoff = 0.3, title = "") {

    if (is.null(gseaRes)) {
        gseaRes <- msigdb_gsea(ind, PCAmodel, category = category, n = n, pvalueCutoff = pvalueCutoff)
        gseaRes <- as.data.frame(gseaRes)
    } else {
        gseaRes <- as.data.frame(gseaRes)
    }

    if (nrow(gseaRes) == 0) return("No pathway is enriched")   # Handle empty dataframes

    # Geneset similarity matrix
    genesets <- vector(mode = "list", length = nrow(gseaRes))
    names(genesets) <- gseaRes$ID
    for (i in seq_along(genesets)) {
        genes <- stringr::str_split(gseaRes$core_enrichment[i], "/")
        genes <- unlist(genes)
        genesets[[i]] <- as.integer(genes)
    }

    genesets.mat <- sapply(genesets, function(x) {
        sapply(genesets, function(y,x) {
            if (similarity_metric == "jaccard_similarity") .jaccard_similarity(x, y)
            else if (similarity_metric == "overlap_similarity") .overlap_similarity(x, y)
            else stop(.format_str("{1} is an invalid metric", similarity_metric))
        }, x)
    })

    m <- as.matrix(genesets.mat)

    # Sparsity settings
    m[m < similarity_cutoff] <- 0

    # Similarity matrix to weighted network
    inet <- igraph::graph.adjacency(m, mode="undirected", weighted=TRUE, diag=FALSE)

    # igraph to visnet
    vnet <- visNetwork::toVisNetworkData(inet)
    nodes <- vnet$nodes
    edges <- vnet$edges

    # Add edge weights
    edges$value <- vnet$edges$weight

    # Add node scaled sizes based on genset size
    size.scaler <- function(x) (x-min(x))/(max(x)-min(x))*30
    node.sizes <- sapply(igraph::V(inet), function(x) gseaRes[x, "setSize"])
    nodes$size <-  size.scaler(node.sizes)+20

    val <- "qvalues"
    nodes$title <- sapply(igraph::V(inet), function(x) {
        paste(val, gseaRes[x, val], sep=": ")
    })

    # Add node scaled weights based on significance
    weight.scaler <- function(x) (x-max(x))/(min(x)-max(x))
    node.weights <- sapply(igraph::V(inet), function(x) gseaRes[x, val])
    nodes$color.border <- "rgb(0,0,0)"
    nodes$color.highlight <- "rgba(199,0,57,0.9)"
    nodes$color.background <- sapply(weight.scaler(node.weights), function(x) {
        if (is.na(x)) {
            return("rgba(199,0,57,0)")
        } else{
            return(paste("rgba(199,0,57,", round(x, 3), ")", sep=""))
        }
    })

    visNetwork(nodes, edges, main=list(text=title, style="font-family:Helvetica")) %>%
        visNodes(borderWidth=1, borderWidthSelected=0) %>%
        visEdges(color="rgb(88,24,69)") %>%
        visOptions(highlightNearest=TRUE) %>%
        visInteraction(multiselect=TRUE, tooltipDelay=300) %>%
        visIgraphLayout(layout="layout_nicely")
}
