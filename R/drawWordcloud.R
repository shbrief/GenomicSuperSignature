#' Extract the list of PCs in a cluster
#'
#' @param x A PCAGenomicSignatures object
#' @param ind An index of PCcluster
#'
#' @export
PCinCluster <- function(x, ind) {
    cluster <- S4Vectors::metadata(x)$cluster
    k <- which(cluster == ind)
    out <- names(k)
    return(out)
}


#' Build a two-column word/frequency table
#'
#' @import dplyr
#'
#' @param x A PCAGenomicSignatures object
#' @param ind An index of PCcluster
#' @param rm.noise An integer. Any MeSH term found less than the given value here
#' will be excluded from wordcloud. If \code{rm.noise = 0}, all the MeSH terms
#' in PCcluster will be used.
#' @param weighted A logical. If \code{TRUE}, MeSH terms from each study are
#' weighted based on the variance explained by the principle component of the
#' study contributing a give PCcluster.
#'
#' @return A table with two columns, \code{word} and \code{freq}. MeSH terms in
#' the defined PCcluster (by \code{ind} argument) is ordered based on their frequency.
#'
#' @export
meshTable <- function(x, ind, rm.noise, weighted) {

    ### Create a 'universe' for bag-of-words model
    bow <- metadata(x)$MeSH_freq %>% unlist  # frequency of the `name` in the background
    bow <- bow[which(bow > rm.noise)]   # remove rare terms

    ### Variance explained by PC
    if (weighted == FALSE) {
        study_id <- studies(x)[[ind]]   # a list of studies in PCcluster
        recount2_MeSH <- mesh(x)   # all the MeSH data
        mesh_subset <- recount2_MeSH[study_id]   # subset to the participating studies

        ### Combine all MeSH words
        d <- list()
        for (i in seq_along(mesh_subset)) {
            d <- c(d, mesh_subset[[i]]$name)
        }

        ### Build a term-frequency table
        summary <- unlist(d) %>% table(.)
        summary <- summary[names(summary) %in% names(bow)]

    ### Weighted, counting on variance explained by
    } else {
        PCs <- PCinCluster(x, ind)
        varAll <- Reduce(cbind, PCAsummary(x))
        var <- varAll[,PCs,drop=FALSE]
        study_id <- gsub("\\.PC.*$", "", PCs)
        recount2_MeSH <- mesh(x)   # all the MeSH data
        mesh_subset <- sapply(recount2_MeSH[study_id], function(x) {x$name})

        weight <- as.data.frame(matrix(nrow = 0, ncol = 3))
        colnames(weight) <- c("PCs", "var", "mesh")

        for (i in seq_along(PCs)) {
            new <- data.frame(PCs = PCs[i],
                             var = var["Variance", i],
                             mesh = mesh_subset[[i]])
            weight <- rbind(weight, new)
        }

        weight <- weight %>% group_by(mesh) %>% summarise(total_var = sum(var))

        ### Build a term-frequency table
        summary <- setNames(weight$total_var, as.character(weight$mesh))
        summary <- summary[names(summary) %in% names(bow)]
    }

    for (i in seq_along(summary)) {
        summary[i] <- summary[i]/bow[names(summary[i])]
    }
    all <- as.data.frame(matrix(NA, ncol = 0, nrow = length(summary)))
    all$word <- names(summary)
    all$freq <- summary
    all <- all[order(all$freq, decreasing = TRUE),]

    return(all)
}


#' @title Draw wordcloud using the collection of PCclusters' MeSH terms
#' @description Plot a word cloud using the remaining MeSH terms in the selected
#' PCcluster after user-defined filtering.
#'
#' @importFrom wordcloud wordcloud
#'
#' @param x PCAGenomicSignatures object
#' @param ind An index of the PCcluster you want to draw wordcloud.
#' @param rm.noise An integer. Under the default condition (\code{rm.noise=NULL}),
#' if cluster size (= \code{s}) is smaller than 8, \code{rm.noise = floor(s*0.5)}.
#' For clusters with >= 8 PCs, \code{rm.noise = 4}. If \code{rm.noise = 0}, all
#' the MeSH terms in PCcluster will be used to draw wordcloud.
#' @param scale A \code{scale} argument for \code{\link[wordcloud]{wordcloud}} function
#' @param weighted A logical. If \code{TRUE} (default), MeSH terms from each study are
#' weighted based on the variance explained by the principle component of the
#' study contributing a give PCcluster.
#' @param seed Random seed. If it is not specified, \code{set.seed(1234)} will be used.
#'
#' @return A word cloud
#'
#' @export
drawWordcloud <- function(x, ind, rm.noise = NULL, scale = c(3, 0.5),
                         weighted = TRUE, seed = NULL) {

    if (is.null(rm.noise)) {
        s <- metadata(x)$size[ind]
        if (s < 8) {rm.noise = floor(s*0.5)}
        else if (s >= 8) {rm.noise = 4}

        ## Minimum rm.noise version
        # s <- metadata(x)$size[ind]
        # rm.noise = floor(s*0.2)
        # if (rm.noise > 4) {rm.noise = 4}
        # else if (rm.noise == 0) {rm.noise = 1}
    }

    # MeSH word table
    all <- meshTable(x, ind, rm.noise = rm.noise, weighted = weighted)

    # generate the word cloud
    if (is.null(seed)) {set.seed(1234)}
    wordcloud(words = all$word, freq = all$freq, scale = scale,
              max.words = Inf, random.order = FALSE,
              colors = brewer.pal(8, "Dark2"))
}
