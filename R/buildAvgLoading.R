#' Map adaptive clusters to ordered, numeric cluster number
#'
#' This function takes a named character vector and maps its values to cluster
#' integers based on a lookup table containing `raw_cl_num` and `cluster` columns.
#'
#' @param char_vector A named character vector from adaptive clustering where 
#' values are cluster identifiers(e.g., "013.001", "012.011") that need to be 
#' mapped to cluster integers.
#' @param mapping_table A data frame containing at least two columns:
#'   \itemize{
#'     \item raw_cl_num: Character values matching those in char_vector
#'     \item cluster: Integer cluster numbers to map to
#'   }
#'
#' @return A named integer vector with the same names as char_vector but with
#' values replaced by their corresponding cluster numbers. Values not found in 
#' the mapping table will be set to NA.
#'
#' @examples
#' # Create example named vector
#' my_vector <- c(
#'   DRP000987.PC1 = "013.001",
#'   DRP000987.PC2 = "012.011",
#'   DRP000987.PC3 = "002.016"
#' )
#' 
#' # Create mapping table
#' mapping <- data.frame(
#'   raw_cl_num = c("001.001", "013.001", "012.011", "002.016"),
#'   Freq = c(39, 25, 22, 20),
#'   cluster = c(1, 2, 3, 4)
#' )
#' 
#' # Map values to clusters
#' result <- map_to_clusters(my_vector, mapping)
#' print(result)
#' # DRP000987.PC1 DRP000987.PC2 DRP000987.PC3 
#' #             2             3             4
#'
.mapToClusters <- function(char_vector, mapping_table) {
    # Create a named vector for fast lookup
    lookup <- setNames(mapping_table$cluster, mapping_table$raw_cl_num)
    
    # Map the values
    result <- lookup[char_vector]
    
    # Preserve the original names
    names(result) <- names(char_vector)
    
    return(result)
}


#' Extract Principal Component Counts by Study
#'
#' This function identifies columns containing principal components (PCs) in a dataset,
#' extracts the study names, and counts the number of PCs for each study.
#'
#' @param data A data frame or matrix containing columns with names in the format
#'   "studyName.PC#" where # is a number (e.g., "GSE26682.GPL96_eset.PC1").
#'
#' @return A named integer vector where names are study names (part before ".PC#")
#'   and values are the total number of PCs per study.
#'
#' @examples
#' # Create example data
#' example_data <- data.frame(
#'   GSE26682.GPL96_eset.PC1 = rnorm(10),
#'   GSE26682.GPL96_eset.PC2 = rnorm(10),
#'   NHS.HPFS_eset.PC1 = rnorm(10)
#' )
#' 
#' # Extract PC counts
#' result <- extract_pc_counts(example_data)
#' print(result)
#' # GSE26682.GPL96_eset      NHS.HPFS_eset 
#' #                  2                  1
#'
.extractPCCounts <- function(data) {
    # Get column names that contain ".PC"
    pc_columns <- grep("\\.PC\\d+$", colnames(data), value = TRUE)
    
    # Extract study names (everything before .PC)
    study_names <- sub("\\.PC\\d+$", "", pc_columns)
    
    # Count PCs per study and return as named vector
    pc_counts <- as.vector(table(study_names))
    names(pc_counts) <- names(table(study_names))
    
    return(pc_counts)
}



#' Calculate Silhouette Information of RAVs
#'
#' @description The silhouette value is a measure of how similar an object is
#' to its own cluster (cohesion) compared to other clusters (separation). The
#' silhouette width ranges from -1 to +1, where a high value indicates that
#' the object is well matched to its own cluster and poorly matched to
#' neighboring clusters.
#'
#' @param dat A matrix with all the top PCs from training data to be clustered.
#' @param kmeansRes Output from \code{stats::kmeans}.
#'
#' @return Silhouette-class object, which is an n x 3 matrix with attributes.
#' 
#' @keywords internal
#'
#' @seealso \code{\link[stats]{kmeans}}
#'
.calculateSilhouetteWidth <- function (dat, kmeansRes) {
    swRes <- cluster::silhouette(x = kmeansRes$cluster,
                                 dist = cluster::daisy(dat))
    SW <- summary(swRes)
    return(SW)
}



#' Calculate average loadings of each cluster
#'
#' @param dat A data frame. Each row represents principal components from
#' different training datasets. Columns are genes used for PCA analysis.
#' @param cluster Provide pre-defined cluster membership of your data.
#' @param study Under default (\code{TRUE}), studies involved in each cluster
#' will be added to the output.
#'
#' @return
#' A named list of 6 elements is returned. It contains:
#' \describe{
#'    \item{\code{cluster}}{A numeric vector on cluster membership of PCs}
#'    \item{\code{size}}{A integer vector on the size of clusters}
#'    \item{\code{avgLoading}}{A matrix of average loadings. Columns for
#'    clusters and rows for genes}
#'    \item{\code{k}}{The number of clusters}
#'    \item{\code{n}}{The number of top PCs used for clustering}
#'    \item{\code{studies}}{A list of character vector containing studies in
#'    each cluster}
#' }
#'
#' @examples
#' data(miniAllZ)
#' data(res_hcut)
#' res <- buildAvgLoading(miniAllZ, cluster = res_hcut$cluster)
#'
#' @export
buildAvgLoading <- function(dat, cluster = NULL, study = TRUE) {

    # Input validation
    if (!is.null(cluster)) {
        k <- length(unique(cluster))
        x <- table(cluster) %>% # assign re-ordered column name
            as.data.frame() %>%
            dplyr::rename(raw_cl_num = cluster) %>%
            dplyr::mutate(cluster = row_number())
        cluster <- .mapToClusters(cluster, x)
            
        res <- list(cluster = cluster, size = x$Freq)
    } else {
        stop("Cluster membership of elements should be provided 
             through 'cluster' argument.")
    }
    stopifnot(length(study) == 1L, !is.na(study), is.logical(study))

    # Extract the number of PCs per study
    pcnum <- .extractPCCounts(dat)
    if (length(unique(pcnum)) == 1) {
        n <- unique(pcnum)
    } else {n <- pcnum}

    # Separate the PC table into each cluster
    cl_ls <- vector(mode = "list", length = k)
    for (i in seq_len(k)) {
        datName <- paste0("Cl", k, "_",
                          formatC(i, width = 2, format = "d", flag = "0"))
        cl_ls[[i]] <- dat[, res$cluster==i, drop=FALSE] %>% t
        names(cl_ls)[i] <- datName
    }

    # the number of unique datasets in each cluster
    unique_sets <- vector(length = k)
    for (i in seq_len(k)) {
        dat <- cl_ls[[i]]
        dataSetName <- gsub(".PC\\d+$", "", rownames(dat))
        uniqueDataSetName <- length(unique(dataSetName))
        unique_sets[i] <- uniqueDataSetName
    }

    # Calculate the average of loadings in each cluster
    names(cl_ls) <- paste0(names(cl_ls), " (", res$size, "/", unique_sets, ")")
    l <- ncol(dat)   # the number of genes for average loading
    avg.loadings <- vapply(cl_ls, colMeans, FUN.VALUE = numeric(l))
    
    # Save results
    res$avgLoading <- as.matrix(avg.loadings) # avgLoading
    res$k <- k # number of clusters
    res$n <- n # number of PCs per study
    if (study) {res$studies <- findStudiesInCluster(res)} # studies in cluster

    # # Silhouette Width <<<<<<<<< now done in `06_Clustering.R`
    # sw <- .calculateSilhouetteWidth(dat, res)
    # res$sw <- sw

    return(res)
}

