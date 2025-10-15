library(dynamicTreeCut)

#' Recursive hierarchical clustering with size constraints
#' 
#' @param data Data matrix/data frame to cluster
#' @param dist_matrix Distance matrix (optional, will be calculated if NULL)
#' @param hclust_obj Hierarchical clustering object (optional, will be calculated if NULL)
#' @param minCl Minimum cluster size
#' @param maxCl Maximum cluster size
#' @param cluster_prefix Prefix for cluster naming (for internal recursion)
#' @param max_depth Maximum recursion depth to prevent infinite loops
#' @param current_depth Current recursion depth (for internal use)
#' @param pad_width Width for zero-padding cluster numbers (default: 3)
#' 
#' @return Named vector of cluster assignments (names are row names from data)
#' 
#' 
recursiveDynamicCut <- function(data, 
                                  dist_matrix = NULL,
                                  hclust_obj = NULL,
                                  minCl, 
                                  maxCl,
                                  cluster_prefix = "",
                                  max_depth = 10,
                                  current_depth = 0,
                                  pad_width = 3) {
    
    # Safety check for recursion depth
    if (current_depth >= max_depth) {
        warning("Maximum recursion depth reached. Some clusters may exceed maxCl.")
        result <- rep(ifelse(cluster_prefix == "", 
                             sprintf("%0*d", pad_width, 1), 
                             paste0(cluster_prefix, ".", sprintf("%0*d", pad_width, 1))), 
                      nrow(data))
        names(result) <- rownames(data)
        return(result)
    }
    
    # Calculate distance matrix if not provided
    if (is.null(dist_matrix)) {
        dist_matrix <- factoextra::get_dist(data, method = "spearman")
    }
    
    # Calculate hierarchical clustering if not provided
    if (is.null(hclust_obj)) {
        hclust_obj <- hclust(dist_matrix, method = "ward.D")
    }
    
    # Apply dynamic tree cut
    clusters <- cutreeDynamic(hclust_obj, 
                              distM = as.matrix(dist_matrix), 
                              method = "hybrid", 
                              deepSplit = 4,
                              minClusterSize = minCl)
    
    # Create cluster names with prefix and zero-padding
    padded_clusters <- sprintf("%0*d", pad_width, clusters)
    
    if (cluster_prefix == "") {
        cluster_names <- padded_clusters
    } else {
        cluster_names <- paste0(cluster_prefix, ".", padded_clusters)
    }
    
    # Create named vector
    names(cluster_names) <- rownames(data)
    
    # Check which clusters exceed maxCl
    cluster_sizes <- table(cluster_names)
    need_recluster <- names(cluster_sizes)[cluster_sizes > maxCl]
    
    # If no clusters need reclustering, we're done
    if (length(need_recluster) == 0) {
        message(paste("Depth", current_depth, ": All clusters are smaller than maxCl =", maxCl))
        return(cluster_names)
    }
    
    message(paste("Depth", current_depth, ":", length(need_recluster), 
                  "cluster(s) need reclustering (sizes:",
                  paste(cluster_sizes[need_recluster], collapse = ", "), ")"))
    
    # Recluster each oversized cluster
    for (cluster_id in need_recluster) {
        # Get indices for this cluster
        subset_indices <- which(cluster_names == cluster_id)
        subset_data <- data[subset_indices, , drop = FALSE]
        
        # Recursively cluster this subset
        sub_clusters <- recursiveDynamicCut(
            data = subset_data,
            dist_matrix = NULL,  # Recalculate for subset
            hclust_obj = NULL,   # Recalculate for subset
            minCl = minCl,
            maxCl = maxCl,
            cluster_prefix = cluster_id,  # Use current cluster name as prefix
            max_depth = max_depth,
            current_depth = current_depth + 1,
            pad_width = pad_width
        )
        
        # Update cluster assignments
        cluster_names[subset_indices] <- sub_clusters
    }
    
    return(cluster_names)
}

# # Example usage:
# # Set parameters
# minCl <- floor(nrow(allStudies) * 0.01)  # minimum cluster size (1%)
# maxCl <- floor(nrow(allStudies) * 0.1)   # maximum cluster size (10%)
# 
# # Calculate initial distance and clustering
# res.dist <- factoextra::get_dist(allStudies, method = "spearman")
# res.hclust <- hclust(res.dist, method = "ward.D")
# 
# # Run recursive clustering with zero-padding
# final_clusters <- recursiveDynamicCut(
#     data = allStudies,
#     dist_matrix = res.dist,
#     hclust_obj = res.hclust,
#     minCl = minCl,
#     maxCl = maxCl,
#     max_depth = 10,
#     pad_width = 3  # Adjust based on expected max clusters (3 = up to 999)
# )
# 
# # Check results
# head(final_clusters)  # Named vector with hierarchical cluster IDs
# table(final_clusters)
# max(table(final_clusters)) <= maxCl  # Should be TRUE
# 
# # Now sorting works correctly:
# sort(unique(final_clusters))
# # Output will be: "001", "002", "002.001", "002.002", "002.010", etc.
# # Instead of: "1", "1.1", "1.10", "1.11", "1.2", etc.
# 
# # Example usage:
# res_adaptive <- recursiveDynamicCut(data = all,
#                                     dist_matrix = res.dist,
#                                     hclust_obj = res.hclust,
#                                     minCl = 5 ,
#                                     maxCl = 53,
#                                     cluster_prefix = "",
#                                     max_depth = 10,
#                                     current_depth = 0,
#                                     pad_width = 3)
