library(doParallel)
library(parallel)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
library(BiocParallel)
library(data.table)
library(foreach)
#' @importFrom dplyr mutate
#' @importFrom magrittr %>%
#' @importFrom parallel stopCluster
#' @importFrom utils globalVariables

utils::globalVariables(c("group_vector", "cluster_id", "p_val"))
# Define a function to perform differential expression analysis for all clusters
cytoKSCE_clusters_Fpsrf <- function(sce,
                                    group_factor, covars, residuals,
                                    initial_permutations, thresholds, perm_increments,
                                    residual_tolerance, random_seed,
  name_assays_expression = "logcounts",
  name_cluster = "cluster_id",
  name_sample = "sample_id", n_cores, FDR_method) {

    # Input validation
    stopifnot(
      (is(sce, "SummarizedExperiment") | is(sce, "SingleCellExperiment")),
      is.character(name_assays_expression), length(name_assays_expression) == 1L,
      is.character(name_cluster), length(name_cluster) == 1L,
      is.character(name_sample), length(name_sample) == 1L
    )



    # Count matrix
    sel = which(names(assays(sce)) == name_assays_expression)
    if (length(sel) == 0) {
      message("'", name_assays_expression, "' not found in names(assays(x))")
      return(NULL)
    }
    if (length(sel) > 1) {
      message("'", name_assays_expression, "' found multiple times in names(assays(x))")
      return(NULL)
    }
    counts = assays(sce)[[sel]]

    if (any(is.na(counts)) | any(is.null(counts))) {
      message("'assays(sce)$", name_assays_expression, "' contains NA or NULL values")
      return(NULL)
    }

    # Remove rows with 0 counts
    #if (nrow(counts) > 1.5) { # Only if more than 1 gene (row): otherwise matrix is transformed into a vector
     # counts = counts[rowSums(counts > 0) > 0, ]
   # }

    # Check if counts are sparse matrix: if not, turn counts into Sparse object
    if (!is(counts, "dgCMatrix")) {
      counts = Matrix::Matrix(data = counts,
                      sparse = TRUE)
    }

  # cluster ids:
  sel = which(names(colData(sce)) == name_cluster)
  if( length(sel) == 0 ){
    message("'", name_cluster,"' not found in names(colData(sce))")
    return(NULL)
  }
  if( length(sel) > 1 ){
    message("'", name_cluster,"' found multiple times in names(colData(sce))")
    return(NULL)
  }
  cluster_ids = unique(factor(colData(sce)[[sel]]))
  # Define the loop
  # Create a parallel backend with the desired number of workers
  #cl <- parallel::makeCluster(parallel::detectCores())
  #cl<- n_cores
  cl <- parallel::makeCluster(n_cores)
  # Register the parallel backend
  doParallel::registerDoParallel(cl)
  #res21 <- foreach(cluster_id = unique(colData(sce)$cell), .combine = rbind) %do% {
  res21 <- foreach(cluster_id = cluster_ids, .combine = rbind) %do% {
    # Subset the SingleCellExperiment object based on the cluster ID value
    #cluster_id<- "CD4 T cells"
    sce_subset <- sce[, cluster_ids == cluster_id]
    ##########.  sce subset ####################################

    # Count matrix
    sel = which(names(assays(sce_subset)) == name_assays_expression)
    if (length(sel) == 0) {
      message("'", name_assays_expression, "' not found in names(assays(x))")
      return(NULL)
    }
    if (length(sel) > 1) {
      message("'", name_assays_expression, "' found multiple times in names(assays(x))")
      return(NULL)
    }
    counts = assays(sce_subset)[[sel]]

    if (any(is.na(counts)) | any(is.null(counts))) {
      message("'assays(sce_subset)$", name_assays_expression, "' contains NA or NULL values")
      return(NULL)
    }

    # Remove rows with 0 counts
    #if (nrow(counts) > 1.5) { # Only if more than 1 gene (row): otherwise matrix is transformed into a vector
    # counts = counts[rowSums(counts > 0) > 0, ]
    # }

    # Check if counts are sparse matrix: if not, turn counts into Sparse object
    if (!is(counts, "dgCMatrix")) {
      counts = Matrix::Matrix(data = counts,
                              sparse = TRUE)
    }


    ################################. sce subset ####################################################################








    # Perform differential expression analysis using CytoKsceProc
    #sample_id <- sce_subset$sample_id
    # sample ids:
    sel = which(names(colData(sce_subset)) == name_sample)
    if( length(sel) == 0 ){
      message("'", name_sample,"' not found in names(colData(sce_subset))")
      return(NULL)
    }
    if( length(sel) > 1 ){
      message("'", name_sample,"' found multiple times in names(colData(sce_subset))")
      return(NULL)
    }
    sample_id = factor(colData(sce_subset)[[sel]])
    #res2 <- CytoKsceProc(counts, group_factor,
        #              sample_id,bandwidth, covars,
        #              lowerRho, upperRho, gridRho, n_cores)



    res2 <- CytoKsceProc_Fpsrf(counts, group_factor, sample_id, covars, residuals=TRUE,
                               initial_permutations=100, thresholds = c(0.1, 0.01, 0.001),
                               perm_increments = c(500, 2000, 10000),
                               random_seed = 999, residual_tolerance = 0.5, n_cores, FDR_method)
    res2$gene <- rownames(sce_subset)
    res2$cluster_id <- rep(cluster_id, length(res2$gene))
    res2
  }

  # Reorder columns to have 'cluster_id' as the first column
  res21 <- res21[, c("gene", "cluster_id", setdiff(colnames(res21), c("cluster_id", "gene")))]
library(dplyr)
  #res21 <- res21 %>%
    #group_by(cluster_id) %>%
    #mutate(AdjPvalue_clusters = p.adjust(vec_pValue, method = "BH"))

  res21 <- res21 %>%
    #group_by(cluster_id) %>%
    mutate(p_adj.glb = adjust_p_values(p_val, FDR_method))


  # Stop and close the parallel cluster
  stopCluster(cl)
  registerDoSEQ()
  #print(cluster_id)
  # Return the differential expression results as a data frame
  res21
}



