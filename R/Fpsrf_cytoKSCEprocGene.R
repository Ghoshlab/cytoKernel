library(doParallel)
library(SingleCellExperiment)
library(BiocParallel)
library(data.table)

adjust_p_values <- function(p_values, FDR_method) {
  # Validate method input
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni",
                     "BH", "BY", "fdr", "none")
  if (!FDR_method %in% valid_methods) {
    stop(paste("Invalid FDR_method. Choose from:", paste(valid_methods, collapse = ", ")))
  }

  # Adjust p-values using the specified method.arg
  adjusted_p_values <- stats::p.adjust(p_values, method = FDR_method)

  # Return the adjusted p-values
  return(adjusted_p_values)
}

CytoKsceProc_Fpsrf <- function(object, group_factor, sample_id, covars, residuals,
                         initial_permutations, thresholds, perm_increments,
                         residual_tolerance, random_seed, n_cores,FDR_method) {
  BPPARAM <- BiocParallel::MulticoreParam(workers=n_cores)

  pValue_features <- BiocParallel::bplapply(seq_len(dim(object)[1]),
                                            function(h){.KSTsce_Fpsrf(object[h,], group_factor, sample_id, covars, residuals,
                                                                      initial_permutations, thresholds, perm_increments,
                                                                      residual_tolerance, random_seed)},
                                            BPPARAM= BPPARAM)
  vec_pValue<- unlist(pValue_features)
  AdjPvalue_features<-  adjust_p_values(vec_pValue, FDR_method)
  return(data.frame(p_val=vec_pValue,
                    p_adj.loc = AdjPvalue_features))
}
