#' @useDynLib cytoKernel
#' @importFrom Rcpp evalCpp
#' @exportPattern '^[[:alpha:]]+'
NULL

# Load necessary libraries
library(Matrix)

.calculate_density <- function(data1, data2, points = 1024) {
  cur_range <- range(c(data1, data2), na.rm = TRUE, finite = TRUE)
  density1 <- density(data1, from = cur_range[1], to = cur_range[2], n = points, na.rm = TRUE)
  density2 <- density(data2, from = cur_range[1], to = cur_range[2], n = points, na.rm = TRUE)
  res <- list(density1 = density1$y, density2 = density2$y)
  return(res)
}

.calculate_jsd <- function(dist1, dist2, min_threshold = 1e-10) {
  dist1[which(dist1 < min_threshold)] <- min_threshold
  dist2[which(dist2 < min_threshold)] <- min_threshold
  mean_dist <- (dist1 + dist2) / 2
  - (sum(dist1 * (log(mean_dist / dist1)), na.rm = TRUE) + sum(dist2 * (log(mean_dist / dist2)), na.rm = TRUE)) / 2
}

.JSD_element <- function(data1, data2, points = 1024) {
  res <- .calculate_density(data1, data2, points = points)
  dist1 <- res$density1 + 1 / points
  dist2 <- res$density2 + 1 / points
  dist1 <- dist1 / sum(dist1)
  dist2 <- dist2 / sum(dist2)
  dist_val <- .calculate_jsd(dist1, dist2)
  return(dist_val)
}


# Main .KSTsce function with adapted PERMANOVA method
.KSTsce_Fpsrf <- function(featureVec, group_factor, sample_id, covars, residuals,
                     initial_permutations, thresholds, perm_increments,
                     residual_tolerance, random_seed) {

  x <- featureVec
  n <- length(group_factor)
  unique_sample_id<- unique(sample_id)
  JSD_matrix <- matrix(0, n, n) # Replace with your distance calculation

  for (i in 1:n) {
    for (j in i:n) {
      # Extract indices for the current pair of samples

      # Extract indices for the current pair of samples
      indices_i <- .getMatchIndices2(sample_id, unique_sample_id[i])
      indices_j <- .getMatchIndices2(sample_id, unique_sample_id[j])

      # Calculate the distance between the two samples
      distance <- .JSD_element(x[indices_i], x[indices_j], points = 1024)


      # Calculate the distance between the two samples
      #distance <- JSD_element(x[indices_i], x[indices_j], points = 1024)

      # Assign the calculated distance to the matrix
      JSD_matrix[i, j] <- sqrt(distance)
      JSD_matrix[j, i] <- JSD_matrix[i, j]  # Ensure the matrix is symmetric
      #k[i,i] = 1
    }
  }
  JSD_matrix[is.na(JSD_matrix)] <- 0  # Replace NA with 0
  JSD_matrix[is.infinite(JSD_matrix)] <- max(JSD_matrix[!is.infinite(JSD_matrix)], na.rm = TRUE)  # Replace Inf with max finite value

  # Set the diagonal elements to zero
  diag(JSD_matrix) <- 0
#print(JSD_matrix)
  # Execute PERMANOVA with or without residuals
  pValue <- .CytoKSCE_psrf_pVal3(JSD_matrix, group_factor, covars,
                         residuals,
                         initial_permutations, thresholds,
                         perm_increments,
                         random_seed,
                         response_y = "binary", residual_tolerance)

  return(pValue)
}
