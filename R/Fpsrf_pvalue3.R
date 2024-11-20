# Load necessary libraries
library(Matrix)
library(stats)


# Define the custom PERMANOVA function designed to replicate ideas::permanova
.CytoKSCE_psrf_pVal3 <- function(JSD_matrix, group_factor, covars = NULL,
                        residuals = FALSE,
                        initial_permutations, thresholds = c(0.1, 0.01, 0.001),
                        perm_increments = c(500, 2000, 10000),
                        random_seed = 999,
                        response_y = "binary", residual_tolerance = 0.5) {

  # Set random seed for reproducibility
  set.seed(random_seed)

  # Handle group_factor based on response_y
  if (response_y == "binary") {
    if (is.character(group_factor) || is.factor(group_factor)) {
      group_factor <- as.numeric(as.factor(group_factor)) - 1
    }
  }

  if (length(group_factor)>9) {
    #group_vector<- group_factor
    # Handle additional predictors as covariates
    if (!is.null(covars)) {
      covariate_matrix <- model.matrix(~ ., data = as.data.frame(covars))
    } else {
      covariate_matrix <- NULL
    }

    # Permutation setup
    x_perm <- replicate(initial_permutations, sample(group_factor))

    # Calculate observed F-statistic
    if (is.null(covariate_matrix)) {
      if (response_y == "binary") {
        F_observed <- .F_manova_comp(JSD_matrix, group_factor)
        F_permutations <- apply(x_perm, 2, function(perm) .F_manova_comp(JSD_matrix, perm))
      } else {
        F_observed <- .F_permanova_comp(JSD_matrix, residualsZ = group_factor)
        F_permutations <- apply(x_perm, 2, function(perm) .F_permanova_comp(JSD_matrix, residualsZ = perm))
      }
    } else {
      if (residuals) {

        if (response_y == "binary") {
          # Step 1: Fit logistic model to obtain residuals
          model <- glm(group_factor ~ -1 + covariate_matrix, family = binomial(link = "logit"))
          fitted_values <- fitted(model)
          residuals_x <- group_factor - fitted_values  # Compute residuals manually
          sd_e <- sd(residuals_x)  # Expected standard deviation of residuals

          # Initialize matrix for residualized permutations
          residual_permutations <- matrix(NA, nrow = length(group_factor), ncol = initial_permutations)
          ip <- 0  # index for usable permutations
          id <- 0  # working index

          # Step 2: Generate permutations with residuals that match observed sd within residual_tolerance tolerance
          while (ip < initial_permutations) {
            id <- id + 1
            perm_x <- rbinom(length(fitted_values), 1, prob = fitted_values)

            # Fit logistic model to permuted group_factors
            perm_model <- tryCatch(glm(perm_x ~ -1 + covariate_matrix, family = binomial(link = "logit")),
                                   warning = function(w) { NULL },
                                   error = function(e) { NULL },
                                   finally = {})

            if (!is.null(perm_model)) {
              perm_resid <- perm_x - fitted(perm_model)
              sd_resid_i <- sd(perm_resid)

              # Check if permuted residual sd is within residual_tolerance tolerance of observed residual sd
              if (sd_resid_i > (1 - residual_tolerance) * sd_e && sd_resid_i < (1 + residual_tolerance) * sd_e) {
                residual_permutations[, ip + 1] <- perm_resid
                ip <- ip + 1
              }
            }
          }

          F_observed <- .F_permanova_comp(JSD_matrix, residuals_x, covariate_matrix)
          F_permutations <- apply(residual_permutations, 2, function(resid) .F_permanova_comp(JSD_matrix, resid, covariate_matrix))

        } else {
          # Step 1: Fit linear model to obtain residuals for continuous variable
          model <- lm(group_factor ~ -1 + covariate_matrix)
          fitted_values <- fitted(model)
          residuals_x <- group_factor - fitted_values

          # Generate continuous residual permutations
          residual_permutations <- matrix(NA, nrow = length(group_factor), ncol = initial_permutations)
          for (ip in 1:initial_permutations) {
            perm_x <- sample(residuals_x, size = length(residuals_x)) + fitted_values
            perm_model <- lm(perm_x ~ -1 + covariate_matrix)
            residual_permutations[, ip] <- perm_x - fitted(perm_model)
          }

          F_observed <- .F_permanova_comp(JSD_matrix, residuals_x, covariate_matrix)
          F_permutations <- apply(residual_permutations, 2, function(resid) .F_permanova_comp(JSD_matrix, resid, covariate_matrix))
        }
      }
      else {
        # No residualization, include covariate adjustment
        F_observed <- .F_permanova_comp(JSD_matrix, group_factor, covariate_matrix)
        F_permutations <- apply(x_perm, 2, function(perm) .F_permanova_comp(JSD_matrix, perm, covariate_matrix))
      }
    }

    # Calculate p-value
    p_value <- mean(F_permutations >= F_observed)
    # -----------------------------------------------------------------
    # Adaptive Permutation Test with Dynamic Permutations
    # -----------------------------------------------------------------
    #JSD_matrix<- JSD_matrix
    p_value <- (sum(F_permutations >= F_observed) + 1) / (initial_permutations + 1)
    permutations <- initial_permutations

    # Increase permutations based on thresholds
    for (i in seq_along(thresholds)) {
      if (p_value <= thresholds[i] && permutations < perm_increments[i]) {
        additional_permutations <- perm_increments[i] - permutations
        extra_perms <- replicate(additional_permutations, {
          perm_group <- sample(group_vector)
          .F_permanova_comp(JSD_matrix, perm_group, covariate_matrix)
        })

        F_permutations <- c(F_permutations, extra_perms)
        permutations <- perm_increments[i]
        p_value <- (sum(F_permutations >= F_observed) + 1) / (permutations + 1)
      }
      # Stop if p-value exceeds threshold after increasing permutations
      if (p_value > thresholds[i]) break
    }

    pValue<- p_value
    return(pValue)
  } else {

    # Generate covariate matrix if provided
    if (!is.null(covars)) {
      covariate_matrix <- model.matrix(~ ., data = as.data.frame(covars))
    } else {
      #covariate_matrix <- NULL
      covariate_matrix<- rep(1, length(group_factor))
    }

    # Set random seed for reproducibility
    set.seed(random_seed)

    # -----------------------------------------------------------------
    # Calculate observed F-statistic
    # -----------------------------------------------------------------

    if (is.null(covariate_matrix)) {
      F_observed <- .F_manova_comp(JSD_matrix, group_factor)
    } else {
      if (residuals) {
        # Residualize group variable with additional predictors
        model <- suppressWarnings(glm(group_factor ~ -1 + covariate_matrix, family = binomial(link = "logit")))
        fitted_values <- fitted(model)
        residuals_x <- group_factor - fitted_values

        # Calculate residualized F-statistic for observed data
        F_observed <- .F_permanova_comp(JSD_matrix, residuals_x, covariate_matrix)
      } else {
        # No residualization, directly compute F-statistic
        F_observed <- .F_permanova_comp(JSD_matrix, group_factor, covariate_matrix)
      }
    }

    # -----------------------------------------------------------------
    # Exact Permutation Test for Small Samples
    # -----------------------------------------------------------------

    n <- length(group_factor)
    F_permutations <- numeric()

    # -----------------------------------------------------------------
    # Approximate Permutation Test with Adaptive Permutations
    # -----------------------------------------------------------------

    F_permutations <- replicate(initial_permutations, {
      perm_group <- sample(group_factor)
      .F_permanova_comp(JSD_matrix, perm_group, covariate_matrix)
    })

    p_value <- (sum(F_permutations >= F_observed) + 1) / (initial_permutations + 1)
    permutations <- initial_permutations

    # Increase permutations based on thresholds if needed
    for (i in seq_along(thresholds)) {
      if (p_value <= thresholds[i] && permutations < perm_increments[i]) {
        additional_permutations <- perm_increments[i] - permutations
        extra_perms <- replicate(additional_permutations, {
          perm_group <- sample(group_factor)
          .F_permanova_comp(JSD_matrix, perm_group, covariate_matrix)
        })

        F_permutations <- c(F_permutations, extra_perms)
        permutations <- perm_increments[i]
        p_value <- (sum(F_permutations >= F_observed) + 1) / (permutations + 1)
      }
      # Stop if p-value exceeds threshold after increasing permutations
      if (p_value > thresholds[i]) break
    }

    pValue<- p_value
    return(pValue)
  }

}
# -----------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------



# Function to center and project a distance matrix for PERMANOVA
.G_comp <- function(D, epsilon = 1e-6) {
  n <- nrow(D)
  centerD <- diag(n) - 1 / n
  G <- -0.5 * centerD %*% (D * D) %*% centerD

  # Perform eigen decomposition
  eG <- eigen(G, symmetric = TRUE)

  # Set negative eigenvalues to zero and add a small positive value to ensure positive definiteness
  adjusted_values <- pmax(eG$values, 0) + epsilon

  # Reconstruct the positive definite matrix
  G <- eG$vectors %*% diag(adjusted_values) %*% t(eG$vectors)
  return(G)
}

.calTrace_sqrt2 <- function(G, H) {
  # Compute eigen decomposition of G
  eigen_G <- eigen(G, symmetric = TRUE)

  # Construct square root of G using only real parts of eigenvalues (if they are non-negative)
  sqrt_eigenvalues <- sqrt(pmax(Re(eigen_G$values), 0))
  real_sqrt_G <- eigen_G$vectors %*% diag(sqrt_eigenvalues) %*% t(eigen_G$vectors)

  # Calculate trace as sum of element-wise product without forming full matrix product
  trace_value <- sum(H * (real_sqrt_G %*% H))

  return(trace_value)
}


# Function to calculate F-statistic for PERMANOVA
.F_permanova_comp <- function(JSD_matrix, residualsZ, Z = NULL, regularization = 1e-6) {
  n <- nrow(JSD_matrix)
  G <- .G_comp(JSD_matrix)

  if (is.null(Z)) {
    XZ <- as.matrix(residualsZ)
  } else {
    XZ <- cbind(residualsZ, Z)
  }

  H <- XZ %*% solve(t(XZ) %*% XZ + diag(regularization, ncol(XZ))) %*% t(XZ)
  IH <- diag(n) - H

  t1 <- .calTrace_sqrt2(G, H)
  t2 <- .calTrace_sqrt2(G, IH)

  F_pseudo_sqrt <- t1 / t2
  return(F_pseudo_sqrt)
}

# Calculate F-statistic for MANOVA
.F_manova_comp <- function(JSD_matrix, group_factor) {
  unique_group_factor <- unique(group_factor)
  a <- length(unique_group_factor)
  n <- length(group_factor)

  d2 <- JSD_matrix^2
  SST <- sum(d2) / n

  denomVector <- matrix(0, n, n)
  for (i in 1:a) {
    group_index <- which(group_factor == unique_group_factor[i])
    denomVector[group_index, group_index] <- 1 / length(group_index)
  }

  SSW <- sum(d2 * as.vector(denomVector))
  F_pseudo_sqrt <- ((SST - SSW) * (n - a)) / (SSW * (a - 1))
  return(F_pseudo_sqrt)
}




