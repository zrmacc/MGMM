# Purpose: Master estimation function for MGMM.
# Updated: 20/07/19

#' Estimate Multivariate Normal Mixture
#'
#' Given an \eqn{n \times d} matrix of random vectors, estimates the parameters
#' of a Gaussian Mixture Model (GMM). Accommodates arbitrary patterns of missingness
#' at random (MAR) in the input vectors.
#'
#' Initial values for the cluster means, covariances, and proportions are
#' specified using \code{M0}, \code{S0}, and \code{pi0}, respectively. If the
#' data contains complete observations, i.e. observations with no missing
#' elements, then \code{fit.GMM} will attempt to initialize these parameters
#' internally using K-means. If the data contains no complete observations, then
#' initial values are required for \code{M0}, \code{S0}, and \code{pi0}.
#'
#' @param data Numeric data matrix.
#' @param k Number of mixture components. Defaults to 1.
#' @param init_means Optional list of initial mean vectors.
#' @param fix_means Fix the means to their starting value? Must provide initial
#'   values.
#' @param init_covs Optional list of initial covariance matrices.
#' @param init_props Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return For a single component model \eqn{k=1}, a list is returned,
#'   containing the estimated mean, covariance, and final EM objective. For a
#'   multi-component model \eqn{k>1}, an object of class \code{mix}, containing
#'   the estimated means, covariances, cluster proportions, cluster
#'   responsibilities, and observation assignments.
#'
#' @export
#' @seealso See \code{\link{rGMM}} for data generation, and \code{\link{ChooseK}} for selecting
#' the number of clusters.
#'
#' @examples
#' \donttest{
#' # Single component without missingness
#' # Bivariate normal observations
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' data <- rGMM(n = 1e3, d = 2, k = 1, means = c(2, 2), covs = sigma)
#' fit <- fit.GMM(data, k = 1)
#' 
#' # Single component with missingness
#' # Trivariate normal observations
#' mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
#' sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = sigma)
#' fit <- fit.GMM(data, k = 2)
#' 
#' # Two components without missingness
#' # Trivariate normal observations
#' mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
#' sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = sigma)
#' fit <- fit.GMM(data, k = 2)
#' 
#' # Four components with missingness
#' # Bivariate normal observations
#' # Note: Fitting is slow.
#' mean_list <- list(c(2, 2), c(2, -2), c(-2, 2), c(-2, -2))
#' sigma <- 0.5 * diag(2)
#' data <- rGMM(
#' n = 1000, 
#' d = 2, 
#' k = 4, 
#' pi = c(0.35, 0.15, 0.15, 0.35), 
#' m = 0.1, 
#' means = mean_list, 
#' covs = sigma)
#' fit <- fit.GMM(data, k = 4)
#' }

fit.GMM <- function(
  data, 
  k = 1, 
  init_means = NULL, 
  fix_means = FALSE, 
  init_covs = NULL, 
  init_props = NULL, 
  maxit = 100, 
  eps = 1e-6, 
  report = TRUE
) {
  
  # Check data.
  if (!is.matrix(data)) {
    stop("Data are expected as a numeric matrix.")
  }
  d <- ncol(data)
  
  # Check initial values.
  ## Mean vectors:
  if (!is.null(init_means)) {
    
    ## Object class.
    if (!is.list(init_means)) {
      stop("If provided, init_means is expected as a list of numeric vectors, one for each mixture component.")
    }
    
    ## Initial mean for each component.
    if (length(init_means) != k) {
      stop("If initial means are provided, one is required for each mixture component.")
    }
    
    ## Dimensional consistency
    mean_dims <- unique(unlist(lapply(init_means, length)))
    if ((length(mean_dims) > 1) | (mean_dims != d)) {
      stop("Each vector in init_means must have length of ncol(data).")
    }
    
  }
  
  if (fix_means & is.null(init_means)) {
    stop("If means are fixed, then initial values are required.")
  }
  
  ## Covariance matrices:
  if (!is.null(init_covs)) {
    
    ## Object class
    if (!is.list(init_covs)) {
      stop("If provided, init_covs is expected as a list of numeric matrices, one for each mixture component.")
    }
    
    ## Initial covariance for each component
    if (length(init_covs) != k) {
      stop("If initial covariances are provided, one is required for each mixture component.")
    }
    
    ## Dimensional consistency
    cov_dims <- unique(unlist(lapply(init_covs, dim)))
    if ((length(cov_dims) > 1) | (cov_dims != d)) {
      stop("Each matrix in init_means must have dimensions of ncol(data) by ncol(data).")
    }
  }
  
  ## Cluster proportions
  if (!is.null(init_props)) {
    # Object type
    if (!is.numeric(init_props)) {
      stop("If provided, init_props is expected as a numeric vector.")
    }
    # Initial proportion for each component
    if (length(init_props) != k) {
      stop("If initial proportions are provided, one is required for each mixture component.")
    }
  } 

  # Estimation. 
  if (k == 1) {
    
    ## Case 1: Single mixture component.
    out <- fit.mvn(
      data,
      init_means[[1]],
      fix_means,
      init_covs[[1]],
      maxit, 
      eps,
      report
    )
    
  } else {
    
    ## Case 2: Multiple mixture components.
    out <- fit.mix(
      data,
      k,
      init_means,
      fix_means,
      init_covs,
      init_props,
      maxit,
      eps,
      report
    )
    
  }
  
  # Output.
  return(out)
}
