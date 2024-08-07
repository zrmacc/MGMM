# Purpose: Master estimation function for MGMM.
# Updated: 2021-07-24

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
#' @param lambda Optional ridge term added to covariance matrix to ensure 
#'   positive definiteness.
#' @param init_props Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return 
#' \itemize{
#'   \item For a single component, an object of class \code{mvn}, containing
#'   the estimated mean and covariance, the final objective function, and the
#'   imputed data. 
#'   \item For a multicomponent model \eqn{k>1}, an object of class \code{mix}, 
#'   containing the estimated means, covariances, cluster proportions, cluster
#'   responsibilities, and observation assignments.
#' }
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
#' fit <- FitGMM(data, k = 1)
#' 
#' # Single component with missingness
#' # Trivariate normal observations
#' mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
#' sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = sigma)
#' fit <- FitGMM(data, k = 2)
#' 
#' # Two components without missingness
#' # Trivariate normal observations
#' mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
#' sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = sigma)
#' fit <- FitGMM(data, k = 2)
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
#' fit <- FitGMM(data, k = 4)
#' }

FitGMM <- function(
  data, 
  k = 1, 
  init_means = NULL, 
  fix_means = FALSE, 
  init_covs = NULL, 
  lambda = 0,
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
    
    ## Dimensional consistency.
    mean_lengths <- unique(unlist(lapply(init_means, length)))
    if ((length(mean_lengths) > 1) | (mean_lengths != d)) {
      stop("Each vector in init_means must have length of ncol(data).")
    }
    
  }
  
  ## Fixed means.
  if (fix_means & is.null(init_means)) {
    stop("If means are fixed, then initial values are required.")
  }
  
  ## Covariance matrices:
  if (!is.null(init_covs)) {
    
    ## Object class.
    if (!is.list(init_covs)) {
      stop("If provided, init_covs is expected as a list of numeric matrices, one for each mixture component.")
    }
    
    ## Initial covariance for each component.
    if (length(init_covs) != k) {
      stop("If initial covariances are provided, one is required for each mixture component.")
    }
    
    ## Dimensional consistency.
    cov_dims <- unique(unlist(lapply(init_covs, dim)))
    if ((length(cov_dims) > 1) | (cov_dims != d)) {
      stop("Each matrix in init_means must have dimensions of ncol(data) by ncol(data).")
    }
  }
  
  ## Cluster proportions
  if (!is.null(init_props)) {
    
    ## Object type
    if (!is.numeric(init_props)) {
      stop("If provided, init_props is expected as a numeric vector.")
    }
    
    ## Initial proportion for each component.
    if (length(init_props) != k) {
      stop("If initial proportions are provided, one is required for each mixture component.")
    }
  } 

  # Estimation. 
  ## Case 1: Single mixture component.
  if (k == 1) {
    out <- FitMVN(
      data = data,
      init_mean = init_means[[1]],
      fix_mean = fix_means,
      init_cov = init_covs[[1]],
      lambda = lambda,
      maxit = maxit, 
      eps = eps,
      report = report
    )

  ## Case 2: Multiple mixture components. 
  } else {
    out <- FitMix(
      data = data,
      k = k,
      init_means = init_means,
      fix_means = fix_means,
      init_covs = init_covs,
      lambda = lambda,
      init_props = init_props,
      maxit = maxit,
      eps = eps,
      report = report
    )

  }
  
  # Output.
  return(out)
}
