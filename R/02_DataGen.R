# Purpose: Generate data from Gaussian Mixture Models.
# Updated: 2021-07-24

#' Generate Data from Gaussian Mixture Models
#'
#' Generates an \eqn{n\times d} matrix of multivariate normal random vectors
#' with observations (examples) as rows. If \eqn{k=1}, all observations belong to the same
#' cluster. If \eqn{k>1} the observations are generated via two-step procedure.
#' First, the cluster membership is drawn from a multinomial distribution, with
#' mixture proportions specified by \code{pi}. Conditional on cluster
#' membership, the observation is drawn from a multivariate normal distribution,
#' with cluster-specific mean and covariance. The cluster means are provided
#' using \code{means}, and the cluster covariance matrices are provided using
#' \code{covs}. If \eqn{miss>0}, missingness is introduced, completely at random, by
#' setting that proportion of elements in the data matrix to \code{NA}.
#'
#' @param n Observations (rows).
#' @param d Observation dimension (columns).
#' @param k Number of mixture components. Defaults to 1.
#' @param pi Mixture proportions. If omitted, components are assumed
#'   equiprobable.
#' @param miss Proportion of elements missing, \eqn{miss\in[0,1)}.
#' @param means Either a prototype mean vector, or a list of mean vectors. Defaults
#'   to the zero vector.
#' @param covs Either a prototype covariance matrix, or a list of covariance
#'   matrices. Defaults to the identity matrix.
#' @return Numeric matrix with observations as rows. Row numbers specify the
#'   true cluster assignments.
#'
#' @export
#' @seealso For estimation, see \code{\link{FitGMM}}.
#'
#' @examples
#' set.seed(100)
#' # Single component without missingness.
#' # Bivariate normal observations.
#' cov <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' data <- rGMM(n = 1e3, d = 2, k = 1, means = c(2, 2), covs = cov)
#' 
#' # Single component with missingness.
#' # Trivariate normal observations.
#' mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
#' cov <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = cov)
#' 
#' # Two components without missingness.
#' # Trivariate normal observations.
#' mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
#' cov <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = cov)
#' 
#' # Four components with missingness.
#' # Bivariate normal observations.
#' mean_list <- list(c(2, 2), c(2, -2), c(-2, 2), c(-2, -2))
#' cov <- 0.5 * diag(2)
#' data <- rGMM(
#' n = 1000, 
#' d = 2, 
#' k = 4, 
#' pi = c(0.35, 0.15, 0.15, 0.35), 
#' miss = 0.1, 
#' means = mean_list, 
#' covs = cov)

rGMM <- function(
  n, 
  d = 2, 
  k = 1, 
  pi = NULL, 
  miss = 0, 
  means = NULL, 
  covs = NULL
) {
  
  # Input checks.
  ## Mixture proportions.
  if (is.null(pi)) {
    pi <- rep(1, k) / k
  } else {
    pi <- pi / sum(pi)
  }
  
  ## Means.
  if (is.null(means)) {
    means <- rep(0, d)
  }
  
  if (!is.list(means)) {
    means <- list(means)
    means <- rep(means, k)
  }
  
  ## Covariance matrices.
  if (is.null(covs)) {
    covs <- diag(d)
  }
  
  if (!is.list(covs)) {
    covs <- list(covs)
    covs <- rep(covs, k)
  }
  
  ## Check dimensional consistency.
  if (k != length(pi)) {
    stop("pi must have length k.")
  }
  
  mean_dims <- unique(unlist(lapply(means, length)))
  if (mean_dims != d) {
    stop("Each vector in means must have length d.")
  }
  
  cov_dims <- unique(unlist(lapply(covs, dim)))
  if (cov_dims != d) {
    stop("Each matrix in covs must have dimension d x d.")
  }
  
  ## Check input validity.
  if (miss < 0 | miss >= 1) {
    stop("Missingness must fall in [0,1).")
  }

  if ((length(pi) > 1) & (min(pi) <= 0 | max(pi) >= 1)) {
    stop("Elements of pi must reside in (0,1).")
  }

  # Data generation.
  if (k == 1) {
    
    # Case 1: Single mixture component.
    mu <- means[[1]]
    sigma <- covs[[1]]
    data <- mvnfast::rmvn(n = n, mu = mu, sigma = sigma)
    row.names(data) <- rep(1, n)
    
  } else {
    
    # Case 2: Multiple mixture components.
    comp <- sample(x = seq_len(k), size = n, replace = TRUE, prob = pi)
    
    # Loop over mixture components. 
    data <- lapply(seq_len(k), function(i){
      ni <- sum(comp == i)
      
      # Output only if component is non-empty
      if (ni > 0) {
        mu <- means[[i]]
        sigma <- covs[[i]]
        data <- mvnfast::rmvn(n = ni, mu = mu, sigma = sigma)
        row.names(data) <- rep(i, ni)
        return(data)
      }
    })

    data <- do.call(rbind, data)
    
    # Permute row order
    data <- data[sample(nrow(data), replace = FALSE), , drop = FALSE]
  }

  # Introduce missingness.
  if (miss > 0) {
    ele <- length(data)
    n_miss <- floor(ele * miss)
    draw <- sort(sample(x = ele, size = n_miss, replace = FALSE))
    data[draw] <- NA
  }

  # Output.
  colnames(data) <- paste0("y", seq_len(d))
  return(data)
}
