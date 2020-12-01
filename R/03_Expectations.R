# Purpose: Calculate E-step expectations.
# Updated: 20/07/16. 

#-------------------------------------------------------------------------------
# Responsibilities
#-------------------------------------------------------------------------------

#' Evaluate the Density of an Incomplete Observations
#' 
#' @param y Vector with missing elements.
#' @param means List of mean vectors.
#' @param covs List of covariances matrices.
#' @param pi Vector of cluster proportions.
#' @return Numeric density of observed elements.

Responsibility.eval_dens_incomp <- function(
  y, 
  means,
  covs,
  pi
) {
  
  # Mixture components. 
  k <- length(pi)
  
  # Observed elements
  is_obs <- !is.na(y)
  obs_ele = y[is_obs]
  
  # Loop over components of mixture:
  aux <- function(j) {
    obs_mean <- means[[j]][is_obs]
    obs_cov <- covs[[j]][is_obs, is_obs]
    return(dmvn(X = obs_ele, mu = obs_mean, sigma = obs_cov) * pi[j])
  }
  
  # Density evaluation of observed elements
  obs_dens <- unlist(lapply(seq(1:k), aux))
  return(obs_dens)
}


#' Responsibilities
#'
#' Calculates the posterior probability of cluster membership given the observed
#' data.
#'
#' @param split_data Data partitioned by missingness.
#' @param means List of mean vectors.
#' @param covs List of covariances matrices.
#' @param pi Vector of cluster proportions.
#' @return List containing:
#' \itemize{
#'   \item k Number of mixture components.
#'   \item Density evaluations `dens_eval0` and responsibilities `gamma0` 
#'     for complete cases.
#'   \item Density evaluations `dens_eval1` and responsibilities `gamma1` 
#'     for incomplete cases. 
#' }
#'
#' @importFrom plyr aaply

Responsibility <- function(
  split_data,
  means, 
  covs, 
  pi
) {
  
  # Dimensions
  n0 <- split_data$n0
  n1 <- split_data$n1
  
  # Clusters
  k <- length(pi)

  # Output structure
  out <- list()
  out$k <- k

  # Multivariate normal density evaluation for complete observations
  if (n0 > 0) {
    dens_eval0 <- lapply(1:k, function(j) {
      dmvn(X = split_data$data_comp,
           mu = means[[j]],
           sigma = covs[[j]]
           ) * pi[j]
      }
    )
    dens_eval0 <- do.call(cbind, dens_eval0)

    # Normalize by row.
    gamma0 <- aaply(
      .data = dens_eval0, 
      .margins = 1, 
      .fun = function(x) {x / sum(x)}, 
      .drop = FALSE
    )
    
    # Format
    colnames(dens_eval0) <- colnames(gamma0) <- paste0("k", seq(1:k))
    rownames(dens_eval0) <- rownames(gamma0) <- seq(1:n0)
    out$dens_eval0 <- dens_eval0
    out$gamma0 <- gamma0
  }

  ## Incomplete observations
  if (n1 > 0) {
    data_incomp <- split_data$data_incomp
    rownames(data_incomp) = NULL
    dens_eval1 <- aaply(
      .data = data_incomp,
      .margins = 1,
      .fun = function(x){Responsibility.eval_dens_incomp(x, means, covs, pi)},
      .drop = FALSE
    )

    # Normalize
    gamma1 <- aaply(
      .data = dens_eval1, 
      .margins = 1, 
      .fun = function(x) {x / sum(x)}, 
      .drop = FALSE
    )
    
    # Format
    colnames(dens_eval1) <- colnames(gamma1) <- paste0("k", seq(1:k))
    rownames(dens_eval1) <- rownames(gamma1) <- seq(1:n1)
    out$dens_eval1 <- dens_eval1
    out$gamma1 <- gamma1
  }
  
  return(out)
}


#-------------------------------------------------------------------------------
# Working Response Vectors. 
#-------------------------------------------------------------------------------

#' Generate Working Response
#' 
#' @param y Vector with missing elements.
#' @param mean Numeric mean.
#' @param cov Numeric covariance.
#' @param gamma Numeric vector of responsibilities.
#' @return Numeric working response vector. 

CalcWorkResp <- function(
  y,
  mean,
  cov,
  gamma
) {
  
  # Current observation
  is_mis <- is.na(y)
  is_obs <- !is_mis
  
  # Permutation
  perm <- c(which(is_obs), which(is_mis))
  rev_perm <- order(perm)
  
  # Partition covariance
  cov_mis_obs <- cov[is_mis, is_obs, drop = FALSE]
  var_mis_mis <- cov[is_obs, is_obs, drop = FALSE]
  var_mis_mis_inv <- matInv(var_mis_mis)
  
  # Observed components
  obs_ele <- matrix(y[is_obs], ncol = 1)
  
  # Conditional expectation of missing compoments
  mis_ele <- mean[is_mis] + MMP(cov_mis_obs, MMP(var_mis_mis_inv, obs_ele - mean[is_obs]))
  
  # Working response
  working_response <- rbind(obs_ele, mis_ele)
  
  # Return
  return(gamma * working_response[rev_perm])
}


#' Working Response Vectors
#'
#' Calculate the working response vectors.
#'
#' @param data_incomp Incomplete observations.
#' @param mean Numeric mean.
#' @param cov Numeric covariance.
#' @param gamma Numeric vector of responsibilities.
#' @return Numeric vector, the responsibility-weighted cumulative working
#'   response vector.

WorkResp <- function(
  data_incomp, 
  mean, 
  cov, 
  gamma = NULL) {
  
  # Dimensions
  n <- nrow(data_incomp)

  # Reponsibilities
  if (is.null(gamma)) {
    gamma <- rep(1, n)
  }

  # Body of "for loop"
  aux <- function(i){CalcWorkResp(data_incomp[i, ], mean, cov, gamma[i])}

  # Loop over observations
  loop <- lapply(seq(1:n), aux)
  out <- do.call(rbind, loop)
  return(out)
}


#-------------------------------------------------------------------------------
# Expected Residual Outer Product. 
#-------------------------------------------------------------------------------

#' Expected Residual Outer Product
#'
#' Calculates the expected residual outer product.
#'
#' @param data_incomp Data for observations with missingness.
#' @param new_mean New mean.
#' @param old_mean Initial mean.
#' @param old_cov Initial covariance.
#' @param gamma Responsibilities.
#'
#' @return Numeric matrix, the responsibility-weighted, cumulative,
#'  expected residual outer product.

ExpResidOP <- function(
  data_incomp, 
  new_mean, 
  old_mean, 
  old_cov, 
  gamma = NULL) {
  
  # Dimensions
  d <- ncol(data_incomp)
  n <- nrow(data_incomp)

  # Reponsibilities
  if (is.null(gamma)) {
    gamma <- rep(1, n)
  }

  # Body of "for loop"
  aux <- function(i) {
    # Current observation
    y <- data_incomp[i, ]
    is_mis <- is.na(y)
    is_obs <- !is_mis

    # Permutation
    perm <- c(which(is_obs), which(is_mis))
    rev_perm <- order(perm)
    n_obs <- sum(is_obs)

    # Partition covariance
    var_mis_mis <- old_cov[is_mis, is_mis, drop = FALSE]
    cov_mis_obs <- old_cov[is_mis, is_obs, drop = FALSE]
    var_obs_obs <- old_cov[is_obs, is_obs, drop = FALSE]

    # Inverses
    var_obs_obs_inv <- matInv(var_obs_obs)
    pre_mis_mis_inv <- SchurC(var_mis_mis, var_obs_obs, cov_mis_obs)

    # Observed components
    obs_ele <- matrix(y[is_obs], ncol = 1)

    # Conditional expectation of missing compoments
    miss_ele <- old_mean[is_mis] + MMP(cov_mis_obs, MMP(var_obs_obs_inv, obs_ele - old_mean[is_obs]))

    # Working response
    working_response <- rbind(obs_ele, miss_ele)

    # Residual
    residual <- working_response - new_mean[perm]

    # Residual outer product
    resid_OP <- matOP(residual, residual)

    # Add correction
    idx <- seq(from = n_obs + 1, to = d)
    resid_OP[idx, idx] <- resid_OP[idx, idx] + pre_mis_mis_inv

    # Recover initial order
    resid_OP <- resid_OP[rev_perm, rev_perm]

    # Return
    return(gamma[i] * resid_OP)
  }

  loop <- lapply(seq(1:n), aux)
  out <- Reduce("+", loop)
  
  # Output
  return(out)
}
