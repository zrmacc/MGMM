# Purpose: Calculate E-step expectations.
# Updated: 2021-07-24

#-------------------------------------------------------------------------------
# Responsibilities
#-------------------------------------------------------------------------------

#' Evaluate the Density of an Incomplete Observation
#'
#' Computes the marginal density of the observed elements under the mixture
#' model for a single observation vector that may contain missing values.
#'
#' @param y Vector with missing elements.
#' @param means List of mean vectors.
#' @param covs List of covariance matrices.
#' @param pi Vector of cluster proportions.
#' @return Numeric density of observed elements.
#' @noRd

EvalDensIncompObs <- function(
  y, 
  means,
  covs,
  pi
) {
  
  # Mixture components. 
  k <- length(pi)
  
  # Observed elements.
  is_obs <- !is.na(y)
  obs_ele = y[is_obs]
  
  # Density evaluation of observed elements.
  obs_dens <- lapply(seq_len(k), function(j) {
    obs_mean <- means[[j]][is_obs]
    obs_cov <- covs[[j]][is_obs, is_obs]
    out <- mvnfast::dmvn(X = obs_ele, mu = obs_mean, sigma = obs_cov) * pi[j]
    return(out)
  })
  
  obs_dens <- unlist(obs_dens)
  return(obs_dens)
}


#' Responsibilities
#'
#' Posterior probability of cluster j for observation i given observed data;
#' paper eq. (3) and (8): \eqn{\gamma_{ij} = f(y_{obs}|\mu_j, \Sigma_j) \pi_j / \sum_l f(y_{obs}|\mu_l, \Sigma_l) \pi_l}{gamma_ij = f(y_obs|mu_j, Sigma_j) pi_j / sum_l f(y_obs|mu_l, Sigma_l) pi_l}.
#'
#' @param split_data Data partitioned by missingness.
#' @param means List of mean vectors.
#' @param covs List of covariance matrices.
#' @param pi Vector of cluster proportions.
#' @return List containing:
#' \itemize{
#'   \item k Number of mixture components.
#'   \item Density evaluations `dens_eval0` and responsibilities `gamma0` 
#'     for complete cases.
#'   \item Density evaluations `dens_eval1` and responsibilities `gamma1` 
#'     for incomplete cases. 
#' }
#' @noRd

Responsibility <- function(
  split_data,
  means, 
  covs, 
  pi
) {

  # Number of complete and incomplete observations.
  n0 <- split_data$n0
  n1 <- split_data$n1
  
  # Clusters.
  k <- length(pi)

  # Output structure.
  out <- list()
  out$k <- k

  if (n0 > 0) {
    dens_eval0 <- do.call(cbind, lapply(seq_len(k), function(j) {
      mvnfast::dmvn(X = split_data$data_comp, mu = means[[j]], sigma = covs[[j]]) * pi[j]
    }))
    gamma0 <- sweep(dens_eval0, 1L, rowSums(dens_eval0), "/")
    colnames(dens_eval0) <- colnames(gamma0) <- paste0("k", seq_len(k))
    rownames(dens_eval0) <- rownames(gamma0) <- seq_len(n0)
    out$dens_eval0 <- dens_eval0
    out$gamma0 <- gamma0
  }

  if (n1 > 0) {
    data_incomp <- split_data$data_incomp
    rownames(data_incomp) <- NULL
    dens_eval1 <- plyr::aaply(
      .data = data_incomp, .margins = 1L,
      .fun = function(x) EvalDensIncompObs(x, means, covs, pi),
      .drop = FALSE
    )
    gamma1 <- sweep(dens_eval1, 1L, rowSums(dens_eval1), "/")
    colnames(dens_eval1) <- colnames(gamma1) <- paste0("k", seq_len(k))
    rownames(dens_eval1) <- rownames(gamma1) <- seq_len(n1)
    out$dens_eval1 <- dens_eval1
    out$gamma1 <- gamma1
  }
  
  return(out)
}


#-------------------------------------------------------------------------------
# Working Response Vectors. 
#-------------------------------------------------------------------------------

#' Working Response Vector
#'
#' Working response \eqn{\hat{y}_{ij}}{y_ij} per paper eq. (4): (y_obs, E[y_miss | z_ij=1, y_obs]),
#' with \eqn{E[y_{miss} | \ldots] = \mu_{miss} + \Sigma_{miss,obs} \Sigma_{obs,obs}^{-1} (y_{obs} - \mu_{obs})}{E[y_miss | ...] = mu_miss + Sigma_miss,obs Sigma_obs,obs^{-1} (y_obs - mu_obs)}.
#' 
#' @param y Vector with missing elements.
#' @param mean Numeric mean.
#' @param cov Numeric covariance.
#' @param gamma Numeric scalar; responsibility (posterior probability) for
#'   this observation under the given component.
#' @return Numeric working response vector (observed elements unchanged;
#'   missing elements filled with conditional expectation). 
#' @noRd

WorkRespIndiv <- function(
  y,
  mean,
  cov,
  gamma
) {
  
  # Current observation.
  is_mis <- is.na(y)
  is_obs <- !is_mis
  if (!any(is_mis)) {return(y)}
  
  # Permutation to place observed elements first (paper: y_obs then y_miss).
  perm <- c(which(is_obs), which(is_mis))
  rev_perm <- order(perm)
  
  # Partition covariance: Sigma_obs,obs, Sigma_miss,obs per paper eq. (4).
  cov_mis_obs <- cov[is_mis, is_obs, drop = FALSE]
  var_obs_obs <- cov[is_obs, is_obs, drop = FALSE]
  var_obs_obs_inv <- matInv(var_obs_obs)
  
  # Observed components.
  obs_ele <- matrix(y[is_obs], ncol = 1)
  
  # Conditional expectation of missing given observed: mu_miss + Sigma_miss,obs Sigma_obs,obs^{-1} (y_obs - mu_obs).
  mis_ele <- mean[is_mis] + MMP(cov_mis_obs, MMP(var_obs_obs_inv, obs_ele - mean[is_obs]))
  
  # Working response.
  working_response <- rbind(obs_ele, mis_ele)
  
  # Output.
  out <- gamma * working_response[rev_perm]
  return(out)
}


#' Working Response Vectors
#'
#' Calculate the working response vectors for all incomplete examples.
#'
#' @param data_incomp Incomplete observations.
#' @param mean Numeric mean.
#' @param cov Numeric covariance.
#' @param gamma Numeric vector of responsibilities (one per row of
#'   \code{data_incomp}); if \code{NULL}, defaults to 1 for each row.
#' @return Numeric matrix with the same dimensions as \code{data_incomp};
#'   missing elements are filled with conditional expectations.
#' @noRd

WorkResp <- function(
  data_incomp, 
  mean, 
  cov, 
  gamma = NULL
) {
  
  # Dimensions.
  n <- nrow(data_incomp)

  # Responsibilities.
  if (is.null(gamma)) {
    gamma <- rep(1, n)
  }

  # Loop over observations.
  out <- lapply(seq_len(n), function(i) {
    WorkRespIndiv(data_incomp[i, ], mean, cov, gamma[i])
  })
  out <- do.call(rbind, out)
  return(out)
}


#-------------------------------------------------------------------------------
# Expected Residual Outer Product. 
#-------------------------------------------------------------------------------

#' Expected Residual Outer Product
#'
#' Working residual outer product per paper eq. (5): \eqn{\gamma_{ij}}{gamma_ij} [ \eqn{(\hat{y}_{ij} - \mu_j)(\hat{y}_{ij} - \mu_j)'}{(y_ij - mu_j)(y_ij - mu_j)'}
#' + (0 0; 0 \eqn{\Lambda^{-1}_{j,TT}}{Lambda^{-1}}) ], where \eqn{\Lambda^{-1}}{Lambda^{-1}} is the conditional covariance of missing
#' given observed (Schur complement).
#'
#' @param data_incomp Data for observations with missingness.
#' @param new_mean New mean \eqn{\mu_j}{mu_j} (used in residual).
#' @param old_mean Previous mean (used for conditional expectation).
#' @param old_cov Previous covariance (used for working response and Schur complement).
#' @param gamma Responsibilities \eqn{\gamma_{ij}}{gamma_ij}.
#' @return Numeric matrix, the responsibility-weighted sum of working residual outer products.
#' @noRd

ExpResidOP <- function(
  data_incomp, 
  new_mean, 
  old_mean, 
  old_cov, 
  gamma = NULL
) {
  d <- ncol(data_incomp)
  n <- nrow(data_incomp)
  if (is.null(gamma)) {
    gamma <- rep(1, n)
  }

  out <- lapply(seq_len(n), function(i) {
    y <- data_incomp[i, ]
    is_mis <- is.na(y)
    is_obs <- !is_mis
    perm <- c(which(is_obs), which(is_mis))
    rev_perm <- order(perm)
    n_obs <- sum(is_obs)

    var_mis_mis <- old_cov[is_mis, is_mis, drop = FALSE]
    cov_mis_obs <- old_cov[is_mis, is_obs, drop = FALSE]
    var_obs_obs <- old_cov[is_obs, is_obs, drop = FALSE]
    var_obs_obs_inv <- matInv(var_obs_obs)
    # Conditional covariance of missing given observed (Schur complement) = Lambda^{-1} in paper eq. (5).
    cond_cov_miss <- SchurC(var_mis_mis, var_obs_obs, cov_mis_obs)

    obs_ele <- matrix(y[is_obs], ncol = 1)
    miss_ele <- old_mean[is_mis] + MMP(cov_mis_obs, MMP(var_obs_obs_inv, obs_ele - old_mean[is_obs]))
    working_response <- rbind(obs_ele, miss_ele)
    residual <- working_response - new_mean[perm]
    resid_OP <- matOP(residual, residual)
    idx <- seq(from = n_obs + 1, to = d)
    resid_OP[idx, idx] <- resid_OP[idx, idx] + cond_cov_miss
    resid_OP <- resid_OP[rev_perm, rev_perm]
    gamma[i] * resid_OP
  })
  Reduce("+", out)
}
