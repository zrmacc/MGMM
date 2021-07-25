# Purpose: Fits a multivariate normal mixture in the presence of missingness.
# Updated: 2021-07-24

#------------------------------------------------------------------------------

#' Parameter Initialization for a Gaussian Mixture Model
#' 
#' @param split_data Data partitioned by missingness.
#' @param k Number of mixture components. 
#' @param init_means Optional list of initial mean vectors.
#' @param init_covs Optional list of initial covariance matrices.
#' @param init_props Optional vector of initial cluster proportions.
#' @return List containing initial parameter values.
#' @noRd

MixInit <- function(
  split_data,
  k, 
  init_means,
  init_covs,
  init_props
) {
  
  # Unpack.
  theta0 <- list()
  n0 <- split_data$n0
  d <- split_data$n_col
  data_comp <- split_data$data_comp
  
  # Case 1: All parameters provided.
  if (all(!is.null(init_means), !is.null(init_covs), !is.null(init_props))) {
    
    theta0$means <- init_means
    theta0$covs <- init_covs
    theta0$pi <- init_props
    
  # Case 2: Initial values partially missing.
  } else {
    
    if (n0 == 0) {
      stop("If no observations are complete, initial values are required for all parameters.")
    }
    
    # If complete cases are available, apply kmeans.
    k_means <- stats::kmeans(data_comp, k, iter.max = 100, nstart = 100)
    
    # Cluster assignments.
    cluster_assignment <- k_means$cluster
    
    # Initialize means.
    if (is.null(init_means)) {
      means <- k_means$centers
      theta0$means <- lapply(seq_len(k), function(i){means[i, ]})
    } else {
      theta0$means <- init_means
    }
    
    # Initialize covariances.
    if (is.null(init_covs)) {
      theta0$covs <- lapply(seq_len(k), function(i) {
        clust <- data_comp[cluster_assignment == i, , drop = FALSE]
        return(matCov(clust, clust))
      })
    } else {
      theta0$covs <- init_covs
    }
    
    # Initialize proportions
    if (is.null(init_props)) {
      theta0$pi <- as.numeric(table(cluster_assignment)) / n0
    } else {
      theta0$pi <- init_props
    }
  } 
  
  # Check that estimated covariances are positive definite. 
  eigen_values <- unlist(lapply(theta0$covs, FUN = eigSym))
  
  if (min(eigen_values) <= 0) {
    stop("Initial covariance matrices are not all positive definite.")
  }
  
  # Initial responsibilities.
  theta0$gamma <- Responsibility(
    split_data, 
    theta0$means, 
    theta0$covs, 
    theta0$pi
  )
  
  # Output.
  return(theta0)
}

#------------------------------------------------------------------------------

#' Cluster Sizes for a Gaussian Mixture Model
#' 
#' Sum of the cluster responsibilities across examples.
#' 
#' @param split_data Data partitioned by missingness.
#' @param gamma List of component responsibilities. 
#' @return K x 1 numeric vector.
#' @noRd 

MixClusterSizes <- function(
  split_data,
  gamma
) {
  
  # Unpack.
  n0 <- split_data$n0
  n1 <- split_data$n1
  k <- gamma$k
  
  # Cluster sizes.
  cluster_sizes <- rep(0, k)
  
  ## Complete cases.
  if (n0 > 0) {
    gamma0 <- gamma$gamma0
    cluster_sizes <- cluster_sizes + apply(gamma0, 2, sum)
  }
  
  ## Incomplete cases. 
  if (n1 > 0) {
    gamma1 <- gamma$gamma1
    cluster_sizes <- cluster_sizes + apply(gamma1, 2, sum)
  }
  
  # Output.
  return(cluster_sizes)
}


#------------------------------------------------------------------------------

#' Expected Residual Outer Product for a Gaussian Mixture Model
#' 
#' @param split_data Data partitioned by missingness.
#' @param new_means List of updated means.
#' @param old_means List of previous means. 
#' @param covs List of component covariances.
#' @param gamma List of component responsibilities. 
#' @return List of k responsibility-weighted expected residual outer products. 
#' @noRd

MixResidOP <- function(
  split_data,
  new_means,
  old_means,
  covs,
  gamma
) {
  
  # Unpack.
  n0 <- split_data$n0
  n1 <- split_data$n1
  d <- split_data$n_col 
  k <- gamma$k
  
  # Loop over mixture components. 
  out <- lapply(seq_len(k), function(j) {
    
    resid_op <- array(0, dim = c(d, d))
    
    ## Complete cases.
    if (n0 > 0) {
      
      # Residuals.
      mean_mat <- matrix(data = new_means[[j]], nrow = n0, ncol = d, byrow = TRUE)
      resid <- split_data$data_comp - mean_mat
      
      # Responsibility-weighted OP.
      resid_op <- resid_op + matIP(resid, gamma$gamma0[, j] * resid)
    }
    
    ## Incomplete cases.
    if (n1 > 0) {
      
      # Responsibility-weighted OP
      resid_op <- resid_op + ExpResidOP(
        split_data$data_incomp, 
        new_means[[j]], 
        old_means[[j]], 
        covs[[j]], 
        gamma$gamma1[, j])
    }
    
    # Return residual outer product. 
    return(resid_op)
  })
  
  return(out)
}


#------------------------------------------------------------------------------

#' EM Objective for a Gaussian Mixture Model
#' 
#' @param cluster_sizes Cluster sizes.
#' @param pi Cluster proportions
#' @param covs List of component covariances. 
#' @param resid_ops List of residual outer products.
#' @return Numeric objective value.
#' @noRd

MixEMObj <- function(
  cluster_sizes,
  pi,
  covs,
  resid_ops
) {
  
  # Pi term.
  k <- length(pi)
  pi_term <- sum(cluster_sizes * log(pi))
  
  # Determinant term.
  det_term <- lapply(1:k, function(j) {
    cluster_sizes[j] * log(det(covs[[j]]))
  })
  det_term <- do.call(sum, det_term)
  
  # Trace term.
  trace_term <- lapply(1:k, function(j) {
    tr(MMP(matInv(covs[[j]]), resid_ops[[j]]))
  })
  trace_term <- do.call(sum, trace_term)
  
  # Objective.
  obj <- pi_term - det_term - trace_term
  return(obj)
}


#------------------------------------------------------------------------------

#' Mean Update for Mixture of MVNs with Missingness.
#' 
#' @param split_data Data partitioned by missingness.
#' @param means List of component means.
#' @param covs List of component covariances.
#' @param gamma List of component responsibilities. 
#' @return List containing the updated component means. 

MixUpdateMeans <- function(
  split_data,
  means,
  covs,
  gamma
) {
  
  # Unpack.
  n0 <- split_data$n0
  n1 <- split_data$n1
  k <- length(means)
  
  # Cluster sizes.
  cluster_sizes <- MixClusterSizes(
    split_data,
    gamma
  )

  # Update means.
  new_means <- lapply(seq_len(k), function(j) {
    total <- 0
    
    ## Complete cases.
    if (n0 > 0) {
      total <- total + apply(gamma$gamma0[, j] * split_data$data_comp, 2, sum)
    }
    
    ## Incomplete cases.
    if (n1 > 0) {
      working_response <- WorkResp(
        split_data$data_incomp,
        means[[j]],
        covs[[j]],
        gamma$gamma1[, j]
      )
      total <- total + apply(working_response, 2, sum)
    }
    
    # Normalize mean vector.
    new_mean <- (total) / cluster_sizes[j]
    names(new_mean) <- split_data$orig_col_names
    return(new_mean)
  })
  
  return(new_means)
}


#------------------------------------------------------------------------------

#' Parameter Update for Mixutre of MVNs with Missingness. 
#' 
#' @param split_data Data partitioned by missingness.
#' @param theta List containing the current `means`, `covs`, `pi`, and `gamma`.
#' @param fix_means Fix the mean to its starting value? Must initialize. 
#' @return List containing:
#' \itemize{
#'   \item The updated `mean`, `cov`, `pi`, and `gamma`.
#'   \item The initial `old_obj` and final `new_obj` EM objective. 
#'   \item The increase in the EM objective `delta`. 
#' }
#' @noRd 

MixUpdate <- function(
  split_data,
  theta,
  fix_means
) {
  
  # Previous parameters.
  old_means <- theta$means
  old_covs <- theta$covs
  old_pi <- theta$pi
  old_gamma <- theta$gamma
  
  # Cluster sizes.
  old_cluster_sizes <- MixClusterSizes(
    split_data,
    old_gamma
  )  
  
  # Old residual outer products.
  old_resid_ops <- MixResidOP(
    split_data,
    old_means,
    old_means,
    old_covs,
    old_gamma
  )
  
  # Initial objective.
  old_obj <- MixEMObj(
    old_cluster_sizes,
    old_pi,
    old_covs,
    old_resid_ops
  )

  # Update means.
  if (fix_means){
    new_means <- old_means
  } else {
    new_means <- MixUpdateMeans(
      split_data,
      old_means,
      old_covs,
      old_gamma 
    )
  }
  
  # Update covariances.
  ## Update outer products
  new_resid_ops <- MixResidOP(
    split_data,
    new_means,
    old_means,
    old_covs,
    old_gamma
  )
  
  ## Normalize covariance matrices. 
  k <- theta$gamma$k
  new_covs <- lapply(seq_len(k), function(j) {
    new_cov <- new_resid_ops[[j]] / old_cluster_sizes[[j]]
    dimnames(new_cov) <- list(split_data$orig_col_names, split_data$orig_col_names)
    return(new_cov)
  })
  
  ## Update responsibilities
  new_gamma <- Responsibility(
    split_data,
    new_means,
    new_covs,
    old_pi
  )
  
  # Update cluster proportions.
  new_pi <- old_cluster_sizes / sum(old_cluster_sizes)
  
  # New EM objective. 
  new_obj <- MixEMObj(
    old_cluster_sizes,
    new_pi,
    new_covs,
    new_resid_ops
  )
  
  # Objective increment.
  delta <- new_obj - old_obj
  
  # Output.
  out <- list()
  out$means <- new_means
  out$covs <- new_covs
  out$pi <- new_pi
  out$gamma <- new_gamma
  out$new_obj <- new_obj
  out$old_obj <- old_obj
  out$delta <- delta
  return(out)
}


#------------------------------------------------------------------------------

#' Cluster Assignment for Mixutre of MVNs with Missingness. 
#' 
#' @param split_data Data partitioned by missingness.
#' @param theta List containing the current `means`, `covs`, `pi`, and `gamma`.
#' @return List containing:
#' \itemize{
#'   \item Matrix of cluster `Assignments`.
#'   \item Matrix of `Density` evaluations. 
#'   \item Matrix of cluster `Responsibilities`.
#' }
#' @noRd

MixClusterAssign <- function(
  split_data,
  theta
) {
  
  # Unpack.
  n2 <- split_data$n2
  d <- split_data$n_col 
  k <- theta$gamma$k
  
  # Responsibilities.
  resp <- rbind(theta$gamma$gamma0, theta$gamma$gamma1)
  
  # Density evaluations.
  dens <- rbind(theta$gamma$dens_eval0, theta$gamma$dens_eval1)
  
  # Assignments.
  map_assign <- apply(resp, 1, which.max)
  if (n2 > 0) {
    map_assign <- c(map_assign, rep(NA, n2))
  }
  
  # Recover initial order.
  init_order <- split_data$init_order
  map_assign <- map_assign[order(init_order)]
  names(map_assign) <- split_data$orig_row_names
  
  # Responsibilities.
  if (n2 > 0) {
    resp <- rbind(resp, array(NA, dim = c(n2, k)))
  }
  
  # Recover initial order.
  resp <- resp[order(init_order), ]
  rownames(resp) <- NULL
  
  # Entropy.
  entropy <- plyr::aaply(
    .data = resp,
    .margins = 1,
    .fun = function(x){-sum(x * log(x)) / log(k)}
  )
  
  # Density evaluations.
  if (n2 > 0) {
    dens <- rbind(dens, array(NA, dim = c(n2, k)))
  }
  
  # Recover initial order.
  dens <- dens[order(init_order), ]
  rownames(dens) <- NULL
  
  # Assignment matrix.
  assign <- cbind(
    Assignments = map_assign,
    Entropy = entropy
  )
  rownames(assign) <- split_data$orig_row_names
  
  # Responsibility matrix.
  rownames(resp) <- split_data$orig_row_names
  
  # Density matrix
  rownames(dens) <- split_data$orig_row_names
  
  # Output.
  out <- list()
  out$Assignments <- assign
  out$Responsibilities <- resp
  out$Density <- dens
  return(out)
}


#------------------------------------------------------------------------------

#' Imputation for Gaussian Mixture Models.
#' 
#' @param split_data Data partitioned by missingness.
#' @param theta List containing the current `means`, `covs`, `pi`, and `gamma`.
#' @return Numeric matrix in the same order as the original data, with missing values
#'   imputed to their expectations. 
#' @noRd

MixImpute <- function(
  split_data,
  theta
) {
  
  # Unpack.
  n0 <- split_data$n0
  n1 <- split_data$n1
  n2 <- split_data$n2
  d <- split_data$n_col 
  k <- theta$gamma$k
  
  # Output structure. 
  out <- matrix(NA, nrow = 0, ncol = d)
  
  ## Complete cases. 
  if (n0 > 0) {
    out <- rbind(out, split_data$data_comp)
  }
  
  ## Incomplete cases. 
  if (n1 > 0) {
    data_imp <- lapply(seq_len(k), function(j) {
      working_response <- WorkResp(
        split_data$data_incomp,
        theta$means[[j]],
        theta$covs[[j]],
        theta$gamma$gamma1[, j]
      )
      return(working_response)
    })
    data_imp <- Reduce("+", data_imp)
    out <- rbind(out, data_imp)
  }
  
  ## Empty cases. 
  if (n2 > 0) {
    data_imp <- lapply(seq_len(k), function(j) {
      return(theta$means[[j]] * theta$pi[j])
    })
    data_imp <- Reduce("+", data_imp)
    data_imp <- matrix(data = data_imp, nrow = n2, ncol = d, byrow = TRUE)
    out <- rbind(out, data_imp)
  }
  
  # Output.
  init_order <- split_data$init_order
  out <- out[order(init_order), ]
  rownames(out) <- split_data$orig_row_names
  colnames(out) <- split_data$orig_col_names
  return(out)
}

#------------------------------------------------------------------------------
# Main Function
#------------------------------------------------------------------------------

#' Fit Multivariate Mixture Distribution
#'
#' Given a matrix of random vectors, estimates the parameters for a mixture of
#' multivariate normal distributions. Accommodates arbitrary patterns of
#' missingness, provided the elements are missing at random (MAR).
#'
#' @param data Numeric data matrix.
#' @param k Number of mixture components. Defaults to 2.
#' @param init_means Optional list of initial mean vectors.
#' @param fix_means Fix means to their starting values? Must initialize. 
#' @param init_covs Optional list of initial covariance matrices.
#' @param init_props Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return Object of class \code{mix}.

FitMix <- function(
  data, 
  k = 2, 
  init_means = NULL, 
  fix_means = FALSE, 
  init_covs = NULL, 
  init_props = NULL, 
  maxit = 100, 
  eps = 1e-6, 
  report = FALSE
) {

  # Partition data.
  split_data <- PartitionData(data)

  # Initialization.
  theta0 <- MixInit(split_data, k, init_means, init_covs, init_props)

  # Maximzation.
  Update <- function(theta) {MixUpdate(split_data, theta, fix_means)}
  theta1 <- Maximization(theta0, Update, maxit, eps, report)

  # Cluster assignments.
  assign <- MixClusterAssign(
    split_data,
    theta1
  )
  
  # Imputation.
  imputed <- MixImpute(split_data, theta1)

  # Output
  out <- methods::new(
    Class = "mix", 
    Assignments = assign$Assignments, 
    Completed = imputed,
    Components = k, 
    Covariances = theta1$covs, 
    Density = assign$Density, 
    Means = theta1$means, 
    Objective = theta1$new_obj,
    Proportions = theta1$pi, 
    Responsibilities = assign$Responsibilities
  )
  return(out)
}
