# Purpose: Fitting function for fitting a single multivariate normal
# distribution in the presence of missingness. 
# Updated: 2021-07-24

#-------------------------------------------------------------------------------

#' Estimation for a Multivariate Normal Distribution with No Missingness.
#'
#' @param data Numeric data matrix.
#' @param init_mean Optional initial mean vector.
#' @param fix_mean Fix the means to their starting value? Must initialize
#' @return An object of class \code{mvn}.
#' @noRd

FitMVNComplete <- function(
  data,
  init_mean = NULL,
  fix_mean = FALSE
) {
  
  # Dimensions.
  n <- nrow(data)
  d <- ncol(data)
  orig_col_names <- colnames(data)
  
  # Mean.
  if (fix_mean) {
    new_mean <- init_mean
  } else {
    new_mean <- apply(data, 2, mean)
  }
  
  # Covariance.
  new_cov <- matCov(data, data)
  new_cov_inv <- matInv(new_cov)
  dimnames(new_cov) <- list(orig_col_names, orig_col_names)
  
  # Objective (log-likelihood).
  n <- nrow(data)
  mean_mat <- matrix(data = new_mean, nrow = n, ncol = d, byrow = TRUE)
  resid <- data - mean_mat
  objective <- -n * log(det(new_cov)) - tr(MMP(matInv(new_cov), matIP(resid, resid)))
  
  # Output.
  out <- methods::new(
    Class = "mvn",
    Completed = data,
    Covariance = new_cov,
    Data = data,
    Mean = new_mean,
    Objective = objective
  )
  return(out)
}


#-------------------------------------------------------------------------------

#' Parameter Initialization for Multivariate Normal Distribution with Missingness 
#' 
#' @param data_comp Complete cases. 
#' @param init_mean Optional initial mean vector.
#' @param init_cov Optional initial covariance matrix.
#' @return List containing the initialized `mean` and `cov`. 
#' @noRd

MVNMissInit <- function(
  data_comp,
  init_mean,
  init_cov
) {
  
  # Output structure. 
  theta0 <- list()
  n <- nrow(data_comp)
  
  # Case 1: Both init_means and init_cov provided.
  if (!is.null(init_mean) & !is.null(init_cov)) {
    
    theta0$mean <- init_mean
    theta0$cov <- init_cov
    
  # Case 2: At least one of init_means or S0 is null
  } else {
    
    # Check for complete obs
    if (n == 0) {
      stop("If there are no complete cases, initial values are required for all parameters.")
    }
    
    # Initialize mean if null. 
    if (is.null(init_mean)) {
      theta0$mean <- apply(data_comp, 2, mean)
    } else {
      theta0$mean <- init_mean
    }
    
    # Initialize covariance if null.
    if (is.null(init_cov)) {
      theta0$cov <- matCov(data_comp, data_comp)
    } else {
      theta0$cov <- init_cov
    }
    
  }
  
  return(theta0)
}


#-------------------------------------------------------------------------------

#' Parameter Update for a Multivariate Normal Distribution with Missingness
#' 
#' @param split_data Data partitioned by missingness.
#' @param theta List containing the current `mean` and `cov`. 
#' @param fix_mean Fix the mean to its starting value? Must initialize. 
#' @return List containing:
#' \itemize{
#'   \item The updated `mean` and `cov`.
#'   \item The initial `old_obj` and final `new_obj` EM objective. 
#'   \item The increase in the EM objective `delta`. 
#' }
#' @noRd

MVNMissUpdate <- function(
  split_data,
  theta,
  fix_mean = FALSE
) {
  
  # Unpack
  data_comp <- split_data$data_comp
  n0 <- split_data$n0
  d <- split_data$n_col
  
  data_incomp <- split_data$data_incomp
  n1 <- split_data$n1
  n <- n0 + n1
  
  old_mean <- theta$mean
  old_cov <- theta$cov
  
  # Initial objective.
  resid_cov <- array(0, c(d, d))
  
  ## Complete cases.
  if (n0 > 0){
    mean_mat <- matrix(data = old_mean, nrow = n0, ncol = d, byrow = TRUE)
    resid_mat <- data_comp - mean_mat
    resid_cov <- resid_cov + matIP(resid_mat, resid_mat)
  }
  
  ## Incomplete cases.
  if (n1 > 0){
    resid_cov <- resid_cov + ExpResidOP(
      data_incomp,
      old_mean,
      old_mean,
      old_cov
    )
  }
 
  ## Initial EM objective. 
  old_obj <- (-n) * log(det(old_cov)) - tr(MMP(matInv(old_cov), resid_cov))
  
  # Update mean:
  if (fix_mean) {
    new_mean <- old_mean
  } else {
    new_mean <- 0
    
    ## Complete cases.
    if (n0 > 0) {
      new_mean <- new_mean + apply(data_comp, 2, sum)
    }
    
    ## Incomplete cases.
    if (n1 > 0) {
      working_response <- WorkResp(data_incomp, old_mean, old_cov)
      new_mean <- new_mean + apply(working_response, 2, sum)
    }
    
    ## Normalize mean vector.
    new_mean <- (new_mean / n)
  }
  
  # Update covariance.
  new_cov <- array(0, c(d, d))
  
  ## Complete observations.
  if (n0 > 0) {
    mean_mat <- matrix(data = new_mean, nrow = n0, ncol = d, byrow = TRUE)
    resid_mat <- data_comp - mean_mat
    resid_cov <- matIP(resid_mat, resid_mat)
    new_cov <- new_cov + resid_cov
  }
  
  ## Incomplete observations.
  if (n1 > 0) {
    new_cov <- new_cov + ExpResidOP(
      data_incomp,
      new_mean,
      old_mean,
      old_cov
    )
  }
  
  ## Normalize covariance matrix.
  new_cov <- (new_cov / n)
  dimnames(new_cov) <- list(split_data$orig_col_names, split_data$orig_col_names)
  
  ## Final EM objective.
  new_obj <- (-n) * log(det(new_cov)) - (d * n)
  
  # Increment.
  delta <- new_obj - old_obj
  
  # Output.
  out <- list()
  out$mean <- new_mean
  out$cov <- new_cov
  out$new_obj <- new_obj
  out$old_obj <- old_obj
  out$delta <- delta
  return(out)
}


#-------------------------------------------------------------------------------

#' Maximization Iteration.
#'
#' Used by both \code{\link{FitMVN}} and \code{\link{FitMix}}.
#'
#' @param theta0 List containing the initial values of `mean` and `cov`.
#' @param Update Function to iteratively update theta0.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return List containing the updated `mean` and `cov`, the final EM objective
#'   `new_obj`. 
#' @noRd

Maximization <- function(
  theta0,
  Update,
  maxit,
  eps,
  report 
) {
  
  # Update until max iterations or tolerance are reached. 
  for (i in 1:maxit) {
    
    # Update:
    theta1 <- Update(theta0)
   
    # Accept if increment in EM objective is positive.
    if (theta1$delta > 0) {
      theta0 <- theta1
      if (report) {
        cat("Objective increment: ", signif(theta1$d, digits = 3), "\n")
      }
    }
    
    # Terminate if increment is below tolerance.
    if (theta1$delta < eps) {
      # If EM failes to perform any updates, e.g. the iteration is initialized
      # at the MLE, keep initial objective.
      if (i == 1) {
        theta0$new_obj <- theta1$old_obj
      }
      break
    }
  }
  
  # Fitting report
  if (report) {
    if (i < maxit) {
      cat(paste0(i - 1, " update(s) performed before reaching tolerance limit.\n\n"))
    } else {
      cat(paste0(i, " update(s) performed without reaching tolerance limit.\n\n"))
    }
  }
  
  # Output.
  return(theta0)
}


#-------------------------------------------------------------------------------

#' Imputation for a Multivariate Normal Distribution with Missingness
#' 
#' @param split_data Data partitioned by missingness.
#' @param theta List containing the `mean` and `cov`.
#' @return Data.matrix, in the same order as the original data, with missing values
#'   imputed to their expectations. 
#' @noRd

MVNMissImpute <- function(
  split_data,
  theta
) {
  
  d <- split_data$n_col 
  out <- matrix(NA, nrow = 0, ncol = d)
  
  # Complete data.
  n0 <- split_data$n0
  if (n0 > 0) {
    out <- rbind(out, split_data$data_comp)
  }
  
  # Impute data.
  n1 <- split_data$n1
  if (n1 > 0) {
    data_imp <- WorkResp(split_data$data_incomp, theta$mean, theta$cov)
    out <- rbind(out, data_imp)
  }
  
  # Empty data.
  n2 <- split_data$n2
  if (n2 > 0) {
    data_empty <- matrix(data = theta$mean, nrow = n2, ncol = d, byrow = TRUE)
    out <- rbind(out, data_empty)
  }
  
  # Restore initial order. 
  init_order <- split_data$init_order
  out <- out[order(init_order), , drop = FALSE]
  
  # Output
  rownames(out) <- split_data$orig_row_names
  colnames(out) <- split_data$orig_col_names
  return(out)
}


#-------------------------------------------------------------------------------

#' Estimation for a Multivariate Normal Distribution with Missingness
#'
#' @param data Numeric data matrix.
#' @param init_mean Optional initial mean vector.
#' @param fix_mean Fix the means to their starting value? Must initialize.
#' @param init_cov Optional initial covariance matrix.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return An object of class \code{mvn}.
#' @noRd

FitMVNMiss <- function(
  data,
  init_mean,
  fix_mean,
  init_cov,
  maxit,
  eps,
  report 
) {
  
  # Partitioned data
  split_data <- PartitionData(data)
  
  # Initialize
  theta0 <- MVNMissInit(
    split_data$data_comp,
    init_mean,
    init_cov
  )
  
  # Maximization
  Update <- function(theta) {MVNMissUpdate(split_data, theta, fix_mean)}
  theta1 <- Maximization(theta0, Update, maxit, eps, report)
  
  # Imputation.
  imputed <- MVNMissImpute(split_data, theta1)
  
  # Output.
  out <- methods::new(
    Class = "mvn",
    Completed = imputed,
    Covariance = theta1$cov,
    Data = data,
    Mean = theta1$mean,
    Objective = theta1$new_obj
  )
  return(out)
}


#-------------------------------------------------------------------------------
# Main Function
#-------------------------------------------------------------------------------

#' Fit Multivariate Normal Distribution
#'
#' Given a matrix of n x d-dimensional random vectors, possibly containing
#' missing elements, estimates the mean and covariance of the best fitting
#' multivariate normal distribution.
#'
#' @param data Numeric data matrix.
#' @param init_mean Optional initial mean vector.
#' @param fix_mean Fix the mean to its starting value? Must initialize. 
#' @param init_cov Optional initial covariance matrix.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @importFrom stats complete.cases
#' @return An object of class \code{mvn}.

FitMVN <- function(
  data, 
  init_mean = NULL, 
  fix_mean = FALSE, 
  init_cov = NULL, 
  maxit = 100, 
  eps = 1e-6, 
  report = TRUE
  ) {
  
  # Check for missingness.
  is_mis <- sum(is.na(data)) > 0

  # Case of no missingness.
  if (!is_mis) {
    out <- FitMVNComplete(
      data,
      init_mean,
      fix_mean
    )
    
  # Case with missingness.
  } else {
    out <- FitMVNMiss(
      data,
      init_mean,
      fix_mean,
      init_cov,
      maxit,
      eps,
      report
    )
  }
  return(out)
}
