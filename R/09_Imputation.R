# Purpose: Generate imputation from a fitted GMM.
# Updated: 2021-12-09

# -----------------------------------------------------------------------------

#' Impute Empty
#' 
#' Generates an imputation for the empty cases.
#' 
#' @param split_data Split
#' @param fit Fitted model.
#' @return Numeric matrix.
#' @importFrom methods is
#' @noRd

ImputeEmpty <- function(split_data, fit) {
  if (split_data$n2 > 0) {
    if (is(fit, "mvn")) {
      out <- rGMM(
        n = split_data$n2,
        d = split_data$n_col,
        k = 1, 
        means = fit@Mean,
        covs = fit@Covariance
      )
    } else {
      out <- rGMM(
        n = split_data$n2,
        d = split_data$n_col,
        k = fit@Components,
        pi = fit@Proportions,
        miss = 0.0,
        means = fit@Means,
        covs = fit@Covariances
      )
    }
  } else {
    out <- NULL
  }
  return(out)
}


# -----------------------------------------------------------------------------

#' Conditional Mean and Covariance
#' 
#' Calculate the conditional mean and covariance of the elements in 
#' index set A given the elements in index set B. 
#' 
#' @param y Outcome vector.
#' @param idx_a Elements whose distribution is sought.
#' @param idx_b Elements being conditioned on. 
#' @param mu Mean vector.
#' @param sigma Covariance matrix.
#' @return List containing the conditional mean and covariance.
#' @noRd

CalcCondDist <- function(y, idx_a, idx_b, mu, sigma) {
  
  # Split outcome.
  y_b <- y[idx_b]
  
  # Split mean.
  mu_a <- mu[idx_a]
  mu_b <- mu[idx_b]
  
  # Split covariance.
  sigma_aa <- sigma[idx_a, idx_a, drop = FALSE]
  sigma_ab <- sigma[idx_a, idx_b, drop = FALSE]
  sigma_bb <- sigma[idx_b, idx_b, drop = FALSE]
  
  # Calculate conditional mean and covariance.
  mu_cond <- mu_a - sigma_ab %*% solve(sigma_bb, y_b - mu_b)
  sigma_cond <- sigma_aa - sigma_ab %*% solve(sigma_bb, t(sigma_ab))
  
  out <- list(
    mu = mu_cond,
    sigma = sigma_cond
  )
  return(out)
}



# -----------------------------------------------------------------------------

#' Impute Incomplete
#' 
#' Generates an imputation for the incomplete cases. The imputation is
#' conditional on the observed data. 
#' 
#' @param split_data Split
#' @param fit Fitted model.
#' @return Numeric matrix.
#' @importFrom methods is
#' @noRd

ImputeIncomplete <- function(split_data, fit) {
  
  if (split_data$n1 > 0) {
    data_incomp <- split_data$data_incomp
    idx <- seq_len(nrow(data_incomp))
    
    if (is(fit, "mvn")) {
      
      sample_comps <- rep(1, length(idx))
      means <- list(fit@Mean)
      covs <- list(fit@Covariance)
      
    } else {
      
      # Sample mixture component with probability proportional 
      # to the responsibility.
      resp <- fit@Responsibilities[split_data$idx_incomp, ]
      split_resp <- split(resp, idx)
      sample_comps <- sapply(split_resp, function(x) {
        sample(seq_len(fit@Components), size = 1, prob = x)
      })
      
      means <- fit@Means
      covs <- fit@Covariances
      
    }
    
    # Impute.
    split_data_imputed <- lapply(idx, function(i) {
      z <- sample_comps[i]
      mu <- means[[z]]
      sigma <- covs[[z]]
      y <- data_incomp[i, ]
      
      # Conditional distribution.
      idx_a <- which(is.na(y))
      idx_b <- which(!is.na(y))
      cond_dist <- CalcCondDist(y, idx_a, idx_b, mu, sigma)
      y[is.na(y)] <- mvnfast::rmvn(
        n = 1,
        mu = cond_dist$mu,
        sigma = cond_dist$sigma
      )
      return(y)
    })
    out <- do.call(rbind, split_data_imputed)
  } else {
    out <- NULL
  }
  return(out)
}

# -----------------------------------------------------------------------------

#' Generate Imputation
#' 
#' Generates a stochastic imputation of a data set from a fitted data set.
#' 
#' @param fit Fitted model.
#' @return Numeric matrix with missing values imputed.
#' @export 
#' @examples
#' set.seed(100)
#' 
#' # Generate data and introduce missingness.
#' data <- rGMM(n = 25, d = 2, k = 1)
#' data[1, 1] <- NA
#' data[2, 2] <- NA
#' data[3, ] <- NA 
#' 
#' # Fit GMM.
#' fit <- FitGMM(data)
#' 
#' # Generate imputation.
#' imputed <- GenImputation(fit)

GenImputation <- function(fit) {
  
  # Check that missing data are present.
  data <- fit@Data
  any_na <- any(is.na(data))
  if (!any_na) {return(data)}
  
  # Split data by missingness pattern.
  split_data <- PartitionData(data)
  split_out <- split_data
  
  # Impute incomplete cases.
  split_out$data_incomp <- ImputeIncomplete(split_data, fit)
  
  # Impute empty cases.
  split_out$data_empty <- ImputeEmpty(split_data, fit)
  
  # Output.
  out <- ReconstituteData(split_out)
  return(out)
}

# -----------------------------------------------------------------------------

#' Combine Multiple Imputations
#' 
#' Combines point estimates and standard errors across multiple imputations.
#' 
#' @param points List of point estimates, potentially vector valued. 
#' @param covs List of sampling covariances, potentially matrix valued.
#' @return List containing the final point estimate (`point`) and 
#'   sampling covariance (`cov`).
#' @export
#' @examples 
#' set.seed(100)
#' 
#' # Generate data and introduce missingness.
#' data <- rGMM(n = 25, d = 2, k = 1)
#' data[1, 1] <- NA
#' data[2, 2] <- NA
#' data[3, ] <- NA 
#' 
#' # Fit GMM.
#' fit <- FitGMM(data)
#' 
#' # Lists to store summary statistics.
#' points <- list()
#' covs <- list()
#' 
#' # Perform 50 multiple imputations.
#' # For each, calculate the marginal mean and its sampling variance.
#' for (i in seq_len(50)) {
#'   imputed <- GenImputation(fit)
#'   points[[i]] <- apply(imputed, 2, mean)
#'   covs[[i]] <- cov(imputed) / nrow(imputed)
#' }
#' 
#' # Combine summary statistics across imputations.
#' results <- CombineMIs(points, covs)

CombineMIs <- function(points, covs) {
  m <- length(points)
  point <- Reduce("+", points) / m
  within_cov <- Reduce("+", covs) / m
  between_cov <- lapply(points, function(x) {
    resid <- (x - point)
    return(resid %*% t(resid))
  })
  between_cov <- Reduce("+", between_cov) / (m - 1)
  overall_cov <- within_cov + (1 + 1 / m) * between_cov
  out <- list(
    point = point,
    cov = overall_cov
  )
}
