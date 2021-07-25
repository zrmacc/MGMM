# Purpose: Functions to evaluate clustering quality and choose k. 
# Updated: 2021-07-24
# Note: Quality metrics must accommodate the case of empty clusters.

#------------------------------------------------------------------------------
# Quality Metrics
#------------------------------------------------------------------------------

#' Partition Data by Cluster Assignment
#'
#' @param data Numeric data.matrix.
#' @param assign Cluster assignments.
#' @return List of numeric data.matrices, separated by cluster assignment.
#' @noRd

PartitionByClust <- function(
  data,
  assign
) {
  labs <- sort(unique(assign))
  k <- length(labs)
  out <- lapply(seq_len(k), function(j){
    return(data[assign == labs[j], , drop = FALSE])
  })
  return(out)
}


#------------------------------------------------------------------------------

#' Within Cluster Dispersion
#' 
#' Used by \code{\link{CalHar}}.
#' 
#' @param clust Numeric data.matrix.
#' @param mean Cluster mean. 
#' @return d x d within-cluster dispersion matrix.
#' @noRd

WithinClusterDisp <- function(
  clust,
  mean
) {
  n <- nrow(clust)
  d <- ncol(clust)
  mean_mat <- matrix(data = mean, nrow = n, ncol = d, byrow = TRUE)
  resid <- clust - mean_mat
  out <- matIP(resid, resid)
  return(out)
}


#' Calinski-Harabaz Index
#'
#' Calculates the Calinski-Harabaz index.
#'
#' @param data Observations.
#' @param assign Assignments.
#' @param means List of cluster means.
#' @return Scalar metric.

CalHar <- function(data, assign, means) {
  
  # Split.
  split_data <- PartitionByClust(data, assign)
  
  # Dimensions.
  labs <- sort(unique(assign))
  d <- ncol(data)
  k <- length(labs)
  n <- nrow(data)
  
  # Total within cluster dispersion
  disp_within <- lapply(seq_len(k), function(j) {
    WithinClusterDisp(split_data[[j]], means[[j]])
  })
  disp_within <- Reduce("+", disp_within)
  
  # Grand mean.
  means <- do.call(rbind, means)
  grand_mean <- apply(means, 2, mean)
  
  # Total between cluster dispersion.
  mean_mat <- matrix(data = grand_mean, nrow = k, ncol = d, byrow = TRUE)
  resid <- means - mean_mat
  disp_between <- matIP(resid, resid)
  
  # CH statistic.
  out <- tr(disp_between) / tr(disp_within) * (n - k) / (k - 1)
  return(out)
}


#------------------------------------------------------------------------------

#' Mean Cluster Diameter
#' 
#' Used by \code{\link{DavBou}}.
#' 
#' @param clust Numeric matrix.
#' @param mean Cluster mean. 
#' @return Scalar mean diameter. 
#' @noRd

ClustDiam <- function(
  clust,
  mean
) {
  n <- nrow(clust)
  d <- ncol(clust)
  mean_mat <- matrix(data = mean, nrow = n, ncol = d, byrow = TRUE)
  resid <- clust - mean_mat
  rownames(resid) <- NULL
  diams <- plyr::aaply(
    .data = resid,
    .margins = 1,
    .fun = function(x) {sqrt(sum(x^2))}
  )
  mean_diam <- mean(diams)
  return(mean_diam)
}


#' Davies-Bouldin Index
#'
#' Calculates the Davies-Bouldin index.
#'
#' @param data Observations
#' @param assign Assignments
#' @param means List of cluster means
#' @return Scalar index.

DavBou <- function(
  data,
  assign,
  means
) {
  
  # Split.
  split_data <- PartitionByClust(data, assign)
  
  # Dimensions.
  labs <- sort(unique(assign))
  d <- ncol(data)
  k <- length(labs)
  n <- nrow(data)

  # Diameters.
  diams <- sapply(seq_len(k), function(j) {
    ClustDiam(split_data[[j]], means[[j]])
  })

  # Davies-Bouldin index.
  db_idx <- lapply(seq_len(k), function(j) {
    focus_mean <- means[[j]]
    focus_diam <- diams[[j]]
    
    # Loop over other clusters.
    scores <- c()
    for(l in 1:k){
      
      # Only consider cluster other than the focus. 
      if(l != j){
        mean_diff <- means[[l]] - focus_mean
        mean_sep <- sqrt(sum(mean_diff^2))
        score <- (diams[[l]] + focus_diam) / mean_sep
        scores <- c(scores, score)
      }
    }
    
    # Output max score.
    max_score <- max(scores)
    return(max_score)
  })
  db_idx <- mean(unlist(db_idx))
  
  # Output.
  return(db_idx) 
}


#------------------------------------------------------------------------------  

#' Cluster Quality
#'
#' Evaluates cluster quality. Returns the following metrics:
#' \itemize{
#' \item BIC: Bayesian Information Criterion, lower value indicates better clustering quality.
#' \item CHI: Calinski-Harabaz Index, higher value indicates better clustering quality.
#' \item DBI: Davies-Bouldin, lower value indicates better clustering quality.
#' \item SIL: Silhouette Width, higher value indicates better clustering quality.
#' }
#'
#' @param fit Object of class mix.
#' @return List containing the cluster quality metrics.
#'
#' @export
#' @seealso See \code{\link{ChooseK}} for using quality metrics to choose the cluster number.
#'
#' @examples
#' set.seed(100)
#' 
#' # Data generation
#' mean_list = list(
#' c(2, 2, 2),
#' c(-2, 2, 2),
#' c(2, -2, 2),
#' c(2, 2, -2)
#' )
#' 
#' data <- rGMM(n = 500, d = 3, k = 4, means = mean_list)
#' fit <- FitGMM(data, k = 4)
#' 
#' # Clustering quality
#' cluster_qual <- ClustQual(fit)

ClustQual <- function(fit) {
  
  # Unpack.
  data <- data.matrix(fit@Completed)
  assign <- fit@Assignments[, 1]
  means <- fit@Means
  k <- fit@Components
  d <- ncol(data)
  n <- nrow(data)
  
  # Exclude NA.
  data <- data[!is.na(assign), ]
  assign <- assign[!is.na(assign)]

  # Output structure.
  out <- list()

  ## BIC.
  out$BIC <- log(n) * (d * k) - fit@Objective

  ## Calinski-Harabaz Index.
  out$CHI <- CalHar(
    data,
    assign,
    means
  )

  ## Davies-Bouldin Index.
  out$DBI <- DavBou(
    data,
    assign,
    means
  )

  ## Silhouette.
  sil <- cluster::silhouette(
    x = assign, 
    dist = stats::dist(x = data, method = "euclidean")
  )
  out$SIL <- mean(sil[, 3])

  # Output.
  return(out)
}


#------------------------------------------------------------------------------

#' Attempt Model Fit and Return Quality Metrics. 
#' 
#' @param data Numeric data matrix.
#' @param k Number of clusters.
#' @param init_means Optional list of initial mean vectors.
#' @param fix_means Fix the means to their starting value? Must initialize.
#' @param init_covs Optional list of initial covariance matrices.
#' @param init_props Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @return Numeric vector containing the 4 cluster quality metrics. Returns
#'   null if the model fails to fit.
#' @noRd

ChooseKIter <- function(
  data, 
  k,
  init_means,
  fix_means, 
  init_covs, 
  init_props, 
  maxit, 
  eps
) {
  
  # Attempt fit. 
  fit <- tryCatch(
    expr = FitGMM(
      data, 
      k, 
      init_means, 
      fix_means, 
      init_covs,
      init_props, 
      maxit, 
      eps, 
      report = FALSE
    ),
    error = function(cond) {return(NULL)}
  )
  
  # Clustering quality metrics. 
  if (!is.null(fit)) {
    out <- unlist(ClustQual(fit))
  } else {
    out <- NULL
  }
  
  # Output.
  return(out)
}


#------------------------------------------------------------------------------

#' Summarize Results of Quality Metric Bootstrap. 
#' 
#' @param k Clusters.
#' @param boot_metrics Bootstrapped quality metrics. 
#' @param report Report bootstrap results? 
#' @return Either a data.frame reporting the means and standard errors
#'   of the quality metrics, or NULL if too few fits succeeded to calculate
#'  standard errors. 
#' @noRd

ChooseKSummarize <- function(
  k,
  boot_metrics,
  report
) {
  
  if (is.null(boot_metrics)) {
    
    # Report. 
    if (report) {
      cat("Model fails to fit observed data at cluster size: ", k, ".\n")
    }
    return(NULL)
  } else {
    
    nb <- nrow(boot_metrics)
    # Report. 
    if (report) {
      cat("Cluster size", k, "complete.", nb, "fit(s) succeeded.\n")
    }
    
    # Number of successful fits. 
    if(nb < 3){
      return(NULL)
    } else {
      
      # Summarize bootstrap results
      means <- apply(boot_metrics, 2, mean)
      ses <- apply(boot_metrics, 2, stats::sd)
      results <- data.frame(
        Clusters = k,
        Fits = nb, 
        Metric = names(means), 
        Mean = means, 
        SE = ses
      )
      rownames(results) <- seq_len(nrow(results))
      return(results)
    }
  }
}


#' Bootstrap Quality Metrics
#'
#' @param boot Number of bootstrap replicates.
#' @param data Numeric data matrix.
#' @param k Number of clusters.
#' @param boot Bootstrap replicates.
#' @param init_means Optional list of initial mean vectors.
#' @param fix_means Fix the means to their starting value? Must initialize.
#' @param init_covs Optional list of initial covariance matrices.
#' @param init_props Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @return Numeric matrix of clustering metrics. Returns null if the models
#'   fails to fit the observed data.
#' @noRd

ChooseKBootstrap <- function(
  boot,
  data, 
  k,
  init_means,
  fix_means, 
  init_covs, 
  init_props, 
  maxit, 
  eps
) {
  
  # Fitting procedure. 
  fit_wrapper <- function (data){
    out <- ChooseKIter(
      data, 
      k,
      init_means,
      fix_means, 
      init_covs, 
      init_props, 
      maxit, 
      eps
    )
    return(out)
  }
  
  # Quality metrics for observed data.
  qual_orig <- fit_wrapper(data)
  
  # If model fails to fit observed data at size k, exit. 
  if(is.null(qual_orig)){
    return(NULL)
  } else {
    
    # Else, fit bootstraps.
    n <- nrow(data)
    
    boot_results <- lapply(seq_len(boot), function(b) {
      draw <- sample(x = n, size = n, replace = TRUE)
      boot_data <- data[draw, ]
      qual_boot <- fit_wrapper(boot_data)
      return(qual_boot)
    })
    boot_results <- do.call(rbind, boot_results)
    out <- rbind(qual_orig, boot_results)
    return(out)
  }
}


#------------------------------------------------------------------------------

#' Recommend Cluster Number based on Bootstrap Results.
#' 
#' @param results Data.frame of bootstrap results.
#' @param metric String, metric of interest.
#' @param max_opt Is maximizing the metric optimal? 
#' @return Data.frame containing:
#' \itemize{
#'   \item `Metric` name.
#'   \item Optimal cluster number `k_opt` and metric value `Metric_opt`.
#'   \item 1SE cluster number `k_1se` and metric value `Metric_1se`.
#' }
#' @noRd

ChooseKRecommend <- function(
  results,
  metric,
  max_opt = FALSE
) {
  
  sub <- results[results$Metric == metric, ]
  
  if (max_opt) {
    sub$Mean <- (-1) * sub$Mean
  }
  
  # Optimal k.
  key <- which.min(sub$Mean)
  k_opt <- sub$Clusters[key]
  metric_opt <- sub$Mean[key]
  se_opt <- sub$SE[key]
  
  # 1SE k.
  sub_thresh <- sub[sub$Mean <= metric_opt + se_opt, ]
  key <- which.max(sub_thresh$Mean)
  k_1se <- sub_thresh$Clusters[key]
  metric_1se <- sub_thresh$Mean[key]
  
  if (max_opt) {
    metric_opt <- (-1) * metric_opt
    metric_1se <- (-1) * metric_1se
  }
  
  # Output.
  out <- data.frame(
    Metric = metric,
    k_opt = k_opt,
    Metric_opt = metric_opt,
    k_1se = k_1se,
    Metric_1se = metric_1se
  )
  return(out)
}


#------------------------------------------------------------------------------
# Cluster Number Selection
#------------------------------------------------------------------------------

#' Cluster Number Selection
#'
#' Function to choose the number of clusters k. Examines cluster numbers between
#' k0 and k1. For each cluster number, generates \code{boot} bootstrap data
#' sets, fits the Gaussian Mixture Model (\code{\link{FitGMM}}), and calculates
#' quality metrics (\code{\link{ClustQual}}). For each metric, determines the
#' optimal cluster number \code{k_opt}, and the \code{k_1SE}, the smallest
#' cluster number whose quality is within 1 SE of the optimum.
#'
#' @param data Numeric data matrix.
#' @param k0 Minimum number of clusters.
#' @param k1 Maximum number of clusters.
#' @param boot Bootstrap replicates.
#' @param init_means Optional list of initial mean vectors.
#' @param fix_means Fix the means to their starting value? Must provide initial
#'   values.
#' @param init_covs Optional list of initial covariance matrices.
#' @param init_props Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report bootstrap progress?
#' @return List containing \code{Choices}, the recommended number of clusters
#'   according to each quality metric, and \code{Results}, the mean and standard
#'   error of the quality metrics at each cluster number evaluated.
#'
#' @export
#' @seealso See \code{\link{ClustQual}} for evaluating cluster quality, and \code{\link{FitGMM}}
#' for estimating the GMM with a specified cluster number.
#'
#' @examples
#' \donttest{
#' set.seed(100)
#' mean_list <- list(c(2, 2), c(2, -2), c(-2, 2), c(-2, -2))
#' data <- rGMM(n = 500, d = 2, k = 4, means = mean_list)
#' choose_k <- ChooseK(data, k0 = 2, k1 = 6, boot = 10)
#' choose_k$Choices
#' }

ChooseK <- function(
  data, 
  k0 = 2, 
  k1 = NULL, 
  boot = 100,
  init_means = NULL,
  fix_means = FALSE, 
  init_covs = NULL, 
  init_props = NULL, 
  maxit = 10, 
  eps = 1e-4, 
  report = TRUE
) {
  
  # Check inputs.
  if (k0 < 2) {
    stop("At least 2 clusters are required to calculate quality metrics.")
  }
  
  if (is.null(k1)) {
    k1 <- k0 + 2
  }
  
  # Dimensions.
  n <- nrow(data)
  candidates <- seq(from = k0, to = k1)
  n_k <- length(candidates)
  
  # Loop over candidate cluster sizes.
  results <- lapply(candidates, function(k) {
    boot_metrics <- ChooseKBootstrap(
      boot,
      data, 
      k,
      init_means,
      fix_means, 
      init_covs, 
      init_props, 
      maxit, 
      eps
    )
    boot_results <- ChooseKSummarize(k, boot_metrics, report)
    return(boot_results)
    
  })
  results <- do.call(rbind, results)

  # Check for any results.
  if (is.null(results)) {
    
    cat("Unable to fit sufficient models at the current cluster sizes.")
    return(NULL)
    
  } else {
    
    choices <- list()
    choices[[1]] <- ChooseKRecommend(results, "BIC", max_opt = FALSE)
    choices[[2]] <- ChooseKRecommend(results, "CHI", max_opt = TRUE)
    choices[[3]] <- ChooseKRecommend(results, "DBI", max_opt = FALSE)
    choices[[4]] <- ChooseKRecommend(results, "SIL", max_opt = TRUE)
    choices <- do.call(rbind, choices)
    
    # Output.
    out <- list()
    out$Choices <- choices
    out$Results <- results
    return(out)
  }
}
