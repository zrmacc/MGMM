# Purpose: Package documentation.
# Updated: 2024-07-02

#' Missingness-Aware Gaussian Mixture Models
#'
#' Estimate and classify using Gaussian mixture models (GMMs) when the input
#' contains missing values. Uses an expectation–conditional maximization (ECM)
#' algorithm with support for full component covariances and optional ridge
#' regularization. See McCaw et al. (2022) <doi:10.1186/s12859-022-04740-9>.
#'
#' @name MGMM-package
#' @aliases MGMM MGMM-package _PACKAGE
#' @useDynLib MGMM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"


## usethis namespace: start
## usethis namespace: end
NULL