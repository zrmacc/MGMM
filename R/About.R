# Purpose: Package documentation
# Updated: 19/01/17

#' @useDynLib MGMM
#' @importFrom Rcpp sourceCpp
NULL

#' MGMM: Missingness Aware Gaussian Mixture Models
#' 
#' Parameter estimation and classification via Gaussian Mixture Models, allowing
#' for missingness in the input vectors. See \code{\link{fit.GMM}} for estimating
#' the GMM, and \code{\link{chooseK}} for selecting the number of cluseters. See 
#' \code{\link{rGMM}} for simulating data from a GMM.
#' 
#' @author Zachary R. McCaw
#' @docType package
#' @name MGMM
NULL