# Purpose: Master estimation function for MGMM
# Updated: 19/01/18

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
#' @param Y Numeric data matrix.
#' @param k Number of mixture components. Defaults to 1.
#' @param M0 Optional list of initial mean vectors.
#' @param fix.means Fix the means to their starting value? Must provide initial
#'   values.
#' @param S0 Optional list of initial covariance matrices.
#' @param pi0 Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @param parallel Run in parallel? Only used if \eqn{k>1}. Must register 
#'   parallel backend first.
#' @return For a single component model \eqn{k=1}, a list is returned,
#'   containing the estimated mean, covariance, and final EM objective. For a
#'   multi-component model \eqn{k>1}, an object of class \code{mix}, containing
#'   the estimated means, covariances, cluster proportions, cluster
#'   responsibilities, and observation assignments.
#' 
#' @export
#' @seealso For data generation, see \code{\link{rGMM}}.
#' 
#' @examples 
#' \dontrun{
#' # Single component without missingness
#' # Bivariate normal observations
#' Sigma = matrix(c(1,0.5,0.5,1),nrow=2);
#' Y = rMNMix(n=1e3,d=2,k=1,M=c(2,2),S=Sigma);
#' M = fit.GMM(Y=Y,k=1);
#' 
#' # Single component with missingness
#' # Trivariate normal observations
#' M = list(c(-2,-2,-2),c(2,2,2));
#' Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
#' Y = rMNMix(n=1e3,d=3,k=2,M=M,S=Sigma);
#' M = fit.GMM(Y=Y,k=2);
#' 
#' # Two components without missingness
#' # Trivariate normal observations
#' Means = list(c(-2,-2,-2),c(2,2,2));
#' Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
#' Y = rMNMix(n=1e3,d=3,k=2,M=Means,S=Sigma);
#' M = fit.GMM(Y=Y,k=2);
#' 
#' # Four components with missingness
#' # Bivariate normal observations
#' # Note: Fitting is slow. 
#' M = list(c(2,2),c(2,-2),c(-2,2),c(-2,-2));
#' S = 0.5*diag(2);
#' Y = rMNMix(n=1000,d=2,k=4,pi=c(0.35,0.15,0.15,0.35),m=0.1,M=M,S=S);
#' M = fit.GMM(Y=Y,k=4);
#' }

fit.GMM = function(Y,k=1,M0,fix.means=F,S0,pi0,maxit=100,eps=1e-6,report=T,parallel=F){
  ## Check data
  if(!is.matrix(Y)){stop("A numeric matrix with observations as rows is expected for Y.")};
  d = ncol(Y);
  ## Check initial values
  # Mean vectors
  if(!missing(M0)){
    # Object type
    if(!is.list(M0)){stop("If M0 is provided, a list of initial vectors, one for each component, is required.")};
    # Initial mean for each component
    if(length(M0)!=k){stop("If initial means are provided, one is required for each mixture component.")};
    # Dimensional consistency
    M0.d = unique(unlist(lapply(M0,length)));
    if((length(M0.d)>1)|(M0.d!=d)){stop("Each vector in M0 must have length of ncol(Y).")};
    rm(M0.d);
  } else {
    M0 = NULL;
  };
  if(fix.means&is.null(M0)){stop("If means are fixed, then initial values are required.")};
  # Covariance matrices
  if(!missing(S0)){
    # Object type 
    if(!is.list(S0)){stop("If S0 is provided, a list of initial matrices, one for each component, is required.")};
    # Initial covariance for each component
    if(length(S0)!=k){stop("If initial covariances are provided, one is required for each mixture component.")};
    # Dimensional consistency
    S0.d = unlist(unique(lapply(S0,dim)));
    if((length(S0.d)>2)|!all.equal(S0.d,c(d,d))){stop("Each matrix in S0 must have dimensions of ncol(Y) by ncol(Y).")};
    rm(S0.d);
  } else {
    S0 = NULL;
  };
  # Cluster proportions
  if(!missing(pi0)){
    # Object type
    if(!is.numeric(pi0)){stop("If pi0 is provided, a numeric vector of proportions is required.")};
    # Initial proportion for each component
    if(length(pi0)!=k){stop("If initial proportions are provided, one is required for each mixture component.")};
  } else {
    pi0 = NULL;
  };
  
  ## Case 1: Single mixture component
  if(k==1){
    Out = fit.mvn(Y=Y,m0=M0[[1]],fix.means=fix.means,
                  S0=S0[[1]],maxit=maxit,eps=eps,report=report);
    return(Out);
  } else {
  ## Case 2: Multiple mixture components  
    Out = fit.mix(Y=Y,k=k,M0=M0,fix.means=fix.means,
                  S0=S0,pi0=pi0,maxit=maxit,eps=eps,report=report,parallel=parallel);
    return(Out);
  };
}