# Purpose: Data generation for mixture of normals
# Updated: 180821

#' Data Generation from Multivariate Normal Mixture Models
#' 
#' Generates an \eqn{n\times d} matrix of multivariate normal random vectors
#' with observations as rows. If \eqn{k=1}, all observations belong to the same
#' cluster. If \eqn{k>1} the observations are generated via two-step precoedure.
#' First, the cluster membership is drawn from a multinomial distribution, with
#' mixture proportions specified by \code{pi}. Conditional on cluster
#' membership, the observation is drawn from a multivariate normal distribution,
#' with cluster-specific mean and covariance. The cluster means are provided
#' using \code{M}, and the cluster covariance matrices are provided using
#' \code{S}. If \eqn{m>0}, missingness is introduced, completely at random, by
#' setting that proportion of elements in the data matrix to \code{NA}.
#' 
#' @param n Observations.
#' @param d Observation dimension.
#' @param k Number of mixture components. Defaults to 1.
#' @param pi Mixture proportions. If omitted, components are assumed 
#'   equi-probable.
#' @param m Proportion of elements missing, \eqn{m\in[0,1)}.
#' @param M Either a prototype mean vector, or a list of mean vectors. Defaults
#'   to the zero vector.
#' @param S Either a prototype covariance matrix, or a list of covariance 
#'   matrices. Defaults to the identity matrix.
#' @return Numeric matrix with observations as rows. 
#' @importFrom foreach foreach '%do%'
#' @importFrom mvnfast rmvn
#' @importFrom stats rmultinom
#' @export
#' 
#' @examples 
#' # Single component without missingness
#' # Bivariate normal observations
#' Sigma = matrix(c(1,0.5,0.5,1),nrow=2);
#' Y = rMNMix(n=1e3,d=2,k=1,M=c(2,2),S=Sigma);
#' 
#' # Single component with missingness
#' # Trivariate normal observations
#' M = list(c(-2,-2,-2),c(2,2,2));
#' Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
#' Y = rMNMix(n=1e3,d=3,k=2,M=M,S=Sigma);
#' 
#' # Two components without missingness
#' # Trivariate normal observations
#' Means = list(c(-2,-2,-2),c(2,2,2));
#' Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
#' Y = rMNMix(n=1e3,d=3,k=2,M=Means,S=Sigma);
#' 
#' # Four components with missingness
#' # Bivariate normal observations
#' M = list(c(2,2),c(2,-2),c(-2,2),c(-2,-2));
#' S = 0.5*diag(2);
#' Y = rMNMix(n=1000,d=2,k=4,pi=c(0.35,0.15,0.15,0.35),m=0.1,M=M,S=S);

rMNMix = function(n,d=2,k=1,pi,m=0,M,S){
  ## Input checks
  # Mixture proportions
  if(missing(pi)){pi = rep(1,k)/k} else {pi = pi/sum(pi)};
  # Means
  if(missing(M)){M=rep(0,d)};
  if(!is.list(M)){
    M = list(M);
    M = rep(M,k);
  }
  # Covariance matrices
  if(missing(S)){S = diag(d)};
  if(!is.list(S)){
    S = list(S);
    S = rep(S,k);
  }
  # Check dimensional consistency
  if(k!=length(pi)){stop("pi must have length k.")};
  d.M = unlist(unique(lapply(M,length)));
  if(d.M!=d){stop("Each vector in M must have length d.")};
  d.S = unlist(unique(lapply(S,dim)));
  if(!all.equal(d.S,rep(d,2))){stop("Each matrix in S must have dimension d x d.")};
  # Check input validity
  if(m<0|m>=1){stop("Missingness m must reside in [0,1).")};
  if((length(pi)>1)&(min(pi)<=0|max(pi)>=1)){stop("Elements of pi must reside in (0,1).")};
  
  ## Data generation
  if(k==1){
    # Case 1: Single mixture component
    Mi = M[[1]];
    Si = S[[1]];
    Y = rmvn(n=n,mu=Mi,sigma=Si);
  } else {
    # Case 2: Multiple mixture components
    z = rmultinom(n=n,size=1,prob=pi);
    aux = function(x){which(x==1)};
    z = apply(z,2,aux);
    i = NULL;
    Y = foreach(i=1:k,.combine=rbind) %do% {
      ni = sum(z==i);
      # Output only if component is non-empty
      if(ni>0){
        Mi = M[[i]];
        Si = S[[i]];
        return(rmvn(n=ni,mu=Mi,sigma=Si));
      }
    }
    # Permute row order
    Y = Y[sample(nrow(Y),replace=F),,drop=F];
  }
  
  ## Introduce missingness
  if(m>0){
    # Missing elements
    e = length(Y)
    nm = floor(e*m);
    draw = sort(sample(x=e,size=nm,replace=F));
    Y[draw] = NA;
  }
  
  ## Format
  colnames(Y) = paste0("y",seq(1:d));
  rownames(Y) = seq(1:n);
  return(Y);
}