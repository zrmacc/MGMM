# Purpose: Calculate E-step expectations of missing observations
# Updated: 19/01/17

#' Responsibilities
#' 
#' Calculates the posterior probability of cluster membership given the observed
#' data.
#' 
#' @param M List of mean vectors.
#' @param S List of covariances matrices.
#' @param pi Vector of cluster proportions.
#' @param Y0 Matrix of complete observations.
#' @param Y1 Matrix of incomplete observations.
#' @param parallel Run in parallel? Must register parallel backend first.
#' @return List of numeric matrices, \code{G0} contains the responsibilities for
#'   complete observations, and \code{G1} for incomplete observations. Each 
#'   responsibility matrix is formatted with observations as rows, and mixture 
#'   components as columns.
#'   
#' @importFrom foreach foreach registerDoSEQ '%do%' '%dopar%'
#' @importFrom plyr aaply

Responsibility = function(M,S,pi,Y0=NULL,Y1=NULL,parallel=F){
  # Parallelization
  if(!parallel){registerDoSEQ()};
  
  # Dimensions
  if(!is.null(Y0)){n0 = nrow(Y0)} else {n0 = 0};
  if(!is.null(Y1)){n1 = nrow(Y1)} else {n1 = 0};
  n = n0+n1;
  
  # Clusters
  k = length(pi);
  i = j = 1;
  
  # Output structure
  Out = list();
  
  ## Complete observations
  if(n0>0){
    A0 = foreach(j=1:k,.combine=cbind) %do% {
      return(dmvn(X=Y0,mu=M[[j]],sigma=S[[j]])*pi[j]);
    }
    # Normalize
    G0 = aaply(.data=A0,.margins=1,.fun=function(x){x/sum(x)},.drop=F);
    # Format
    colnames(G0) = paste0("k",seq(1:k));
    rownames(G0) = seq(1:n0);
    Out$G0 = G0;
  }
  
  ## Incomplete observations
  if(n1>0){
    A1 = foreach(i=1:n1,.combine=rbind) %dopar% {
      # Current observation
      y = Y1[i,];
      key = !is.na(y);
      
      # Observed components
      s = y[key];
      
      # Loop over clusters
      Sub = foreach(j=1:k,.combine=c) %do% {
        # Mean of observed components
        mij = M[[j]][key];
        # Covariance of observed components
        Sij = S[[j]][key,key];
        # Density contribution
        return(dmvn(X=s,mu=mij,sigma=Sij)*pi[j]);
      };
      
      # Output
      return(Sub);
    }
    
    # Normalize
    G1 = aaply(.data=A1,.margins=1,.fun=function(x){x/sum(x)},.drop=F);
    # Format
    colnames(G1) = paste0("k",seq(1:k));
    rownames(G1) = seq(1:n1);
    Out$G1 = G1;
  }
  return(Out);
}

#' Sum of Working Response Vectors
#' 
#' Calculate the sum of the working response vectors. 
#'
#' @param Y1 Incomplete observations.
#' @param m0 Previous mean.
#' @param S0 Previous covariance.
#' @param g Responsibilities
#' @return Numeric vector, the responsibility-weighted cumulative working
#'   response vector.
#' 
#' @importFrom foreach foreach '%do%'

SumWorkResp = function(Y1,m0,S0,g=NULL){
  # Dimensions
  d = ncol(Y1);
  n1 = nrow(Y1);
  
  # Reponsibilities
  if(is.null(g)){
    g = rep(1,n1);
  }
  
  # Loop over observations
  i = 1;
  Out = foreach(i=1:n1,.combine="+") %do% {
    # Current observation
    y = Y1[i,];
    key = is.na(y);
    
    # Permutation
    Perm = c(which(!key),which(key));
    Rev_Perm = order(Perm);
    
    # Partition covariance
    Cov_TS = S0[key,!key,drop=F];
    Var_SS = S0[!key,!key,drop=F];
    Var_SS_Inv = matInv(Var_SS);
    
    # Observed components
    s = matrix(y[!key],ncol=1);
    
    # Conditional expectation of missing compoments
    t = m0[key] + MMP(Cov_TS,MMP(Var_SS_Inv,s-m0[!key]));
    
    # Working response
    Yhat = rbind(s,t);
    
    # Return
    return(g[i]*Yhat[Rev_Perm])
  }
  return(Out);
}

#' Expected Residual Outer Product
#' 
#' Calculates the expected residual outer product. 
#' 
#' @param Y1 Data for observations with missingness.
#' @param m1 New mean.
#' @param m0 Old mean.
#' @param S0 Old covariance.
#' @param g Responsibilities.
#' 
#' @return Numeric matrix, the responsibility-weighted, cumulative,
#'  expected residual outer product. 
#'  
#' @importFrom foreach foreach '%do%'

ExpResidOP = function(Y1,m1,m0,S0,g=NULL){
  # Dimensions
  d = ncol(Y1);
  n1 = nrow(Y1);
  
  # Reponsibilities
  if(is.null(g)){
    g = rep(1,n1);
  }
  
  # Output structure
  
  # Loop over observations
  i = 1;
  Out = foreach(i=1:n1,.combine="+") %do% {
    # Current observation
    y = Y1[i,];
    key = is.na(y);
    
    # Permutation
    Perm = c(which(!key),which(key));
    Rev_Perm = order(Perm);
    ds = sum(!key);
    
    # Partition covariance
    Var_TT = S0[key,key,drop=F];
    Cov_TS = S0[key,!key,drop=F];
    Var_SS = S0[!key,!key,drop=F];
    
    # Inverses
    Var_SS_Inv = matInv(Var_SS);
    Lambda_TT_Inv = SchurC(Var_TT,Var_SS,Cov_TS);
    
    # Observed components
    s = matrix(y[!key],ncol=1);
    
    # Conditional expectation of missing compoments
    t = m0[key] + MMP(Cov_TS,MMP(Var_SS_Inv,s-m0[!key]));
    
    # Working response
    Yhat = rbind(s,t);
    
    # Residual
    Resid = Yhat-m1[Perm];
    
    # Residual outer product
    E = matOP(Resid,Resid);
    
    # Add correction
    idx = seq(from=ds+1,to=d);
    E[idx,idx] = E[idx,idx] + Lambda_TT_Inv;
    
    # Recover initial order
    E = E[Rev_Perm,Rev_Perm];
    
    # Return
    return(g[i]*E);
  }
  # Output
  return(Out);
}
