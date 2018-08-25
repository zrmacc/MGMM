# Purpose: Fits a multivariate normal mixture in the presence of missingness
# Updated: 180823

########################
# Auxiliary Functions
########################

#' Responsibilities
#' 
#' Calculates the posterior probability of cluster membership given the observed
#' data.
#' 
#' @param M List of mean vectors.
#' @param S List of covariance matrices.
#' @param pi Vector of cluster proportions.
#' @param Y0 Complete observations.
#' @param Y1 Incomplete observations.
#' @param parallel Run in parallel? Must register parallel backend first.
#' @return List of numeric matrices, \code{G0} contains the responsibilities for
#'   complete observations, and \code{G1} for incomplete observations. Each 
#'   responsibility matrix is formatted with observations as rows, and mixture 
#'   components as columns.
#' @importFrom foreach foreach registerDoSEQ '%do%' '%dopar%' '%:%'
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
  i = j = NULL;
  
  ## Complete observations
  Out = list();
  # Loop over clusters
  if(n0>0){
    A0 = foreach(j=1:k,.combine=cbind) %do% {
      return(dmvn(X=Y0,mu=M[[j]],sigma=S[[j]])*pi[j]);
    }
    # Normalize
    G0 = aaply(.data=A0,.margins=1,.fun=function(x){x/sum(x)});
    # Format
    colnames(G0) = paste0("k",seq(1:k));
    rownames(G0) = seq(1:n0);
    Out$G0 = G0;
  }
  
  ## Incomplete observations
  # Loop over clusters %:% Loop over observations
  if(n1>0){
    A1 = foreach(j=1:k,.combine=cbind) %:% foreach(i=1:n1,.combine=rbind) %dopar% {
      # Current observation
      u = Y1[i,];
      key = !is.na(u);
      s = u[key];
      # Mean of observed components
      mij = M[[j]][key];
      # Covariance of observed components
      Sij = S[[j]][key,key];
      # Density contribution
      return(dmvn(X=s,mu=mij,sigma=Sij)*pi[j]);
    }
    # Normalize
    G1 = aaply(.data=A1,.margins=1,.fun=function(x){x/sum(x)});
    # Format
    colnames(G1) = paste0("k",seq(1:k));
    rownames(G1) = seq(1:n1);
    Out$G1 = G1;
  }
  return(Out);
}

#' Impute Sum of Observations with Missingness
#' 
#' Calculates the expected sum of response vectors, conditional on the
#' observed elements. 
#' 
#' @param Y1 Incomplete observations.
#' @param g Responsibilities.
#' @param m0 Previous mean.
#' @param S0 Previous covariance. 
#' @return Numeric vector, the sum of the individual conditional expectations. 

imputeSum = function(Y1,g,m0,S0){
  # Dimensions
  d = ncol(Y1);
  n1 = nrow(Y1);
  # Output structure
  Out = rep(0,d);
  # Loop over subjects
  i = NULL;
  for(i in 1:n1){
    # Current observation
    u = Y1[i,];
    key = is.na(u);
    # Permutation
    p = order(c(which(!key),which(key)));
    # Partition covariance
    Sts = S0[key,!key,drop=F];
    Sss = S0[!key,!key,drop=F];
    # Inverse covariance for observed elements
    Sssi = fastInv(Sss);
    # Observed components
    s = matrix(u[!key],ncol=1);
    zs0 = s-m0[!key];
    # Imputed observation
    t1 = m0[key]+as.numeric(fastMMp(Sts,fastMMp(Sssi,zs0)));
    yhat = c(s,t1);
    # Recover initial order
    yi = yhat[p];
    # Update
    Out = Out+g[i]*yi;
  }
  # Output
  return(Out);
}

#' Imputes Sum of Residual Outer Products
#' 
#' Calculates the expected sum of residual outer products,
#' conditional on the observed elements. 
#' 
#' @param Y1 Data for observations with missingness.
#' @param g Responsibilities.
#' @param m1 New mean.
#' @param m0 Old mean.
#' @param S0 Old covariance.
#' @return Numeric matrix, the sum of the individual conditional expectations. 

imputeOP = function(Y1,g,m1,m0,S0){
  # Dimensions
  d = ncol(Y1);
  n1 = nrow(Y1);
  # Output structure
  Out = array(0,dim=c(d,d));
  # Loop over subjects
  i = NULL;
  for(i in 1:n1){
    # Current observation
    u = Y1[i,];
    key = is.na(u);
    # Permutation
    p = order(c(which(!key),which(key)));
    ds = sum(!key);
    dt = sum(key);
    # Partition covariance
    Stt = S0[key,key,drop=F];
    Sts = S0[key,!key,drop=F];
    Sss = S0[!key,!key,drop=F];
    # Inverse covariance for observed elements
    Sssi = fastInv(Sss);
    # Conditional covariance of missing elements
    Ltti = SchurC(Stt,Sss,Sts);
    # Observed components
    s = matrix(u[!key],ncol=1);
    zs0 = s-m0[!key];
    # Imputed observation
    t1 = m0[key]+fastMMp(Sts,fastMMp(Sssi,zs0));
    yi = rbind(s,t1);
    # Residual
    e = yi-m1[p];
    # Residual outer product
    E = fastOP(e,e);
    Vi = E+rbind(array(0,dim=c(ds,d)),cbind(array(0,dim=c(dt,ds)),Ltti));
    # Recover initial order
    Vi = Vi[p,p];
    # Update
    Out = Out+g[i]*Vi;
  }
  # Output
  return(Out);
}

#' Update Function
#' 
#' EM parameter update for the normal mixture model with missingness.
#' 
#' @param theta Current parameter list.
#' @param Y0 Complete observations.
#' @param Y1 Incomplete observations.
#' @param parallel Run in parallel? Must register parallel backend first.
#' @return List containing the proposed means, covariances, cluster proportions,
#'   responsibilities, and objective function.

Update.mix = function(theta,Y0=NULL,Y1=NULL,parallel=F){
  # Parallelization
  if(!parallel){registerDoSEQ()};
  # Dimensions
  if(!is.null(Y0)){n0 = nrow(Y0); d=ncol(Y0); lab=colnames(Y0);} else {Y0 = NULL; n0 = 0};
  if(!is.null(Y1)){n1 = nrow(Y1); d=ncol(Y1); lab=colnames(Y1);} else {Y1 = NULL; n1 = 0};
  n = n0+n1;
  # Clusters
  k = length(theta$pi);
  j = NULL;
  # Cluster sizes
  N = rep(0,k);
  # Previous parameters
  M0 = theta$M;
  S0 = theta$S;
  
  ## Previous responsibilities
  if(n0>0){
    G0=theta$G$G0
    N = N + apply(G0,2,sum);
    };
  if(n1>0){
    G1=theta$G$G1
    N = N + apply(G1,2,sum);
    };
  
  ## Update Means
  # Loop over clusters
  M1 = foreach(j=1:k) %dopar% {
    # Contribution of complete cases
    if(n0>0){
      a0 = apply(G0[,j]*Y0,2,sum);
    } else {
      a0 = 0;
    }
    # Contribution of incomplete cases
    if(n1>0){
      a1 = imputeSum(Y1,G1[,j],M0[[j]],S0[[j]]);
    } else {
      a1 = 0;
    }
    # Update
    m1 = (a0+a1)/N[j];
    names(m1) = lab;
    return(m1);
  };
  
  ## Update covariances
  V0s = list();
  V1s = list();
  
  # Loop over clusters
  S1 = foreach(j=1:k) %do% {
    # Contribution of complete cases
    if(n0>0){
      V0 = wCov(Y=Y0,m=M1[[j]],w=G0[,j])*sum(G0[,j]);
    } else {
      V0 = array(0,c(d,d));
    }
    V0s[[j]] = V0;
    # Contribution of incomplete cases
    if(n1>0){
      V1 = imputeOP(Y1=Y1,g=G1[,j],m1=M1[[j]],m0=M0[[j]],S0=S0[[j]]);
    } else {
      V1 = array(0,c(d,d));
    }
    V1s[[j]] = V1;
    # Update
    S1 = (V0+V1)/N[j];
    rownames(S1) = colnames(S1) = lab;
    return(S1);
  }
  
  ## Update cluster proportions
  pi1 = N/n;
  
  ## Update responsibilities
  G2 = Responsibility(M=M1,S=S1,pi=pi1,Y0=Y0,Y1=Y1,parallel=parallel);
  # New objective
  q1 = Q.mix(S=S1,pi=pi1,Y0=Y0,Y1=Y1,G0=G0,G1=G1,V0s=V0s,V1s=V1s);
  
  ## Output
  Out = list("M"=M1,"S"=S1,"pi"=pi1,"G"=G2,"q"=q1);
  return(Out);
}

#' EM Objective Function
#' 
#' @param M Current means. Unnecessary if outer products are provided. 
#' @param S Current covariance matrices.
#' @param pi Current cluster proportions. 
#' @param Y0 Complete observations.
#' @param Y1 Incomplete observations.
#' @param G0 Responsibilities for complete observations.
#' @param G1 Responsibilities for incomplete observations.
#' @param V0s Outer products for complete observations.
#' @param V1s Outer products for incomplete observations. 
#' @return Scalar value of the objective. 
#' @importFrom foreach foreach '%do%' 

Q.mix = function(M,S,pi,Y0=NULL,Y1=NULL,G0,G1,V0s,V1s){
  # Dimensions
  if(!is.null(Y0)){n0 = nrow(Y0)} else {n0 = 0};
  if(!is.null(Y1)){n1 = nrow(Y1)} else {n1 = 0};
  # Clusters
  k = length(pi);
  j = NULL;
  # Cluster sizes
  N = rep(0,k);
  if(n0>0){
    N = N + apply(G0,2,sum);
    # Calculate outer products, if missing
    if(missing(V0s)){
      V0s = foreach(j=1:k) %do% {
        return(wCov(Y=Y0,m=M[[j]],w=G0[,j])*sum(G0[,j]));
      };
    };
  };
  if(n1>0){
    N = N + apply(G1,2,sum);
    # Calculate outer products, if missing
    if(missing(V1s)){
      V1s = foreach(j=1:k) %do% {
        return(imputeOP(Y1=Y1,g=G1[,j],m1=M[[j]],m0=M[[j]],S0=S[[j]]));
      };
    };
  };
  # Objective
  Out = 0;
  
  ## Objective
  # Loop over components
  for(j in 1:k){
    # Invert precision
    Lj = fastInv(S[[j]]);
    # Contribution of pi
    Out = Out+N[j]*log(pi[j]);
    # Contribution of determinant
    Out = Out+(N[j]/2)*log(fastDet(Lj));
    ## Contribution of outer products
    # Complete observations
    if(n0>0){
      Out = Out-0.5*tr(fastMMp(Lj,V0s[[j]]));
    }
    # Incomplete observations
    if(n1>0){
      Out = Out-0.5*tr(fastMMp(Lj,V1s[[j]]));
    }
  }
  # Output
  names(Out) = c("Q");
  return(Out);
}

########################
# Main Function
########################

#' Fit Multivariate Mixture Distribution
#' 
#' Given a matrix of random vectors, estimates the parameters for a mixture of
#' multivariate normal distributions. Accommodates arbitrary patterns of 
#' missingness, provided the elements are missing at random (MAR).
#' 
#' @param Y Numeric data matrix.
#' @param k Number of mixture components. Defaults to 2.
#' @param M0 Optional list of initial mean vectors.
#' @param S0 Optional list of initial covariance matrices.
#' @param pi0 Optional vector of initial cluster proportions. 
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective. 
#' @param report Report fitting progress?
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @return Object of class \code{mix} containing the estimated 
#' @importFrom methods new
#' @importFrom mvnfast dmvn
#' @importFrom plyr aaply
#' @importFrom stats cov kmeans

fit.mix = function(Y,k=2,M0=NULL,S0=NULL,pi0=NULL,maxit=100,eps=1e-6,report=F,parallel=F){
  # Dimensions
  d = ncol(Y);
  ID = seq(1:nrow(Y));
  # Partition
  ind = complete.cases(Y);
  # Complete observations
  Y0 = Y[ind,];
  n0 = nrow(Y0);
  ID0 = ID[ind];
  if(sum(!ind)>0){
    # Incomplete observations
    Y1 = Y[!ind,];
    ID1 = ID[!ind];
    # Exclude completely missing observations
    aux1 = function(x){sum(is.na(x))!=d};
    keep = apply(Y1,1,aux1);
    rm(aux1);
    Y1 = Y1[keep,];
    ID2 = ID1[!keep];
    ID1 = ID1[keep];
    n1 = nrow(Y1);
    n = n0+n1;
  } else {
    Y1 = ID1 = ID2 = NULL;
    n = n0;
  }
  
  ## Initialization
  theta0 = theta1 = list();
  if(n0>0){
    # Kmeans
    K = kmeans(x=Y0,centers=k,iter.max=10,nstart=10);
    # Assignments
    A = K$cluster;
    # Initial centers
    M = K$centers;
    theta0$M = lapply(seq(1:k),function(i){M[i,]});
    # Initial scatter matrices
    aux2 = function(i){
      # Subset
      Sub = Y0[A==i,];
      ni = nrow(Sub);
      if(ni>0){
        Si = cov(Sub);
        return(Si);
      };
    };
    theta0$S = lapply(seq(1:k),aux2);
    rm(aux2);
    # Initial cluster proportions
    theta0$pi = as.numeric(table(A))/n0;
    # Clean workspace
    rm(K,A,M);
  } else {
    if(is.null(M0)|is.null(S0)|is.null(pi0)){
      stop("If no observations are complete, initial values are required for all parameters.")
    };
  }
  
  # Overwrite if initial values are provided
  if(!is.null(M0)){theta0$M=M0};
  if(!is.null(S0)){theta0$S=S0};
  if(!is.null(pi0)){theta0$pi=pi0};
  
  # Initial reponsibilities
  theta0$G = Responsibility(M=theta0$M,S=theta0$S,pi=theta0$pi,Y0=Y0,Y1=Y1,parallel=parallel);
  # Initial objective
  theta0$q = Q.mix(M=theta0$M,S=theta0$S,pi=theta0$pi,Y0=Y0,Y1=Y1,G0=theta0$G$G0,G1=theta0$G$G1);
  
  ## EM Iterations
  i = NULL;
  for(i in 1:maxit){
    # Propose
    theta1 = Update.mix(theta0,Y0=Y0,Y1=Y1,parallel=parallel);
    # Check objective increment
    delta = theta1$q-theta0$q;
    # Accept increment if positive
    if(delta>0){
      theta0 = theta1;
      if(report){cat(paste0("Objective increment: ",signif(delta,3),"\n"))};
    };
    # Terminate if increment is insufficient
    if(delta<eps){
      break;
    };
  }; # End EM loop
  
  ## Fitting report
  if(report){
    if(i<maxit){
      cat(paste0(i-1," update(s) performed before tolerance limit."),"\n");
    } else {
      cat(paste0(i," update(s) performed without reaching tolerance limit."));
    }
  };
  
  ## Cluster assignments
  G = rbind(theta0$G$G0,theta0$G$G1);
  A = apply(G,1,which.max);
  A = cbind("ID"=c(ID0,ID1,ID2),"Assignment"=c(A,rep(NA,length(ID2))));
  A = A[order(A[,1]),];
  rownames(A) = seq(1:nrow(A));
  A = data.frame(A);
  # Responsibilities
  G = cbind("ID"=c(ID0,ID1),G);
  G = rbind(G,cbind("ID"=c(ID2),array(NA,dim=c(length(ID2),k))));
  G = G[order(G[,1]),];
  rownames(G) = seq(1:nrow(G));
  G = data.frame(G);
  ## Output
  Out = new(Class="mix",Components=k,Means=theta0$M,Covariances=theta0$S,Proportions=theta0$pi,
            Objective=theta0$q,Responsibilities=G,Assignments=A);
  return(Out);
}
