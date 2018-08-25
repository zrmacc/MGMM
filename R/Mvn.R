# Purpose: Fitting function for a single component 'mixture'
# Updated: 180824

########################
# Auxiliary Functions
########################

#' EM Objective Function
#' 
#' Expectation of the complete data log likelihood, given the
#' observed data, for a multivariate normal distribution with
#' a single component. 
#' 
#' @param m Current mean, unnecessary if outer products are provided.
#' @param S Current covariance matrix.
#' @param Y0 Complete observations.
#' @param Y1 Incomplete observations.
#' @param V0 Outer products for complete observations.
#' @param V1 Outer products for incomplete observations. 
#' @return Scalar value of the objective. 

Q.mvn = function(m,S,Y0=NULL,Y1=NULL,V0,V1){
  # Dimensions
  if(!is.null(Y0)){n0 = nrow(Y0)} else {n0 = 0};
  if(!is.null(Y1)){n1 = nrow(Y1)} else {n1 = 0};
  n = n0+n1;
  
  ## Calculate outer products, if missing
  # Complete observations
  if((n0>0)&missing(V0)){
    M0 = matrix(rep(x=m,times=n0),nrow=n0,byrow=T);
    E0 = Y0-M0;
    V0 = fastIP(E0,E0);
  };
  # Incomplete observations
  if((n1>0)&missing(V1)){
    V1 = imputeOP(Y1=Y1,g=rep(1,n1),m1=m,m0=m,S0=S);
  }
  
  ## Initialize objective
  Out = 0;
  # Precision
  L = fastInv(S);
  
  ## Determinant Term
  Out = Out+(n/2)*log(fastDet(L));
  
  ## Outer Product Term
  # Complete cases
  if(n0>0){
    Out = Out-(0.5)*tr(fastMMp(L,V0));
  }
  # Incomplete cases
  if(n1>0){
    Out = Out-(0.5)*tr(fastMMp(L,V1));
  }
  return(Out);
}

#' Update Function
#' 
#' EM parameter update for a multivariate normal distribution
#' with a single component.
#' 
#' @param theta Current parameter list, containing the mean, covariance, and
#'   objective.
#' @param Y0 Complete observations.
#' @param Y1 Incomplete observations.
#' @return List containing the proposed mean, covariance, and objective. 

Update.mvn = function(theta,Y0=NULL,Y1=NULL){
  # Dimensions
  if(!is.null(Y0)){n0 = nrow(Y0); d=ncol(Y0); lab=colnames(Y0);} else {Y0 = array(dim=c(0,0)); n0 = 0};
  if(!is.null(Y1)){n1 = nrow(Y1); d=ncol(Y1); lab=colnames(Y1);} else {Y1 = array(dim=c(0,0)); n1 = 0};
  n = n0+n1;
  # Previous parameters
  m0 = theta$m;
  S0 = theta$S;
  
  ## Update Mean
  m1 = 0;
  # Complete observations
  if(n0>0){
    m1 = m1+apply(Y0,2,sum);
  };
  # Incomplete observations
  if(n1>0){
    m1 = m1+imputeSum(Y1=Y1,g=rep(1,n1),m0=m0,S0=S0);
  }
  # Normalize
  m1 = (m1/n);
  names(m1) = lab;
  
  ## Update Covariance
  S1 = array(0,c(d,d));
  # Complete observations
  if(n0>0){
    M0 = matrix(rep(x=m1,times=n0),nrow=n0,byrow=T);
    E0 = Y0-M0;
    V0 = fastIP(E0,E0);
    S1 = S1+V0;
  }
  # Incomplete observations
  if(n1>0){
    V1 = imputeOP(Y1=Y1,g=rep(1,n1),m1=m1,m0=m0,S0=S0);
    S1 = S1+V1;
  }
  # Normalize
  S1 = (S1/n);
  rownames(S1) = colnames(S1) = lab;
  
  ## Update Objective
  q1 = Q.mvn(S=S1,Y0=Y0,Y1=Y1,V0=V0,V1=V1);
  
  ## Output
  Out = list("m"=m1,"S"=S1,"q"=q1);
  return(Out);
}

########################
# Main Function
########################

#' Fit Multivariate Normal Distribution
#' 
#' Given a matrix of random vectors, estimates the parameters of a multivariate 
#' normal distribution Accommodates arbitrary patterns of missingness, provided
#' the elements are missing at random (MAR).
#' 
#' @param Y Numeric data matrix.
#' @param m0 Optional initial mean vector.
#' @param S0 Optional initial covariance matrix.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return List containing the estimated mean and covariance. If missingness is 
#'   present, the final value of the EM objective is also provided.
#' @importFrom stats complete.cases

fit.mvn = function(Y,m0=NULL,S0=NULL,maxit=100,eps=1e-6,report=T){
  # Check for missingness
  Miss = sum(is.na(Y))>0;
  lab = colnames(Y);
  ## If missingness is absent
  if(!Miss){
    ## Estimation
    # Mean
    m1 = apply(Y,2,mean);
    names(m1) = lab;
    # Covariance
    S1 = cov(Y);
    rownames(S1) = colnames(S1) = lab;
    
    ## Output
    Out = list("Mean"=m1,"Covariance"=S1);
    return(Out);
    
  } else {
  ## If missingness is present
    # Dimensions
    d = ncol(Y);
    # Partition
    ind = complete.cases(Y);
    # Complete observations
    Y0 = Y[ind,];
    n0 = nrow(Y0);
    # Incomplete observations
    Y1 = Y[!ind,];
    # Exclude completely missing observations
    aux = function(x){sum(is.na(x))!=d};
    keep = apply(Y1,1,aux);
    Y1 = Y1[keep,];
    n1 = nrow(Y1);
    n = n0+n1;
    
    ## Initialization
    theta0 = theta1 = list();
    # Initial values
    if(n0>0){
      theta0$m = apply(Y0,2,mean);
      theta0$S = cov(Y0);
    } else {
      if(is.null(m0)|is.null(S0)){
        stop("If no observations are complete, initial values are required for all parameters.")
      };
    };
    # Overwrite if initial values were provided
    if(!is.null(m0)){theta0$m=m0};
    if(!is.null(S0)){theta0$S=S0};
    
    # Initial objective
    theta0$q = Q.mvn(m=theta0$m,S=theta0$S,Y0=Y0,Y1=Y1);
    
    ## EM Iterations
    i = NULL;
    for(i in 1:maxit){
      # Propose
      theta1 = Update.mvn(theta=theta0,Y0=Y0,Y1=Y1);
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
    
    ## Output
    Out = list("Mean"=theta0$m,"Covariance"=theta0$S,"Objective"=theta0$q);
    return(Out);
  }
}
