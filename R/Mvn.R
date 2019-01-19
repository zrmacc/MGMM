# Purpose: Fitting function for fitting a single multivariate normal with missingness.
# Updated: 19/01/18

########################
# Main Function
########################

#' Fit Multivariate Normal Distribution
#'
#' Given a matrix of random vectors, estimates the mean and covariance of a
#' multivariate normal distribution. Accommodates arbitrary patterns of
#' missingness at random.
#' 
#' @param Y Numeric data matrix.
#' @param m0 Optional initial mean vector.
#' @param fix.means Fix the means to their starting value? Must provide initial
#'   values.
#' @param S0 Optional initial covariance matrix.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report fitting progress?
#' @return List containing the estimated mean, covariance, and objective. 
#' @importFrom stats complete.cases

fit.mvn = function(Y,m0=NULL,fix.means=F,S0=NULL,maxit=100,eps=1e-6,report=T){
  # Check for missingness
  Miss = sum(is.na(Y))>0;
  lab = colnames(Y);
  
  ## If missingness is absent:
  if(!Miss){
    ## Estimation
    if(fix.means){
      m1 = m0;
    } else {
      # Mean
      m1 = apply(Y,2,mean);
      names(m1) = lab;
    };
    # Covariance
    S1 = cov(Y,Y);
    rownames(S1) = colnames(S1) = lab;
    L1 = matInv(S1);
    # Objective
    n = nrow(Y);
    V = matIP(Y-m1,Y-m1);
    Q = -n*log(det(S1))-tr(MMP(L1,V));
    
    ## Output
    Out = list("Mean"=m1,"Covariance"=S1,"Objective"=Q);
    return(Out);
  } else {
    ## If missingness is present:
    # Dimensions
    d = ncol(Y);
    ID = seq(1:nrow(Y));
    
    # Partition
    ind = complete.cases(Y);
    # Complete observations
    Y0 = Y[ind,];
    ID0 = ID[ind];
    n0 = nrow(Y0);
    if(n0>0){
      # Sum of complete observations
      t0 = apply(Y0,2,sum);
      # Incomplete observations
      Y1 = Y[!ind,];
      ID1 = ID[!ind];
      # Exclude completely missing observations
      aux = function(x){sum(is.na(x))!=d};
      keep = apply(Y1,1,aux);
      ID2 = ID1[!keep];
      Y1 = Y1[keep,];
      ID1 = ID1[keep];
      n1 = nrow(Y1);
      n2 = length(ID2);
      n = n0+n1;
    } else {
      if(is.null(m0)|is.null(S0)){
        stop("If no observations are complete, initial values are required for all parameters.")
      };
    };
    
    ## Initialization
    theta0 = theta1 = list();
    # Initial values
    if(n0>0){
      theta0$m = apply(Y0,2,mean);
      theta0$S = cov(Y0,Y0);
    } else {
      if(is.null(m0)|is.null(S0)){
        stop("If no observations are complete, initial values are required for all parameters.")
      };
    };
    # Overwrite if initial values were provided
    if(!is.null(m0)){theta0$m=m0};
    if(!is.null(S0)){theta0$S=S0};
    
    ## Define update function
    Update = function(theta){
      # Previous parameters
      m0 = theta$m;
      S0 = theta$S;
      L0 = matInv(S0);
      
      ## Initial objective
      V0 = matIP(Y0-m0,Y0-m0);
      V1 = ExpResidOP(Y1=Y1,m1=m0,m0=m0,S0=S0);
      Q0 = (-n)*log(det(S0))-tr(MMP(L0,V0+V1));
      
      ## Update mean
      if(fix.means){
        m1 = m0;
      } else {
        m1 = 0;
        # Complete observations
        if(n0>0){
          m1 = m1+t0;
        };
        # Incomplete observations
        if(n1>0){
          m1 = m1+SumWorkResp(Y1=Y1,m0=m0,S0=S0);
        }
        # Normalize
        m1 = (m1/n);
      };
      
      ## Update covariance
      S1 = array(0,c(d,d));
      
      # Complete observations
      if(n0>0){
        V0 = matIP(Y0-m1,Y0-m1);
        S1 = S1+V0;
      }
      # Incomplete observations
      if(n1>0){
        V1 = ExpResidOP(Y1=Y1,m1=m1,m0=m0,S0=S0);
        S1 = S1+V1;
      }
      # Normalize
      S1 = (S1/n);
      rownames(S1) = colnames(S1) = lab;
      
      ## Final objective
      L1 = matInv(S1);
      Q1 = (-n)*log(det(S1))-tr(MMP(L1,V0+V1));
      
      # Increment
      d = Q1-Q0;
      
      ## Output
      Out = list("m"=m1,"S"=S1,"Q1"=Q1,"Q0"=Q0,"d"=d);
      return(Out);
    }
    
    ## Maximzation
    i = 1;
    for(i in 1:maxit){
      # Update
      theta1 = Update(theta0);
      # Accept if increment is positive
      if(theta1$d>0){
        theta0 = theta1;
        if(report){cat("Objective increment: ",signif(theta1$d,digits=3),"\n")}
      }
      # Terminate if increment is below tolerance
      if(theta1$d<eps){
        # If EM failes to perform any updates, keep initial objective
        if(i==1){theta0$Q1=theta1$Q0};
        rm(theta1);
        break;
      };
    };
    
    ## Fitting report
    if(report){
      if(i<maxit){
        cat(paste0(i-1," update(s) performed before tolerance limit.\n\n"));
      } else {
        cat(paste0(i," update(s) performed without reaching tolerance limit.\n\n"));
      };
    };
    
    ## Completed data
    Y_Complete = Y0;
    Y1_Complete = WorkResp(Y1=Y1,m0=theta0$m,S0=theta0$S);
    Y_Complete = rbind(Y_Complete,Y1_Complete);
    if(n2>0){
      Y2_Complete = foreach(i=1:length(ID2),.combine=rbind) %do% {
        return(theta0$m);
      };
      Y_Complete = rbind(Y_Complete,Y2_Complete);
    }
    # Order
    Y_Complete = cbind(c(ID0,ID1,ID2),Y_Complete);
    Y_Complete = Y_Complete[order(Y_Complete[,1]),-1];
    rownames(Y_Complete) = rownames(Y);
    
    ## Output
    Out = list("Mean"=theta0$m,"Covariance"=theta0$S,"Objective"=theta0$Q1,"Completed"=Y_Complete);
    return(Out);
  };
};
