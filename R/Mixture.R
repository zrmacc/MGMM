# Purpose: Fits a multivariate normal mixture in the presence of missingness
# Updated: 19/01/24

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
#' @param fix.means Fix the means to their starting value? Must provide initial
#'   values.
#' @param S0 Optional list of initial covariance matrices.
#' @param pi0 Optional vector of initial cluster proportions. 
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective. 
#' @param report Report fitting progress?
#' @param parallel Run in parallel? Must register parallel backend first. 
#' @return Object of class \code{mix} containing the estimated 
#' 
#' @importFrom methods new
#' @importFrom mvnfast dmvn
#' @importFrom plyr aaply
#' @importFrom stats kmeans

fit.mix = function(Y,k=2,M0=NULL,fix.means=F,S0=NULL,pi0=NULL,maxit=100,eps=1e-6,report=F,parallel=F){
  # Dimensions
  d = ncol(Y);
  lab = colnames(Y);
  ID = seq(1:nrow(Y));
  
  ## Partition
  ind = complete.cases(Y);
  # Complete observations
  Y0 = Y[ind,];
  n0 = nrow(Y0);
  ID0 = ID[ind];
  # Check for incomplete obs
  Miss = sum(!ind);
  
  # Check for incomplete obs
  if(Miss==0){
    Y1 = ID1 = ID2 = NULL;
    n1 = n2 = 0;
    n = n0;
  # If missingness is present: 
  } else {
    # Incomplete obs
    Y1 = Y[!ind,];
    ID1 = ID[!ind];
    # Find completely missing obs
    keep = apply(Y1,1,FUN=function(x){sum(is.na(x))!=d});
    ID2 = ID1[!keep];
    # Separate partially and completely missing obs
    Y1 = Y1[keep,];
    ID1 = ID1[keep];
    n1 = nrow(Y1);
    n2 = length(ID2);
    n = n0+n1;
  }

  ## Initialization
  theta0 = list();
  
  # Case 1: All parameters provided
  if(!is.null(M0)&!is.null(S0)&!is.null(pi0)){
    theta0$M=M0;
    theta0$S=S0;
    theta0$pi=pi0;
  # Case 2: Initial values partially missing
  } else {
    if(n0==0){
      stop("If no observations are complete, initial values are required for all parameters.");
    }
    # If complete cases are available, apply kmeans:
    K = kmeans(x=Y0,centers=k,iter.max=100,nstart=100);
    # Assignments
    Z0 = K$cluster;
    # Initialize means
    if(is.null(M0)){
      M = K$centers;
      theta0$M = lapply(seq(1:k),function(i){M[i,]});
    } else {
      theta0$M = M0;
    }
    # Initialize covariances
    if(is.null(S0)){
      theta0$S = lapply(seq(1:k),function(i){cov(Y0[Z0==i,],Y0[Z0==i,])});
    } else {
      theta0$S = S0;
    }
    # Initialize proportions
    if(is.null(pi0)){
      theta0$pi = as.numeric(table(Z0))/n0;
    } else {
      theta0$pi = pi0;
    }
  } # End Case 2.

  # Check that estimated covariances are positive definite
  aux = function(j){
    E = eigen(x=j,symmetric=T,only.values=T);
    return(E$values);
  }
  evalues = unlist(lapply(X=theta0$S,FUN=aux));
  if(min(evalues)<=0){
    stop("Initial covariance matrices are not all positive definite.");
  };
  rm(aux);
  
  # Initial responsibilities
  theta0$G = Responsibility(M=theta0$M,S=theta0$S,pi=theta0$pi,Y0=Y0,Y1=Y1,parallel);
  
  ## Update Function
  Update = function(theta){
    # Parallelization
    if(!parallel){registerDoSEQ()};
    # Previous parameters
    M0 = theta$M;
    S0 = theta$S;
    p0 = theta$pi;
    
    ## Cluster sizes
    N0 = rep(0,k);
    if(n0>0){
      G0 = theta$G$G0;
      N0 = N0 + apply(G0,2,sum);
    }
    if(n1>0){
      G1 = theta$G$G1;
      N0 = N0 + apply(G1,2,sum);
    }
    
    ## Residual outer products
    aux = function(j){
      # Output structure
      Out = array(0,dim=c(d,d));
      # Complete cases
      if(n0>0){
        # Residuals
        m0 = matrix(data=M0[[j]],nrow=n0,ncol=d,byrow=T);
        E0 = Y0-m0;
        # Responsibility-weighted OP
        Out = Out + matIP(E0,G0[,j]*E0);
      };
      # Incomplete cases
      if(n1>0){
        # Responsibility-weighted OP
        Out = Out + ExpResidOP(Y1=Y1,m1=M0[[j]],m0=M0[[j]],S0=S0[[j]],g=G1[,j]);
      }
      return(Out);
    };
    V0 = lapply(seq(1:k),aux);
    
    ## Initial objective
    # Pi term
    P0 = sum(N0*log(p0));
    # Determinant term
    L = lapply(seq(1:k),function(j){N0[j]*log(det(S0[[j]]))});
    D0 = do.call(sum,L);
    # Trace term
    L = lapply(seq(1:k),function(j){tr(MMP(matInv(S0[[j]]),V0[[j]]))});
    T0 = do.call(sum,L);
    # Objective
    Q0 = P0-D0-T0;
    ## Update Means
    if(fix.means){
      M1 = M0;
    } else {
      # Update function
      aux = function(j){
        # Total
        Total = 0;
        # Complete cases
        if(n0>0){
          Total = Total + apply(G0[,j]*Y0,2,sum);
        };
        # Incomplete cases
        if(n1>0){
          Total = Total + SumWorkResp(Y1=Y1,m0=M0[[j]],S0=S0[[j]],g=G1[,j]);
        };
        # Update
        m1 = (Total)/N0[j];
        names(m1) = lab;
        return(m1);
      };
      # Update means
      M1 = lapply(seq(1:k),aux);
    };
    
    ## Update outer products
    aux = function(j){
      Out = array(0,dim=c(d,d));
      # Complete cases
      if(n0>0){
        # Residuals
        m1 = matrix(data=M1[[j]],nrow=n0,ncol=d,byrow=T);
        E0 = Y0-m1;
        # Responsibility-weighted OP
        Out = Out + matIP(E0,G0[,j]*E0);
      };
      # Incomplete cases
      if(n1>0){
        # Responsibility-weighted OP
        Out = Out + ExpResidOP(Y1=Y1,m1=M1[[j]],m0=M0[[j]],S0=S0[[j]],g=G1[,j]);
      };
      return(Out);
    }
    # Update outer products
    V1 = lapply(seq(1:k),aux);
    
    ## Update covariances
    aux = function(j){
      # Covariances
      S = V1[[j]]/N0[[j]];
      rownames(S) = colnames(S) = lab;
      return(S);
    }
    # Update covariances
    S1 = lapply(seq(1:k),aux);

    ## Update responsibilities
    G = Responsibility(M=M1,S=S1,pi=p0,Y0=Y0,Y1=Y1,parallel=parallel);
    
    # Cluster sizes
    N1 = rep(0,k);
    if(n0>0){
      G0 = theta$G$G0;
      N1 = N1 + apply(G0,2,sum);
    }
    if(n1>0){
      G1 = theta$G$G1;
      N1 = N1 + apply(G1,2,sum);
    }
    # Update cluster proportions
    p1 = N1/n;

    ## Final objective
    # Pi term
    P1 = sum(N1*log(p1));
    # Determinant term
    L = lapply(seq(1:k),function(j){N1[j]*log(det(S1[[j]]))});
    D1 = do.call(sum,L);
    # Trace term
    L = lapply(seq(1:k),function(j){tr(MMP(matInv(S1[[j]]),V1[[j]]))});
    T1 = do.call(sum,L);
    # Objective
    Q1 = P1-D1-T1;
    
    # Increment
    d = Q1-Q0;
    ## Output
    Out = list("M"=M1,"S"=S1,"pi"=p1,"G"=G,"Q1"=Q1,"Q0"=Q0,"d"=d);
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
  
  ## Cluster assignments
  
  # Responsibilities
  G = rbind(theta0$G$G0,theta0$G$G1);
  # Density evaluations
  D = rbind(theta0$G$D0,theta0$G$D1);
  
  # Assignments
  A = apply(G,1,which.max);
  A = cbind("ID"=c(ID0,ID1,ID2),"Assignment"=c(A,rep(NA,length(ID2))));
  A = A[order(A[,1]),];
  rownames(A) = seq(1:nrow(A));
  A = data.frame(A[,2,drop=F]);
  
  # Responsibilities
  G = cbind("ID"=c(ID0,ID1),G);
  G = rbind(G,cbind("ID"=c(ID2),array(NA,dim=c(length(ID2),k))));
  G = G[order(G[,1]),];
  rownames(G) = seq(1:nrow(G));
  G = G[,-1,drop=F];
  
  # Add entropies
  E = aaply(.data=G[],.margins=1,.fun=function(x){-sum(x*log(x))/log(k)});
  E[is.na(E)] = 0;
  A$Entropy = E;
  G = data.frame(G);
  
  # Density evaluations
  D = cbind("ID"=c(ID0,ID1),D);
  D = rbind(D,cbind("ID"=c(ID2),array(NA,dim=c(length(ID2),k))));
  D = D[order(D[,1]),];
  rownames(D) = seq(1:nrow(D));
  D = data.frame(D[,-1,drop=F]);
  
  ## Create completed data
  # Subjects with complete data
  Y_Complete = Y0;
  # Subjects with partial data
  if(n1>0){
    aux = function(j){return(WorkResp(Y1=Y1,m0=theta0$M[[j]],S0=theta0$S[[j]],g=theta0$G$G1[,j]))};
    L = lapply(seq(1:k),aux);
    Y1_Complete = Reduce("+",L);
    Y_Complete = rbind(Y_Complete,Y1_Complete);
  }
  # Subjects with no data
  if(n2>0){
    aux = function(j){return(theta0$M[[j]]*theta0$pi[j])};
    L = lapply(seq(1:k),aux);
    mu = Reduce("+",L);
    Y2_Complete = matrix(data=mu,nrow=n2,ncol=d,byrow=T);
    Y_Complete = rbind(Y_Complete,Y2_Complete);
  }
  # Add ID
  Y_Complete = cbind(c(ID0,ID1,ID2),Y_Complete);
  # Sort
  Y_Complete = Y_Complete[order(Y_Complete[,1]),-1];
  rownames(Y_Complete) = rownames(Y);
  
  # Objective
  Obj = theta0$Q1;
  names(Obj) = NULL;

  ## Output
  Out = new(Class="mix",Components=k,Means=theta0$M,Covariances=theta0$S,Proportions=theta0$pi,Objective=Obj,
            Density=D,Responsibilities=G,Assignments=A,Completed=Y_Complete);
  return(Out);
}

