# Purpose: Functions to evaluate clustering quality and choose k
# Updated: 19/01/24
# Note: Quality metrics must accommodate the case of empty clusters. 

########################
# Quality Metrics
########################

#' Calinski-Harabaz Index
#' 
#' Calculates the Calinski-Harabaz index.
#' 
#' @param Y Observations
#' @param a Assignments
#' @param m List of cluster means
#' @return Scalar.

CalHar = function(Y,a,m){
  # Clusters
  labs = sort(unique(a));
  k = length(labs);
  n = nrow(Y);
  
  ## Total within cluster dispersion
  aux = function(j){
    # Subset obserations in cluster
    key = labs[j];
    Sub = Y[a==key,,drop=F];
    nj = nrow(Sub);
    mj = m[[key]];
    # Loop over obs in cluster
    aux1 = function(i){return(matOP(Sub[i,]-mj,Sub[i,]-mj))};
    L = lapply(seq(1:nj),aux1);
    Wj = Reduce("+",L);
    return(Wj);  
  };
  # Within cluster dispersion
  L = lapply(seq(1:k),aux);
  W = Reduce("+",L);
  
  ## Grand mean
  aux = function(j){
    key = labs[j];
    nj = sum(a==key);
    if(nj>0){return(m[[key]])};
  }
  # Grand total
  L = lapply(seq(1:k),aux);
  Total = Reduce("+",L);
  # Grand mean
  mu = Total/k;
  
  ## Total between cluster disperson
  aux = function(j){
    # Subset obserations in cluster
    key = labs[j];
    nj = sum(a==key);
    mj = m[[key]];
    e = (mj-mu);
    return(matOP(e,e));
  };
  # Between cluster dispersion
  L = lapply(seq(1:k),aux);
  B = Reduce("+",L);
  # CH statistic
  Out = tr(B)/tr(W)*(n-k)/(k-1);
  return(Out);
}

#' Davies-Bouldin Index
#' 
#' Calculates the Davies-Bouldin index.
#' 
#' @param Y Observations
#' @param a Assignments
#' @param m List of cluster means
#' @return Scalar.

DavBou = function(Y,a,m){
  # Clusters
  labs = sort(unique(a));
  k = length(labs);
  
  ## Cluster Diameters
  aux = function(j){
    # Subset obserations in cluster
    key = labs[j];
    Sub = Y[a==key,,drop=F];
    nj = nrow(Sub);
    mj = m[[key]];
    # Loop over obs in cluster
    aux1 = function(i){
      e = Sub[i,]-mj;
      d = sqrt(sum(e^2));
      return(d);
    };
    Dj = unlist(lapply(seq(1:nj),aux1));
    return(mean(Dj));
  }
  # Diameters
  D = unlist(lapply(seq(1:k),aux));
  
  ## Calculate statistic
  aux = function(j){
    # Current cluster
    key1 = labs[j];
    mj = m[[key1]];
    # Index is j, not key1, since D has only length(labs) elements. 
    Dj = D[j];
    aux1 = function(i){
      key2 = labs[i];
      mi = m[[key2]];
      e = (mi-mj);
      d = sqrt(sum(e^2));
      if(d>0){
        return((D[i]+Dj)/d);
      }
    };
    Sij = unlist(lapply(seq(1:k),aux1));
    # Maximum of similarity score
    return(max(Sij));
  }
  Stat = unlist(lapply(seq(1:k),aux));
  # Output
  return(mean(Stat));
}

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
#' @param M Object of class mix. 
#' @return A list containing the cluster quality metrics. 
#' 
#' @importFrom cluster silhouette
#' @importFrom stats dist 
#' @export 
#' @seealso See \code{\link{chooseK}} for using quality metrics to choose the cluster number. 
#' 
#' @examples
#' \dontrun{
#' set.seed(100);
#' # Data generation
#' Y = rGMM(n=500,d=3,k=4);
#' M = fit.GMM(Y=Y,k=3);
#' # Clustering quality
#' Q = clustQual(M);
#' }

clustQual = function(M){
  # Unpack
  a = M@Assignments[,1];
  Y = data.matrix(M@Completed);
  k = length(M@Proportions);
  n = nrow(Y);
  
  # Output structure
  Out = list();
  
  ## BIC
  Out$BIC = log(n)*(3*k)-M@Objective;
  
  ## Calinski-Harabaz Index
  Out$CHI = CalHar(Y=Y,a=a,m=M@Means);
  
  ## Davies-Bouldin Index
  Out$DBI = DavBou(Y=Y,a=a,m=M@Means);
  
  ## Silhouette
  Sil = silhouette(x=a,dist=dist(x=Y,method="euclidean"));
  Out$SIL = mean(Sil[,3]);
  
  # Output
  return(Out);
}

########################
# Cluster Number Selection
########################

#' Cluster Number Selection
#'
#' Function to choose the number of clusters k. Examines cluster numbers between
#' k0 and k1. For each cluster number, generates B bootstrap data sets, fits the
#' Gaussian Mixture Model (\code{\link{fit.GMM}}), and calculates quality
#' metrics (\code{\link{clustQual}}). For each metric, determines the optimal
#' cluster number \code{kopt}, and the \code{k1se}, the smallest cluster number
#'  whose quality is within 1 SE of the optimum.
#'
#' @param Y Numeric data matrix.
#' @param k0 Minimum number of clusters.
#' @param k1 Maximum number of clusters.
#' @param B Bootstrap replicates.
#' @param M0 Optional list of initial mean vectors.
#' @param fix.means Fix the means to their starting value? Must provide initial
#'   values.
#' @param S0 Optional list of initial covariance matrices.
#' @param pi0 Optional vector of initial cluster proportions.
#' @param maxit Maximum number of EM iterations.
#' @param eps Minimum acceptable increment in the EM objective.
#' @param report Report bootstrap progress?
#' @param parallel Run in parallel? Only used if \eqn{k>1}. Must register
#'   parallel backend first.
#' @return List containing \code{Choices}, the recommended number of clusters
#'   according to each quality metric, and \code{Results}, the mean and standard
#'   error of the quality metrics at each cluster number evaluated.
#'   
#' @importFrom foreach foreach registerDoSEQ '%do%' '%dopar%'
#' @importFrom stats var
#' @export 
#' @seealso See \code{\link{clustQual}} for evaluating cluster quality, and \code{\link{fit.GMM}}
#' for estimating the GMM with a specified cluster number. 
#' 
#' @examples 
#' \dontrun{
#' set.seed(100);
#' M = list(c(2,2),c(2,-2),c(-2,2),c(-2,-2));
#' Y = rGMM(n=500,d=2,k=4,M=M);
#' K = chooseK(Y=Y,k0=2,k1=6,B=10,maxit=10,eps=1e-4,report=T);
#' K$Choices;
#' }

chooseK = function(Y,k0=2,k1=NULL,B=100,M0=NULL,fix.means=F,S0=NULL,pi0=NULL,maxit=10,eps=1e-4,report=T,parallel=F){
  # Check inputs
  if(k0<2){stop("At least 2 clusters are required to calculate quality metrics.")};
  if(is.null(k1)){k1=k0+2};
  n = nrow(Y);
  
  # Cluster numbers
  K = seq(from=k0,to=k1);
  nk = length(K);
  
  # Parallelization
  if(!parallel){registerDoSEQ()};
  
  # Loop over cluster numbers
  j = b = 1;
  Results = foreach(j=1:nk,.combine=rbind) %do% {
    # Current cluster number
    k = K[j];
    # Loop over bootstrap replicates
    Boot = foreach(b=1:B,.combine=rbind,.errorhandling="remove") %dopar% {
      if(b==1){
        ## Original sample
        M = tryCatch(expr=fit.GMM(Y=Y,k=k,M0=M0,fix.means=fix.means,S0=S0,pi0=pi0,maxit=maxit,eps=eps,report=F,parallel=F),
                     error=function(cond){return(NULL)});
        if(!is.null(M)){
          Out = unlist(clustQual(M));
        } else {
          Out = rep(NA,4);
        }
      } else {
        ## Bootstrap
        # Draw
        Draw = sort(sample(x=n,size=n,replace=T));
        Yb = Y[Draw,];
        # Fit
        M = tryCatch(expr=fit.GMM(Y=Yb,k=k,M0=M0,fix.means=fix.means,S0=S0,pi0=pi0,maxit=maxit,eps=eps,report=F,parallel=F),
                     error=function(cond){return(NULL)});
        if(!is.null(M)){
          Out = unlist(clustQual(M));
        } else{
          Out = rep(NA,4);
        }
      }; 
      return(Out);
    }; # End bootstrap loop
    
    # Check if any fits were successful
    if(is.null(Boot)){
      cat("No fits succeeded at size",k,"\n");
    } else {
      # Remove NA
      Boot = Boot[complete.cases(Boot),,drop=F];
      nb = nrow(Boot);
      # Report
      if(report){cat("Cluster size",k,"complete.",nb,"fit(s) succeeded.\n")};
      # Output if at least 3 fits were successful
      if(!is.null(nb)&(nb>=3)){
        # Summary statistics
        Means = aaply(.data=Boot,.margins=2,.fun=mean);  
        Vs = aaply(.data=Boot,.margins=2,.fun=var);
        SEs = sqrt(Vs/nb);
        # Output
        Out = data.frame("Clusters"=k,"Fits"=nb,"Metric"=names(Means),"Mean"=Means,"SE"=SEs);
        rownames(Out) = NULL;
        return(Out);
      };
    };
  };
  
  # Check for any results
  if(is.null(Results)){
    stop("Unable to fit sufficient models with the current cluster numbers.");
  };
  
  ## Cluster recommendations
  Metrics = c("BIC","CHI","DBI","SIL");
  Choices = data.frame("Metric"=Metrics,"kopt"=rep(NA,4),"Metric_kopt"=rep(NA,4),
                       "k1se"=rep(NA,4),"Metric_k1se"=rep(NA,4));
  
  # BIC
  BIC = Results[Results$Metric=="BIC",];
  kopt = BIC$Clusters[which.min(BIC$Mean)];
  BICopt = BIC$Mean[BIC$Clusters==kopt];
  se.BICopt = BIC$SE[BIC$Clusters==kopt];
  # Threshold
  BIC = BIC[BIC$Mean<=BICopt+se.BICopt,];
  k1se = min(BIC$Clusters);
  # Store
  Choices[1,2] = kopt;
  Choices[1,3] = BICopt;
  Choices[1,4] = k1se;
  Choices[1,5] = BIC$Mean[BIC$Clusters==k1se];
  rm(kopt,k1se,BIC,BICopt,se.BICopt);
  
  # CHI
  CHI = Results[Results$Metric=="CHI",];
  kopt = CHI$Clusters[which.max(CHI$Mean)];
  CHIopt = CHI$Mean[CHI$Clusters==kopt];
  se.CHIopt = CHI$SE[CHI$Clusters==kopt];
  # Threshold
  CHI = CHI[CHI$Mean>=CHIopt-se.CHIopt,];
  k1se = min(CHI$Clusters);
  # Store
  Choices[2,2] = kopt;
  Choices[2,3] = CHIopt;
  Choices[2,4] = k1se;
  Choices[2,5] = CHI$Mean[CHI$Clusters==k1se];
  rm(kopt,k1se,CHI,CHIopt,se.CHIopt);
  
  # DBI
  DBI = Results[Results$Metric=="DBI",];
  kopt = DBI$Clusters[which.min(DBI$Mean)];
  DBIopt = DBI$Mean[DBI$Clusters==kopt];
  se.DBIopt = DBI$SE[DBI$Clusters==kopt];
  # Threshold
  DBI = DBI[DBI$Mean<=DBIopt+se.DBIopt,];
  k1se = min(DBI$Clusters);
  # Store
  Choices[3,2] = kopt;
  Choices[3,3] = DBIopt;
  Choices[3,4] = k1se;
  Choices[3,5] = DBI$Mean[DBI$Clusters==k1se];
  rm(kopt,k1se,DBI,DBIopt,se.DBIopt);
  
  # SIL
  SIL = Results[Results$Metric=="SIL",];
  kopt = SIL$Clusters[which.max(SIL$Mean)];
  SILopt = SIL$Mean[SIL$Clusters==kopt];
  se.SILopt = SIL$SE[SIL$Clusters==kopt];
  # Threshold
  SIL = SIL[SIL$Mean>=SILopt-se.SILopt,];
  k1se = min(SIL$Clusters);
  # Store
  Choices[4,2] = kopt;
  Choices[4,3] = SILopt;
  Choices[4,4] = k1se;
  Choices[4,5] = SIL$Mean[SIL$Clusters==k1se];
  rm(kopt,k1se,SIL,SILopt,se.SILopt);
  
  # Output
  Out = list();
  Out$Choices = Choices;  
  Out$Results = Results;
  return(Out);
}
