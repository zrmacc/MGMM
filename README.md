---
title: "README"
author: "Zachary McCaw"
date: "2019-01-18"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Introduction](#introduction)
* [Data Generation](#data-generation)
* [Parameter Estimation](#parameter-estimation)

# Introduction

## Overview

Suppose that the data consist of independent random vectors. Each observation belongs to one of several distinct clusters. Conditional on cluster membership, the observation follows a multivariate normal distribution, with cluster-specific mean and covariance. The elements of an observation are subject to arbitrary patterns of missingness at random. Given such data, this package estimates the parameters of a Gaussian Mixture Model (GMM), including the marginal probabilities of cluster membership, and the cluster-specific means and covariances. For each observation, the posterior probabiity of membership to each cluster is calculated, and a maximum a posteriori classification is provided. 

## Model

Suppose that the data consist of $n$ random vectors $y_{i}$ in $\mathbb{R}^{p}$. Each observation arises from one of $k$ distinct clusters. Let $z_{ij} = 1$ if observation $i$ belongs to cluster $j$. The marginal probability of membership the $j$th cluster is $\pi_{j} = P[z_{i}=j]$. Conditional on membership to the $j$th cluster, $y_{i}$ follows a multivariate normal distribution, with mean $\mu_{j}$ and covariance $\Sigma_{j}$. The generative model is:

$$
z_{i} \sim \text{Multinomial}[\pi_{1},\cdots,\pi_{k}]
$$

$$
y_{i}\big|(z_{i}=j) \sim N\big(\mu_{j},\Sigma_{j}\big)
$$

The EM algorithm is used to obtain maximum likelihood estimates (MLEs) of the model parameters. The procedure treats both the unobserved elements of each random vector $y_{i}$, and the unobserved cluster assignments $z_{ij}$ as missing data. From the MLEs, posterior probabilities of cluster membership are calculated as follows. Partition each observation $y_{i}$ as $(s_{i},\ t_{i})$, where $s_{i}$ denotes the observed elements, and $t_{i}$ denotes the missing elements. The posterior probability of membership to cluster $j$, given the observed data $s_{i}$, is:

$$
\gamma_{ij} = P[z_{i}=j|s_{i}] = \frac{f(s_{i}|\ \mu_{j},\Sigma_{j})\pi_{j}}{\sum_{l=1}^{k}f(s_{i}|\ \mu_{l},\Sigma_{l})\pi_{l}}
$$

The maximum a posteriori classificaiton for $y_{i}$ is:

$$
\arg\max_{j}\ \hat{\gamma}_{ij}
$$

# Data Generation

## Description

The function `rMNMix` simulates observations from a Gaussian Mixture Model. The number of observations is specified by `n`, and the dimension of each observation by `d`. The number of clusters is set using `k`, which defaults to one. The marginal probabilities of cluster membership are provided as a numeric vector `pi`, which should contain `k` elements. If $k>0$ but `pi` is omitted, the clusters are taken as equi-probable. The proportion of elements in the $n \times d$ data matrix that are missing is specified by `m`, which defaults to zero. Note that when $m>0$ it is possible for all elements of an observation to go missing. The cluster means `M` are provided as a numeric prototype vector, or a list of such vectors. If a single prototype is provided, that vector is used as the mean for all clusters. By default, the zero vector is adopted as the prototype. The cluster covariances `S` are provided as a numeric matrix, or a list of such matrices. If a single prototype is provided, that matrix is used as the covariance for all clusters. By default, the identity matrix is adopted as the prototype. 

## Examples

### Single Component without Missingness

In this example, $10^{3}$ observations are simulated from a single `k=1` bivariate normal distribution `d=2` without missingness. The mean is $\mu=(2,2)$, and the covariance is an exchangeable correlation structure with off-diagonal $\rho=0.5$.  


```r
# Single component without missingness
Sigma = matrix(c(1,0.5,0.5,1),nrow=2);
Y = rMNMix(n=1e3,d=2,k=1,M=c(2,2),S=Sigma);
```

### Single Component with Missingness 

In this example, $10^{3}$ observations are simulated from a single `k=1` trivariate normal distribution `d=3` with 20% missingness `m=0.2`. The mean defaults to the zero vector, and the covariance to the identity matrix. 


```r
# Single component with missingness
Y = rMNMix(n=1e3,d=3,k=1,m=0.2);
```

### Two Components without Missingness

In this example, $10^{3}$ observations are simulated from a two-coomponent `k=2` trivariate normal distribution `d=3` without missingness. The mean vectors are $\mu_{1}=(-2,-2,-2)$ and $\mu_{2}=(2,2,2)$. The covariance matrices are both exchangeable with off-diagonal $\rho=0.5$. Since `pi` is omitted, the cluster are equi-probable, i.e. $\pi_{1}=\pi_{2}=1/2$. 


```r
# Two-component mixture without missingness
M = list(c(-2,-2,-2),c(2,2,2));
Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
Y = rMNMix(n=1e3,d=3,k=2,M=M,S=Sigma);
```

### Four Components with Missingness

In this example, $10^{3}$ observations are simulated from a four-coomponent `k=4` bivariate normal distribution `d=2` with 10% missingness `m=0.1`. The mean vectors are $\mu_{1}=(-2,-2)$, $\mu_{2}=(-2,2)$, $\mu_{3}=(2,-2)$ and $\mu_{4}=(2,2)$. The covariance matrices are all $0.5*I$. The cluster proportions are (35%, 15%, 15%, 35%) for $(\pi_{1},\pi_{2},\pi_{3},\pi_{4})$, respectively. 


```r
# Four-component mixture with missingness
M = list(c(-2,-2),c(-2,2),c(2,-2),c(2,2));
S = 0.5*diag(2);
Y = rMNMix(n=1e3,d=2,k=4,pi=c(0.35,0.15,0.15,0.35),m=0.1,M=M,S=S);
```

# Parameter Estimation

## Description

The function `fit.MNMix` estimates parameters for the Gaussian Mixture Model. The data are expected as a numeric matrix `Y`, with observations as rows. The number of mixture components is specified using `k`, which defaults to one. Initial values for the mean vectors, covariance matrices, and cluster proportions are provided using `M0`, `S0`, and `pi0`, respectively. If the data `Y` contain complete observations, i.e. observations with no missing elements, `fit.MNMix` will attempt to initialize the model parameters ($\mu,\Sigma,\pi$) via K-means. However, if the data `Y` contain no complete observations, then initial values are required for each of `M0`, `S0`, and `pi0`. Supplying initial values may also result in more accurate estimates when there are relatively few complete observations. The initial means `M0` are provided as a list of vectors, the covariances `S0` as a list of matrices, and the cluster proportions `pi0` as a numeric vector. Note that `M0` and `S0` are expected as lists even if the model only contains a single component `k=1`.

The arguments `maxit`, `eps`, `report`, and `parallel` control the fitting procedure. `maxit` sets the maximum number of EM iterations to attempt. The default is $10^{2}$. `eps` sets the minimum acceptable improvement in the EM objective function. The default is $10^{-6}$. If `report=TRUE`, then fitting progress is displayed. For models with missingness and more than one mixture component `k>1`, setting `parallel=TRUE` may improve run time. The parallel backend must be registered beforehand. 

## Examples

### Single Component without Missingness

In this example, $10^{3}$ observations are simulated with a single bivariate normal distribution without missingness. Since the model contains only a single component, the output is a list containing the estimated mean and covariance. In the case of a single component without missingness, the maximum likelihood estimates are available in closed form.  


```r
# Single component without missingness
Sigma = matrix(c(1,0.5,0.5,1),nrow=2);
Y = rMNMix(n=1e3,d=2,k=1,M=c(2,2),S=Sigma);
M = fit.MNMix(Y=Y,k=1);
cat("Estimated Mean and Covariance:\n");
show(M);
```

```
## Estimated Mean and Covariance:
## $Mean
##       y1       y2 
## 1.997693 2.010779 
## 
## $Covariance
##           y1        y2
## y1 1.0335551 0.5215485
## y2 0.5215485 1.0504981
## 
## $Objective
## [1] -1792.507
```

### Single Component with Missingness

In this example, $10^{3}$ observations are simulated from a single trivariate normal distribution with 20% missingness. Since the model contains only a single component, the output is again a list. However, in the case of missingness, the EM algorithm is used for estimation. In addition to the estimated mean and covariance, the output now contains the final EM objective. The true mean is the zero vector, and the true covariance is identity. For `M1` below, the initial mean and covariance are estimated internally using complete observations. For `M2` below, the mean and covariance are initialized at the truth. The final value of the EM objective is increased by initializing at the truth. 


```r
# Single component with missingness
set.seed(100);
Y = rMNMix(n=1e3,d=3,k=1,m=0.2);
cat("Initial parameter values set internally:\n");
M1 = fit.MNMix(Y=Y,k=1);
cat("\n");
show(M1);
cat("Initial parameter values set manually:\n");
m0 = rep(0,3);
S0 = diag(3);
M2 = fit.MNMix(Y=Y,k=1,M0=list(m0),S0=list(S0));
cat("\n");
show(M2);
cat("Gain in final objective by initializing parameters at the truth:\n")
M2$Objective-M1$Objective;
```

```
## Initial parameter values set internally:
## Objective increment:  5.81 
## Objective increment:  0.662 
## Objective increment:  0.0976 
## Objective increment:  0.0143 
## Objective increment:  0.00147 
## 5 update(s) performed before tolerance limit.
## 
## 
## $Mean
##          y1          y2          y3 
##  0.05363057 -0.03032985  0.02118261 
## 
## $Covariance
##              y1            y2            y3
## y1  0.995409482 -0.0272300548  0.0084440124
## y2 -0.027230055  1.0377546994 -0.0008892942
## y3  0.008444012 -0.0008892942  1.0590200651
## 
## $Objective
## [1] -3064.304
## 
## Initial parameter values set manually:
## 0 update(s) performed before tolerance limit.
## 
## 
## $Mean
## [1] 0 0 0
## 
## $Covariance
##      [,1] [,2] [,3]
## [1,]    1    0    0
## [2,]    0    1    0
## [3,]    0    0    1
## 
## $Objective
## NULL
## 
## Gain in final objective by initializing parameters at the truth:
## numeric(0)
```

### Two Components without Missingness

In this example, $10^{3}$ observations are simulated from a two-component trivariate normal distribution without missingness. Since the model has multiple components, the output is an object of class `mix`. The show method displays the estimated cluster proportions and the final objective. The slots `@Means` and `@Covariances` contain lists of the estimated cluster means and covariances. The posterior probability of membership to each cluster is contained in the `@Responsibilities` slot, and the highest posterior probability classification of each observation is contained in the `@Assignments` slot. 


```r
# Two componets without missingness
Means = list(c(-2,-2,-2),c(2,2,2));
Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
Y = rMNMix(n=1e3,d=3,k=2,M=Means,S=Sigma);
M = fit.MNMix(Y=Y,k=2,eps=1e-8);
cat("\n");
show(M);
cat("Cluster means:\n");
M@Means;
cat("Cluster covariances:\n");
M@Covariances;
cat("Cluster responsibilities:\n");
head(M@Responsibilities);
cat("\nCluster assignments:\n");
head(M@Assignments);
```

```
## Objective increment:  0.396 
## Objective increment:  0.034 
## Objective increment:  0.00965 
## Objective increment:  0.00282 
## Objective increment:  0.000957 
## Objective increment:  0.000395 
## Objective increment:  0.000189 
## Objective increment:  9.75e-05 
## Objective increment:  5.22e-05 
## Objective increment:  2.83e-05 
## Objective increment:  1.54e-05 
## Objective increment:  8.43e-06 
## Objective increment:  4.61e-06 
## Objective increment:  2.52e-06 
## Objective increment:  1.38e-06 
## Objective increment:  7.56e-07 
## Objective increment:  4.14e-07 
## Objective increment:  2.27e-07 
## Objective increment:  1.24e-07 
## Objective increment:  6.79e-08 
## Objective increment:  3.72e-08 
## Objective increment:  2.03e-08 
## Objective increment:  1.11e-08 
## Objective increment:  6.1e-09 
## 23 update(s) performed before tolerance limit.
## 
## 
## Normal Mixture Model with 2 Components. 
## Cluster Proportions:
##    k1    k2 
## 0.518 0.482 
## 
## Final Objective:
##      k1 
## -2982.3 
## 
## Cluster means:
## [[1]]
##        y1        y2        y3 
## -1.949991 -1.912998 -1.982099 
## 
## [[2]]
##       y1       y2       y3 
## 1.964810 1.933682 1.956676 
## 
## Cluster covariances:
## [[1]]
##           y1        y2        y3
## y1 0.9347688 0.5095796 0.4896254
## y2 0.5095796 1.0477300 0.4895118
## y3 0.4896254 0.4895118 0.9365137
## 
## [[2]]
##          y1        y2        y3
## y1 1.070444 0.5310510 0.5222580
## y2 0.531051 1.0483091 0.5286882
## y3 0.522258 0.5286882 1.0013230
## 
## Cluster responsibilities:
##             k1           k2
## 1 9.997767e-01 2.232718e-04
## 2 1.693337e-03 9.983067e-01
## 3 9.999994e-01 6.164588e-07
## 4 9.910511e-01 8.948881e-03
## 5 3.285010e-05 9.999671e-01
## 6 4.350046e-08 1.000000e+00
## 
## Cluster assignments:
##   Assignment
## 1          1
## 2          2
## 3          1
## 4          1
## 5          2
## 6          2
```

### Four Components with Missingness

In this example, $10^{3}$ observations are simulated from a four-component bivariate normal distribution with 10% missingness. Since the model has multiple components, the output is an object of class `mix`.


```r
set.seed(200);
# Four components with missingness
M = list(c(2,2),c(2,-2),c(-2,2),c(-2,-2));
S = 0.5*diag(2);
Y = rMNMix(n=1000,d=2,k=4,pi=c(0.35,0.15,0.15,0.35),m=0.1,M=M,S=S);
M = fit.MNMix(Y=Y,k=4,eps=1e-8);
show(M);
cat("Cluster means:\n");
M@Means;
cat("\nCluster assignments:\n");
head(M@Assignments);
```

```
## Objective increment:  6980 
## Objective increment:  168 
## Objective increment:  90.9 
## Objective increment:  57.4 
## Objective increment:  36.8 
## Objective increment:  26.9 
## Objective increment:  20.6 
## Objective increment:  17.1 
## Objective increment:  15.1 
## Objective increment:  14.1 
## Objective increment:  14.1 
## Objective increment:  14.3 
## Objective increment:  15.1 
## Objective increment:  15.4 
## Objective increment:  15 
## Objective increment:  13.2 
## Objective increment:  10.7 
## Objective increment:  7.68 
## Objective increment:  5.26 
## Objective increment:  3.31 
## Objective increment:  2.02 
## Objective increment:  1.12 
## Objective increment:  0.574 
## Objective increment:  0.207 
## 24 update(s) performed before tolerance limit.
## 
## Normal Mixture Model with 4 Components. 
## Cluster Proportions:
##     k1     k2     k3     k4 
## 0.0863 0.3220 0.2190 0.3730 
## 
## Final Objective:
##       k1 
## -2649.01 
## 
## Cluster means:
## [[1]]
##         y1         y2 
## -0.9088869  0.4385624 
## 
## [[2]]
##       y1       y2 
## 1.994131 1.992714 
## 
## [[3]]
##          y1          y2 
##  0.03938591 -0.03059968 
## 
## [[4]]
##        y1        y2 
## -2.038837 -2.006364 
## 
## 
## Cluster assignments:
##   Assignment
## 1          2
## 2          3
## 3          4
## 4          4
## 5          2
## 6          2
```
