---
title: "README"
author: "Zachary McCaw"
date: "2018-08-25"
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

Suppose that the data consists of multivariate normal random vectors. Each observation is thought to arise from one of several distinct clusters. The collection of observations belonging to a given cluster follow a multivariate normal distribution, with cluster-specific means and covariances. The elements of an observation are subject to arbitrary patterns of missingness. However, whether or not an element is missing is taken as independent of that element’s value, given the observed elements. Provided such data, and the number of clusters in the population, this package estimates the cluster-specific means, covariances, and marginal membership probabilities. In addition, for each observation, the posterior probability of membership is calculated for all clusters. A classification of the observation into the cluster with the highest posterior probability is provided. 

## Model

Suppose that the data consist of random vectors $y_{i}$ in $\mathbb{R}^{p}$. Each observation belongs to one of $k$ clusters. Let $z_{i}\in\{1,\cdots,k\}$ denote the cluster to which $y_{i}$ belongs. The marginal probability of membership the $j$th cluster is $\pi_{j} = P[z_{i}=j]$. Conditional on membership to cluster $j$, $z_{i}$ follows a multivariate normal distribution with mean $\mu_{j}$ and covariance $\Sigma_{j}$. The generative model is:

$$
z_{i} \sim \text{Multinomial}[\pi_{1},\cdots,\pi_{k}]
$$

$$
y_{i}\big|(z_{i}=j) \sim N\big(\mu_{j},\Sigma_{j}\big)
$$

The EM algorithm is used to obtain maximum likelihood estimates of the model parameters. Posterior probabilities of cluster membership are calculated as follows. Partition each observation $y_{i}$ as $(s_{i},\ t_{i})$, where $s_{i}$ denotes the observed elements, and $t_{i}$ the missing elements. The posterior probability of membership to cluster $j$, given the observed data $s_{i}$, is:

$$
\gamma_{ij} = P[z_{i}=j|s_{i}] = \frac{f(s_{i}|\ \mu_{j},\Sigma_{j})\pi_{j}}{\sum_{l=1}^{k}f(s_{i}|\ \mu_{l},\Sigma_{l})\pi_{l}}
$$

The estimated responsibility of the $j$th cluster for the $i$th observation is obtained by substituting in the MLEs. The predicted cluster membership of observation $y_{i}$ is:

$$
\hat{z}_{i} = \arg\max_{j}\ \hat{\gamma}_{ij}
$$

# Data Generation

## Description

The function `rMNMix` simulates observations from a mixture of multivariate normal distributions. The number of observations is specified by `n`, and the dimension of each observation by `d`. The number of clusters is set using `k`, which defaults to one. The marginal probabilities of cluster membership are provided as a numeric vector `pi`, which should contain `k` elements. If $k>0$ but `pi` is omitted, the clusters are taken as equi-probable. The proportion of elements in the $n \times d$ data matrix that are missing is specified by `m`, which defaults to zero. Note that when $m>0$ it is possible for all elements of an observation to go missing. The cluster means `M` are provided as a numeric prototype vector, or a list of such vectors. If a single prototype is provided, that vector is used as the mean for all clusters. By default, the zero vector is adopted as the prototype. The cluster covariances `S` are provided as a numeric matrix, or a list of such matrices. If a single prototype is provided, that matrix is used as the covariance for all clusters. By default, the identity matrix is adopted as the prototype. 

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

The function `fit.MNMix` estimates parameters for multivariate normal mixture model. The data are expected as a numeric matrix `Y`, with observations as rows. The number of mixture components is specified using `k`, which defaults to one. Initial values for the mean vectors, covariance matrices, and cluster proportions are provided using `M0`, `S0`, and `pi0`, respectively. If the data `Y` contain complete observations, i.e. observations with no missing elements, `fit.MNMix` will attempt to initialize the model parameters ($\mu,\Sigma,\pi$) via K-means. However, if the data `Y` contain no complete observations, then initial values are required for each of `M0`, `S0`, and `pi0`. Supplying initial values may also result in more accurate estimates when there are relatively few complete observations. The initial means `M0` are provided as a list of vectors, the covariances `S0` as a list of matrices, and the cluster proportions `pi0` as a numeric vector. Note that `M0` and `S0` are expected as lists even if the model only contains a single component `k=1`.

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
## 2.050852 2.012795 
## 
## $Covariance
##           y1        y2
## y1 0.9965648 0.5194510
## y2 0.5194510 0.9747743
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
## Objective increment: 1.32
## Objective increment: 0.0242
## 2 update(s) performed before tolerance limit. 
## 
## $Mean
##          y1          y2          y3 
##  0.05366774 -0.03125547  0.02210157 
## 
## $Covariance
##             y1           y2          y3
## y1  0.99206733 -0.026448480 0.014700353
## y2 -0.02644848  1.035613830 0.007076264
## y3  0.01470035  0.007076264 1.051829916
## 
## $Objective
## [1] -1526.005
## 
## Initial parameter values set manually:
## Objective increment: 2.99
## 1 update(s) performed before tolerance limit. 
## 
## $Mean
##          y1          y2          y3 
##  0.04401760 -0.02438124  0.01684148 
## 
## $Covariance
##              y1            y2           y3
## y1  0.991860758 -0.0157399011 0.0067010171
## y2 -0.015739901  1.0272907159 0.0002498111
## y3  0.006701017  0.0002498111 1.0439131667
## 
## $Objective
## [1] -1518.475
## 
## Gain in final objective by initializing parameters at the truth:
## [1] 7.529237
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
## Objective increment: 0.198
## 1 update(s) performed before tolerance limit. 
## 
## Normal Mixture Model with 2 Components. 
## Cluster Proportions:
##    k1    k2 
## 0.518 0.482 
## 
## Final Objective:
##        Q 
## -1833.56 
## 
## Cluster means:
## [[1]]
##        y1        y2        y3 
## -1.949637 -1.912012 -1.981330 
## 
## [[2]]
##       y1       y2       y3 
## 1.966730 1.934883 1.958164 
## 
## Cluster covariances:
## [[1]]
##           y1        y2        y3
## y1 0.9359385 0.5095797 0.4897276
## y2 0.5095797 1.0466288 0.4919448
## y3 0.4897276 0.4919448 0.9350531
## 
## [[2]]
##           y1        y2        y3
## y1 1.0647953 0.5277063 0.5179572
## y2 0.5277063 1.0457644 0.5274515
## y3 0.5179572 0.5274515 0.9961200
## 
## Cluster responsibilities:
##   ID           k1           k2
## 1  1 9.997793e-01 2.206736e-04
## 2  2 1.701449e-03 9.982986e-01
## 3  3 9.999994e-01 5.746544e-07
## 4  4 9.912762e-01 8.723782e-03
## 5  5 3.309337e-05 9.999669e-01
## 6  6 4.419866e-08 1.000000e+00
## 
## Cluster assignments:
##   ID Assignment
## 1  1          1
## 2  2          2
## 3  3          1
## 4  4          1
## 5  5          2
## 6  6          2
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
## Objective increment: 0.412
## Objective increment: 1.34
## Objective increment: 0.167
## Objective increment: 0.0244
## Objective increment: 0.00449
## Objective increment: 0.001
## Objective increment: 0.000246
## Objective increment: 6.26e-05
## Objective increment: 1.57e-05
## Objective increment: 3.79e-06
## Objective increment: 8.43e-07
## Objective increment: 1.55e-07
## Objective increment: 1.26e-08
## 13 update(s) performed before tolerance limit. 
## Normal Mixture Model with 4 Components. 
## Cluster Proportions:
##    k1    k2    k3    k4 
## 0.326 0.146 0.153 0.374 
## 
## Final Objective:
##        Q 
## -1570.82 
## 
## Cluster means:
## [[1]]
##       y1       y2 
## 2.046009 1.991048 
## 
## [[2]]
##        y1        y2 
## -1.895473  2.031035 
## 
## [[3]]
##        y1        y2 
##  2.027171 -2.094533 
## 
## [[4]]
##        y1        y2 
## -2.024894 -1.997215 
## 
## 
## Cluster assignments:
##   ID Assignment
## 1  1          3
## 2  2          1
## 3  3          2
## 4  4          2
## 5  5          2
## 6  6          1
```
