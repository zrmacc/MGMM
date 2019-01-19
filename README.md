---
title: "README"
author: "Zachary McCaw"
date: "2019-01-19"
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

This package implements clustering of multivariate normal random vectors with missing elements. Clustering is achieved by fitting a Gaussian Mixture Model (GMM). The parameters are estimated by maximum likelihood, using the Expectation Maximization (EM) algorithm. Our implementation extends existing methods by allowing for missingness in the input data. The EM algorithm addresses missingness of both the cluster assignments and the vector components. The output includes the marginal cluster membership probabilities; the mean and covariance of each cluster; the density of each mixture component evaluated on the observations; the posterior probabilities of cluster membership; maximum *a posteriori* cluster assignments; and a completed version of the input data, with missing values replaced by their posterior expectations. 

## Model

Suppose the data consist of $n$ random vectors $\{y_{i}\}_{i=1}^{n}$ in $\mathbb{R}^{d}$. Each observation $y_{i}$ arises from one of $k$ distinct clusters. Associate with each observation a $k\times 1$ vector of indicators $z_{i}$, where $z_{ij} = 1$ if observation $i$ belongs to cluster $j$. The marginal probability of membership the $j$th cluster is $\pi_{j} = P[z_{ij}=1]$. Conditional on membership to the $j$th cluster, $y_{i}$ follows a multivariate normal distribution, with mean $\mu_{j}$ and covariance $\Sigma_{j}$. The generative model is:

$$
z_{i} \sim \text{Multinomial}[\pi_{1},\cdots,\pi_{k}]
$$

$$
y_{i}\big|(z_{ij}=1) \sim N\big(\mu_{j},\Sigma_{j}\big)
$$

The marginal distribution of the observations is a GMM:
$$
f(y_{i}) = \sum_{j=1}^{k}f(y_{i},z_{ij}=1) = \sum_{j=1}^{k}f(y_{i}|z_{ij}=1)P(z_{ij}=1) = \sum_{j=1}^{k}f(y_{j}|\mu_{j},\Sigma_{j})\pi_{j}
$$

Each element $y_{ik}$ of $y_{i}$ is potentially missing at random (MAR). Associated with each observation a $d\times 1$ vector of indicators $R_{i}$, where $R_{ik}=1$ if element $k$ of observation $i$ is observed. Partition each $y_{i}$ into its observed components $s_{i}$ and its unobserved components $t_{i}$. That is, $y_{ik}$ belongs to $s_{i}$ if $R_{ik}=1$, and belongs to $t_{i}$ if $R_{ik}=0$. The missingness is at random if $(R_{ik}\perp y_{ik})|s_{i}$ for each $y_{ik}$ in $t_{i}$. 

Maximum likelihood estimates (MLEs) for the parameters of the GMM are obtained using the EM algorithm. During the E-step, both the cluster assignments $z_{ij}$ and the unobserved components $t_{i}$ of $y_{i}$ are treated as missing data. Suppose temporarily that all data were observed for observation $i$. The contribution of subject $i$ to the *complete data* log likelihood would be: 
$$
\ell_{i} = \sum_{j=1}^{k}z_{ij}\ln\pi_{j}-\frac{1}{2}\sum_{j=1}^{k}z_{ij}\ln\det(\Sigma_{j})-\frac{1}{2}\sum_{j=1}^{k}z_{ij}(y_{i}-\mu_{j})\Sigma_{j}^{-1}(y_{i}-\mu_{j})
$$

Since $z_{ij}$ is never observed, and $y_{i}$ is incompletely observed, the complete data log likelihood cannot be evaluated. In the E-step, the EM objective is formed by taking the expectation of the complete data log likelihood, conditional on the observed data and the current parameter state:
$$
q_{i}^{(r)} \equiv E[\ell_{i}|s_{i},\vartheta^{(r)}] = \sum_{j=1}^{k}\gamma_{ij}^{(r)}\ln\pi_{j}-\frac{1}{2}\sum_{j=1}^{k}\gamma_{ij}^{(r)}\ln\det(\Sigma_{j})-\frac{1}{2}\sum_{j=1}^{k}\text{tr}(\Sigma_{j}^{-1}V_{ij}^{(r)})
$$

Here $s_{i}$ is the observed data for subject $i$; $\vartheta^{(r)}$ is the current parameter state; $\gamma_{ij}^{(r)}$ is the *responsibility*, defined as $E[z_{ij}|s_{i},\vartheta^{(r)}]$; and $V_{ij}^{(r)}$ is the *expected residual outer product*, defined as $E[z_{ij}(y_{i}-\mu_{j})(y_{i}-\mu_{j})'|s_{i},\vartheta^{(r)}]$. In the M-step, updates of the model parameters $\vartheta$ are obtained by maximizing the EM objective. 

Once the convergence criterion has been achieved, the responsibility, or posterior probability of cluster membership, is calculated as:
$$
\gamma_{ij} = P[z_{ij}=1|s_{i}] = \frac{f(s_{i}|\ \mu_{j},\Sigma_{j})\pi_{j}}{\sum_{l=1}^{k}f(s_{i}|\ \mu_{l},\Sigma_{l})\pi_{l}}
$$

The maximum a posteriori classification for $y_{i}$ is given by:
$$
A_{i} = \arg\max_{j}\ \gamma_{ij}
$$

For observation $y_{i}$, posterior expectation of the missing elements $t_{i}$, given the observed elements $s_{i}$, is:
$$
E[t_{i}|s_{i}] = \sum_{j=1}^{k}E[t_{i}|s_{i},z_{ij}=1]\pi_{j}
$$

# Data Generation

## Description

The function `rGMM` simulates observations from a Gaussian Mixture Model. The number of observations is specified by `n`, and the dimension of each observation by `d`. The number of clusters is set using `k`, which defaults to one. The marginal probabilities of cluster membership are provided as a numeric vector `pi`, which should contain `k` elements. If $k>0$ but `pi` is omitted, the clusters are assumed equi-probable. The proportion of elements in the $n \times d$ data matrix that are missing is specified by `m`, which defaults to zero. Note that when $m>0$ it is possible for all elements of an observation to be missing. The cluster means `M` are provided either as a numeric prototype vector, or a list of such vectors. If a single prototype is provided, that vector is used as the mean for all clusters. By default, the zero vector is adopted as the prototype. The cluster covariances `S` are provided as a numeric prototype matrix, or a list of such matrices. If a single prototype is provided, that matrix is used as the covariance for all clusters. By default, the identity matrix is adopted as the prototype. 

## Examples

### Single Component without Missingness

In this example, $10^{3}$ observations are simulated from a single `k=1` bivariate normal distribution `d=2` without missingness. The mean is $\mu=(2,2)$, and the covariance is an exchangeable correlation structure with off-diagonal $\rho=0.5$. 


```r
set.seed(100);
# Single component without missingness
Sigma = matrix(c(1,0.5,0.5,1),nrow=2);
Y = rGMM(n=1e3,d=2,k=1,M=c(2,2),S=Sigma);
```

### Single Component with Missingness 

In this example, $10^{3}$ observations are simulated from a single `k=1` trivariate normal distribution `d=3` with 20% missingness `m=0.2`. The mean defaults to the zero vector, and the covariance to the identity matrix. 


```r
# Single component with missingness
Y = rGMM(n=1e3,d=3,k=1,m=0.2);
```

### Two Components without Missingness

In this example, $10^{3}$ observations are simulated from a two-component `k=2` trivariate normal distribution `d=3` without missingness. The mean vectors are $\mu_{1}=(-2,-2,-2)$ and $\mu_{2}=(2,2,2)$. The covariance matrices are both exchangeable with off-diagonal $\rho=0.5$. Since `pi` is omitted, the cluster are equi-probable, i.e. $\pi_{1}=\pi_{2}=1/2$. 


```r
# Two-component mixture without missingness
M = list(c(-2,-2,-2),c(2,2,2));
Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
Y = rGMM(n=1e3,d=3,k=2,M=M,S=Sigma);
```

### Four Components with Missingness

In this example, $10^{3}$ observations are simulated from a four-component `k=4` bivariate normal distribution `d=2` with 10% missingness `m=0.1`. The mean vectors are $\mu_{1}=(-2,-2)$, $\mu_{2}=(-2,2)$, $\mu_{3}=(2,-2)$ and $\mu_{4}=(2,2)$. The covariance matrices are all $0.5*I$. The cluster proportions are (35%, 15%, 15%, 35%) for $(\pi_{1},\pi_{2},\pi_{3},\pi_{4})$, respectively. 


```r
# Four-component mixture with missingness
M = list(c(-2,-2),c(-2,2),c(2,-2),c(2,2));
S = 0.5*diag(2);
Y = rGMM(n=1e3,d=2,k=4,pi=c(0.35,0.15,0.15,0.35),m=0.1,M=M,S=S);
```

# Parameter Estimation

## Description

The function `fit.GMM` estimates parameters of the GMM. The data are expected as a numeric matrix `Y`, with observations as rows. The number of mixture components is specified using `k`, which defaults to one. Initial values for the mean vectors, covariance matrices, and cluster proportions are provided using `M0`, `S0`, and `pi0`, respectively. If the data `Y` contain complete observations, i.e. observations with no missing elements, `fit.GMM` will attempt to initialize the model parameters ($\mu,\Sigma,\pi$) via K-means. However, if the data `Y` contain no complete observations, then initial values are required for each of `M0`, `S0`, and `pi0`. Supplying initial values may also result in more accurate estimates when there are relatively few complete observations. The initial means `M0` are provided as a list of vectors, the covariances `S0` as a list of matrices, and the cluster proportions `pi0` as a numeric vector. Note that `M0` and `S0` are expected as lists even if the model only contains a single component `k=1`.

The arguments `maxit`, `eps`, `report`, and `parallel` control the fitting procedure. `maxit` sets the maximum number of EM iterations to attempt. The default is $10^{2}$. `eps` sets the minimum acceptable improvement in the EM objective function. The default is $10^{-6}$. If `report=TRUE`, then fitting progress is displayed. For models with missingness and more than one mixture component `k>1`, setting `parallel=TRUE` may improve run time. The parallel backend must be registered beforehand. 


## Examples

### Single Component without Missingness

In this example, $10^{3}$ observations are simulated with a single bivariate normal distribution without missingness. Since the model contains only a single component, the output is a list containing the estimated mean and covariance. In the case of a single component without missingness, the maximum likelihood estimates are available in closed form.  


```r
# Single component without missingness
Sigma = matrix(c(1,0.5,0.5,1),nrow=2);
Y = rGMM(n=1e3,d=2,k=1,M=c(2,2),S=Sigma);
M = fit.GMM(Y=Y,k=1);
cat("Estimated Mean and Covariance:\n");
show(M);
```

```
## Estimated Mean and Covariance:
## $Mean
##       y1       y2 
## 2.000419 1.975837 
## 
## $Covariance
##           y1        y2
## y1 0.9712543 0.4798659
## y2 0.4798659 1.0434352
## 
## $Objective
## [1] -1756.026
```

### Single Component with Missingness

In this example, $10^{3}$ observations are simulated from a single trivariate normal distribution with 20% missingness. Since the model contains only a single component, the output is again a list. However, in the case of missingness, the EM algorithm is used for estimation, and a completed version of the input data is returned, with missing values replaced by their posterior expectations.The true mean is the zero vector, and the true covariance is identity. For `M1` below, the initial mean and covariance are estimated internally using complete observations. For `M2` below, the mean and covariance are initialized at the truth. The final value of the EM objective is increased by initializing at the truth. 


```r
set.seed(101);
# Single component with missingness
Y = rGMM(n=1e3,d=3,k=1,m=0.2);
cat("Initial parameter values set internally:\n");
M1 = fit.GMM(Y=Y,k=1);
cat("\nEstimated mean:\n");
show(M1$Mean);
cat("\nEstimated covariance:\n");
show(M1$Covariance);
cat("\nFinal objective:\n");
show(M1$Objective);
cat("Initial parameter values set manually:\n");
m0 = rep(0,3);
S0 = diag(3);
M2 = fit.GMM(Y=Y,k=1,M0=list(m0),S0=list(S0));
cat("\nEstimated mean:\n");
show(M1$Mean);
cat("\nEstimated covariance:\n");
show(M1$Covariance);
cat("\nFinal objective:\n");
show(M1$Objective);
cat("\nGain in final objective by initializing parameters at the truth:\n")
M2$Objective-M1$Objective;
```

```
## Initial parameter values set internally:
## Objective increment:  6.41 
## Objective increment:  0.633 
## Objective increment:  0.0946 
## Objective increment:  0.0168 
## Objective increment:  0.00312 
## Objective increment:  0.000556 
## Objective increment:  7.71e-05 
## 7 update(s) performed before tolerance limit.
## 
## 
## Estimated mean:
##          y1          y2          y3 
##  0.02821992 -0.01186953 -0.02684324 
## 
## Estimated covariance:
##            y1           y2           y3
## y1 1.02062945  0.069358137  0.021234977
## y2 0.06935814  1.074999999 -0.005636683
## y3 0.02123498 -0.005636683  0.980448577
## 
## Final objective:
## [1] -3034.353
## Initial parameter values set manually:
## Objective increment:  3.97 
## Objective increment:  0.185 
## Objective increment:  0.00439 
## 3 update(s) performed before tolerance limit.
## 
## 
## Estimated mean:
##          y1          y2          y3 
##  0.02821992 -0.01186953 -0.02684324 
## 
## Estimated covariance:
##            y1           y2           y3
## y1 1.02062945  0.069358137  0.021234977
## y2 0.06935814  1.074999999 -0.005636683
## y3 0.02123498 -0.005636683  0.980448577
## 
## Final objective:
## [1] -3034.353
## 
## Gain in final objective by initializing parameters at the truth:
## [1] 0.09823569
```

### Two Components without Missingness

In this example, $10^{3}$ observations are simulated from a two-component, trivariate normal distribution without missingness. Since the model has multiple components, the output is an object of class `mix`. The *show method* displays the estimated cluster proportions and the final objective. The slots of the output contain the following:
* `@Means` and `@Covariances`: lists of the estimated cluster means and covariances.
* `@Density`: the cluster densities evaluated at the observations.
* `@Responsibilities`: the posterior membership probabilities for each observation.
* `@Assignments`: the maximum a posteriori cluster assignments.
* `@Completed`: a completed version of the input data is returned, with missing values replaced by their posterior expectations


```r
# Two componets without missingness
Means = list(c(-2,-2,-2),c(2,2,2));
Sigma = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3);
Y = rGMM(n=1e3,d=3,k=2,M=Means,S=Sigma);
M = fit.GMM(Y=Y,k=2,maxit=10,eps=1e-8);
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
## Objective increment:  0.456 
## Objective increment:  0.0227 
## Objective increment:  0.00479 
## Objective increment:  0.00161 
## Objective increment:  0.000641 
## Objective increment:  0.000276 
## Objective increment:  0.000124 
## Objective increment:  5.62e-05 
## Objective increment:  2.58e-05 
## Objective increment:  1.18e-05 
## 10 update(s) performed without reaching tolerance limit.
## 
## 
## Normal Mixture Model with 2 Components. 
## Cluster Proportions:
##    k1    k2 
## 0.491 0.509 
## 
## Final Objective:
##       k1 
## -2828.15 
## 
## Cluster means:
## [[1]]
##        y1        y2        y3 
## -1.968374 -1.977109 -2.022991 
## 
## [[2]]
##       y1       y2       y3 
## 2.004663 1.969927 1.976144 
## 
## Cluster covariances:
## [[1]]
##           y1        y2        y3
## y1 0.9882814 0.5352541 0.5061826
## y2 0.5352541 1.0098347 0.5847462
## y3 0.5061826 0.5847462 1.0150855
## 
## [[2]]
##           y1        y2        y3
## y1 0.8619684 0.3846493 0.3556651
## y2 0.3846493 0.9157168 0.4517638
## y3 0.3556651 0.4517638 0.8882217
## 
## Cluster responsibilities:
##             k1           k2
## 1 9.999937e-01 6.322548e-06
## 2 8.252326e-06 9.999917e-01
## 3 8.250334e-01 1.749666e-01
## 4 9.235427e-08 9.999999e-01
## 5 2.004523e-06 9.999980e-01
## 6 9.999979e-01 2.134687e-06
## 
## Cluster assignments:
##   Assignment
## 1          1
## 2          2
## 3          1
## 4          2
## 5          2
## 6          1
```

### Four Components with Missingness

In this example, $10^{3}$ observations are simulated from a four-component bivariate normal distribution with 10% missingness. Since the model has multiple components, the output is an object of class `mix`.


```r
set.seed(200);
# Four components with missingness
M = list(c(2,2),c(2,-2),c(-2,2),c(-2,-2));
S = 0.5*diag(2);
Y = rGMM(n=1000,d=2,k=4,pi=c(0.35,0.15,0.15,0.35),m=0.1,M=M,S=S);
M = fit.GMM(Y=Y,k=4,maxit=10,eps=1e-8);
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
## 10 update(s) performed without reaching tolerance limit.
## 
## Normal Mixture Model with 4 Components. 
## Cluster Proportions:
##    k1    k2    k3    k4 
## 0.170 0.303 0.171 0.356 
## 
## Final Objective:
##       k1 
## -3115.16 
## 
## Cluster means:
## [[1]]
##         y1         y2 
## -0.8214704  0.7002373 
## 
## [[2]]
##       y1       y2 
## 2.013224 2.026105 
## 
## [[3]]
##         y1         y2 
##  0.4390973 -0.5140020 
## 
## [[4]]
##        y1        y2 
## -2.063338 -2.027735 
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

