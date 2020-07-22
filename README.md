## Missingness Aware Gaussian Mixture Models

This package performs estimation and inference for Gaussian Mixture Models (GMMs) where the input data may contain missing values. Rather than imputing missing values before fitting the GMM, this package uses an extended EM algorithm to obtain the true maximum likelihood estimates of all model parameters given the observed data. In particular `MGMM` performs the following tasks:

* Maximum likelihood estimation of cluster means, covariances, and proportions.
* Calculation of cluster membership probabilities and maximum a posteriori classification of the input vectors. 
* ‘Completion’ of the input data, by imputing missing elements to their posterior means. 

## Main Functions

* `fit.GMM` estimates model parameters, performs classification and imputation.
* `rGMM` simulates observations from a GMM, potentially with missingness. 
* `ChooseK` provides guidance on choosing the number of clusters. 

## Compact Example


```r
set.seed(101)
require(MGMM)
```

```
## Loading required package: MGMM
```

```r
# Parameter settings.
mean_list <- list(
  c(1, 1),
  c(-1, -1)
)
cov_list <- list(
  matrix(c(1, -0.5, -0.5, 1), nrow = 2),
  matrix(c(1, 0.5, 0.5, 1), nrow = 2)
)

# Generate data.
data <- rGMM(
  n = 1e3, 
  d = 2, 
  k = 2, 
  miss = 0.1, 
  means = mean_list, 
  covs = cov_list
)

# ChooseK.
choose_k <- ChooseK(
  data,
  k0 = 2,
  k1 = 4,
  boot = 10,
  maxit = 10,
  eps = 1e-4,
  report = T
)
```

```
## Cluster size 2 complete. 11 fit(s) succeeded.
## Cluster size 3 complete. 11 fit(s) succeeded.
## Cluster size 4 complete. 11 fit(s) succeeded.
```

```r
# Cluster number recommendations. 
show(choose_k$Choices)
```

```
##   Metric k_opt   Metric_opt k_1se   Metric_1se
## 1    BIC     4 2246.2227138     4 2246.2227138
## 2    CHI     4    4.6410289     4    4.6410289
## 3    DBI     2    0.7876387     2    0.7876387
## 4    SIL     2    0.4762178     2    0.4762178
```

```r
# Estimation
fit <- fit.GMM(
  data,
  k = 2,
  maxit = 10
)
```

```
## Objective increment:  13.3 
## Objective increment:  5.22 
## Objective increment:  4.16 
## Objective increment:  3.27 
## Objective increment:  2.63 
## Objective increment:  2.12 
## Objective increment:  1.7 
## Objective increment:  1.37 
## Objective increment:  1.11 
## Objective increment:  0.903 
## 10 update(s) performed without reaching tolerance limit.
```

```r
# Estimated means. 
show(fit@Means)
```

```
## [[1]]
##        y1        y2 
## -1.036999 -1.052884 
## 
## [[2]]
##        y1        y2 
## 0.9515979 0.9609898
```

```r
# Estimated covariances. 
show(fit@Covariances)
```

```
## [[1]]
##           y1        y2
## y1 0.9647105 0.5370824
## y2 0.5370824 0.9778174
## 
## [[2]]
##            y1         y2
## y1  0.9986800 -0.4598704
## y2 -0.4598704  0.9707609
```

```r
# Cluster assignments. 
head(fit@Assignments)
```

```
##   Assignments      Entropy
## 1           2 9.921032e-02
## 2           1 8.258841e-01
## 2           1 5.318370e-07
## 2           1 6.179668e-03
## 1           2 9.640208e-02
## 2           1 3.354176e-04
```

```r
# Completed data. 
head(fit@Completed)
```

```
##           y1          y2
## 1  1.6512855  2.60621938
## 2 -0.5721069 -0.15672789
## 2 -2.0045376 -2.31888263
## 2 -0.6229388 -1.51543968
## 1  2.0258413  0.06921658
## 2 -1.3476380 -1.51915826
```

## Vignette

Additional examples and details may be found [here](https://github.com/zrmacc/MGMM/tree/master/vignettes/Vignette.pdf).
