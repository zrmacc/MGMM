# Missingness Aware Gaussian Mixture Models

Zachary McCaw <br>
Updated: 2020-12-01

This package performs estimation and inference for Gaussian Mixture Models (GMMs) where the input data may contain missing values. Rather than imputing missing values before fitting the GMM, this package uses an extended EM algorithm to obtain the true maximum likelihood estimates of all model parameters given the observed data. In particular `MGMM` performs the following tasks:

* Maximum likelihood estimation of cluster means, covariances, and proportions.
* Calculation of cluster membership probabilities and maximum a posteriori classification of the input vectors. 
* Completion of the input data, by imputing missing elements to their posterior means. 

The method is detailed in [MGMM: An R Package for fitting Gaussian Mixture Models on Incomplete Data](https://www.biorxiv.org/content/10.1101/2019.12.20.884551v2).

## Main Functions

* `fit.GMM` estimates model parameters, performs classification and imputation.
* `rGMM` simulates observations from a GMM, potentially with missingness. 
* `ChooseK` provides guidance on choosing the number of clusters. 

## Compact Example


```r
set.seed(101)
require(MGMM)

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
## 1    BIC     4 2285.4269195     4 2285.4269195
## 2    CHI     4    4.5081110     4    4.5081110
## 3    DBI     2    0.7818011     2    0.7818011
## 4    SIL     2    0.4785951     2    0.4785951
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
## Objective increment:  11.1 
## Objective increment:  2.39 
## Objective increment:  1.57 
## Objective increment:  1.08 
## Objective increment:  0.91 
## Objective increment:  0.745 
## Objective increment:  0.617 
## Objective increment:  0.507 
## Objective increment:  0.416 
## Objective increment:  0.34 
## 10 update(s) performed without reaching tolerance limit.
```

```r
# Estimated means. 
show(fit@Means)
```

```
## [[1]]
##        y1        y2 
## -1.056071 -1.070412 
## 
## [[2]]
##        y1        y2 
## 0.9437473 0.9513797
```

```r
# Estimated covariances. 
show(fit@Covariances)
```

```
## [[1]]
##           y1        y2
## y1 0.9447484 0.5201638
## y2 0.5201638 0.9611714
## 
## [[2]]
##            y1         y2
## y1  0.9973684 -0.4489898
## y2 -0.4489898  0.9728258
```

```r
# Cluster assignments. 
head(fit@Assignments)
```

```
##   Assignments      Entropy
## 1           2 7.957793e-02
## 2           1 8.426200e-01
## 2           1 8.629837e-07
## 2           1 7.790894e-03
## 1           2 8.451183e-02
## 2           1 4.521627e-04
```

```r
# Completed data. 
head(fit@Completed)
```

```
##           y1          y2
## 1  1.6512855  2.60621938
## 2 -0.5721069 -0.14379539
## 2 -2.0045376 -2.31888263
## 2 -0.6229388 -1.51543968
## 1  2.0258413  0.06921658
## 2 -1.3476380 -1.51915826
```

## Vignette

Additional examples and details may be found [here](https://github.com/zrmacc/MGMM/tree/master/vignettes/Vignette.pdf).
