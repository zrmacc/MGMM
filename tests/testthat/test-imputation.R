library(MGMM)

test_that("Conditional distribution.", {

  # Case: 2 x 2.
  # Input data. 
  mu <- c(1, 2)
  sigma <- matrix(c(2, 1, 1, 2), nrow = 2)
  y <- c(0, 0)
  
  # Dist 1 | 2.
  idx_a <- 1
  idx_b <- 2
  dist_12 <- CalcCondDist(y, idx_a, idx_b, mu, sigma)
  
  exp_mu <- 1 - 1 * (0 - 2) / 2
  expect_equal(matrix(exp_mu, ncol = 1), dist_12$mu)
  
  exp_sigma <- 2 - 1 / 2
  expect_equal(matrix(exp_sigma, ncol = 1), dist_12$sigma)
  
  # Dist 2 | 1.
  dist_21 <- CalcCondDist(y, idx_b, idx_a, mu, sigma)
  
  exp_mu <- 2 - 1 * (0 - 1) / 2
  expect_equal(matrix(exp_mu, ncol = 1), dist_21$mu)
  
  # Note: covariance of 2 | 1 = covariance of 1 | 2.
  
  # -------------------------------------------------------
  
  # Case: 3 x 3.
  # Input data. 
  mu <- c(1, 2, 3)
  sigma <- matrix(
    c(3, 2, 1,
      2, 3, 2,
      1, 2, 3),
    ncol = 3
  )
  y <- c(0, 0, 0)
  
  # Dist 2 | (1, 3)
  idx_a <- c(2)
  idx_b <- c(1, 3)
  dist_2_13 <- CalcCondDist(y, idx_a, idx_b, mu, sigma)
  
  sub_13 <- matrix(
    c(3, 1, 
      1, 3),
    ncol = 2
  )
  exp_mu <- c(2) - c(2, 2) %*% solve(sub_13, c(0, 0) - c(1, 3))
  expect_equal(matrix(exp_mu, ncol = 1), dist_2_13$mu)
  
  exp_sigma <- c(3) - c(2, 2) %*% solve(sub_13, c(2, 2))
  expect_equal(matrix(exp_sigma, ncol = 1), dist_2_13$sigma)
  
  # Dist (1, 2) | 3
  idx_a <- c(1, 2)
  idx_b <- c(3)
  dist_12_3 <- CalcCondDist(y, idx_a, idx_b, mu, sigma)
  
  exp_mu <- c(1, 2) - matrix(c(1, 2), ncol = 1) * c(0 - 3) / 3
  expect_equal(dist_12_3$mu, exp_mu)
  
  sub_12 <- matrix(
    c(3, 2,
      2, 3),
    ncol = 2
  )
  exp_sigma <- sub_12 - matrix(c(1, 2), ncol = 1) %*% c(1, 2) / 3
  expect_equal(dist_12_3$sigma, exp_sigma)
  
})

# -----------------------------------------------------------------------------

test_that("Test imputation.", {
  withr::local_seed(101)
  
  # Checks that:
  # 1. imputed contains no NA.
  # 2. orig and imputed are identical where orig is not NA.
  is_proper_imputation <- function(orig, imputed) {
    not_na <- !is.na(orig)
    any_na <- any(is.na(imputed))
    out <- all.equal(orig[not_na], imputed[not_na]) & (!any_na)
    return(out)
  }
  
  # Case: single component, no missing data.
  data <- rGMM(n = 10, d = 3, k = 1, miss = 0)
  
  fit <- FitGMM(data, report = FALSE)
  imputed <- GenImputation(fit)
  expect_true(is_proper_imputation(data, imputed))
  
  # Case: single component, missing data.
  data[1, ] <- NA
  data[2, 1:2] <- NA
  
  fit <- FitGMM(data, report = FALSE)
  imputed <- GenImputation(fit)
  expect_true(is_proper_imputation(data, imputed))
  
  # Case: multiple components, no missing data.
  data <- rGMM(n = 20, d = 3, k = 2, miss = 0)
  
  fit <- FitGMM(data, k = 2, report = FALSE)
  imputed <- GenImputation(fit)
  expect_true(is_proper_imputation(data, imputed))
  
  # Case: multiple components, missing data.
  data[1, ] <- NA
  data[2, 1:2] <- NA
  
  fit <- FitGMM(data, k = 2, report = FALSE)
  imputed <- GenImputation(fit)
  expect_true(is_proper_imputation(data, imputed))
  
})

# -----------------------------------------------------------------------------

test_that("Test combination of multiple imputations.", {
  
  # Test 1.
  points <- list()
  covs <- list()
  
  for (i in seq_len(3)) {
    points[[i]] <- 1
    covs[[i]] <- i
  }
  
  exp_point <- 1
  exp_vcov <- mean(c(1, 2, 3))
  out <- CombineMIs(points, covs)
  expect_equal(exp_point, out$point)
  expect_equal(exp_vcov, as.numeric(out$cov))
  
  # Test 2.
  points <- list()
  covs <- list()
  
  for (i in seq_len(3)) {
    points[[i]] <- i
    covs[[i]] <- 2 * i
  }
  
  exp_point <- mean(c(1, 2, 3))
  exp_vcov <- mean(c(2, 4, 6)) + (1 + 1 / 3) * var(c(1, 2, 3))
  out <- CombineMIs(points, covs)
  expect_equal(exp_point, out$point)
  expect_equal(exp_vcov, as.numeric(out$cov))
  
})
