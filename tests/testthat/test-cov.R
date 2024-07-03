library(MGMM)

test_that("Covariance calculation.", {
  withr::local_seed(1010)
  n <- 10
  d <- 3
  a <- matrix(stats::rnorm(n * d), nrow = n)  
  
  # Covariance calculation.
  exp <- stats::cov(a, a)
  obs <- matCov(a, a)
  expect_equal(exp, obs, tolerance = 1e-6)
  
  # Covariance with ridge.
  eps <- 0.001
  exp <- stats::cov(a, a) + eps * diag(d) / (n - 1)
  obs <- matCov(a, a, eps = eps)
  expect_equal(exp, obs, tolerance = 1e-6)
  
  # Calculate correlation
  exp <- stats::cor(a, a)
  obs <- matCov(a, a, corMat = TRUE)
  expect_equal(exp, obs, tolerance = 1e-6)
  
})

