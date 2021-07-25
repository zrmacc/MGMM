library(MGMM)

# -----------------------------------------------------------------------------

test_that("Density evaluation of incomplete obs.", {
  
  means <- list(c(1, 0))
  covs <- list(diag(2))
  g <- function(y) {EvalDensIncompObs(y, means, covs, 1)}
  
  # No elements missing.
  y <- c(0, 0)
  obs <- g(y)
  exp <- prod(stats::dnorm(x = y, mean = c(1, 0)))
  expect_equal(obs, exp)
  
  # First element missing.
  y <- c(NA, 0)
  obs <- g(y)
  exp <- stats::dnorm(x = 0, mean = 0, sd = 1)
  expect_equal(obs, exp)
  
  # Second element missing.
  y <- c(0, NA)
  obs <- g(y)
  exp <- stats::dnorm(x = 0, mean = 1, sd = 1)
  expect_equal(obs, exp)
  
})

# -----------------------------------------------------------------------------

test_that("Working response vector.", {
  
  # Mean and covariance.
  mu <- c(0, 0)
  sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  g <- function(y) {WorkRespIndiv(y, mu, sigma, 1)}
  
  # No missing elements.
  y <- c(0, 0)
  obs <- g(y)
  exp <- y
  expect_equal(obs, exp)
  
  # First element missing.
  y <- c(NA, 2)
  obs <- g(y)
  exp <- c(mu[1] + sigma[1, 2] / sigma[2, 2] * (y[2] - mu[2]), y[2])
  expect_equal(obs, exp)
  
  # Second element missing.
  y <- c(1, NA)
  obs <- g(y)
  exp <- c(y[1], mu[2] + sigma[2, 1] / sigma[1, 1] * (y[1] - mu[1]))
  expect_equal(obs, exp)

})

# -----------------------------------------------------------------------------

test_that("Working residual outer product.", {
  
  new_mean <- c(1, 1)
  old_mean <- c(0, 0)
  old_cov <- matrix(c(2, 0.5, 0.5, 2), nrow = 2)
  g <- function(data) {
    ExpResidOP(data, new_mean, old_mean, old_cov)
  }
  
  # First element missing.
  data <- matrix(c(NA, 0), nrow = 1)
  obs <- g(data)
  res <- WorkRespIndiv(data, old_mean, old_cov, 1) - new_mean
  correction <- old_cov[1, 1] - old_cov[1, 2] / old_cov[2, 2] * old_cov[2, 1]
  exp <- res %*% t(res) + diag(c(correction, 0))
  expect_equal(obs, exp)
  
  # Second element missing.
  data <- matrix(c(1, NA), nrow = 1)
  obs <- g(data)
  res <- WorkRespIndiv(data, old_mean, old_cov, 1) - new_mean
  correction <- old_cov[2, 2] - old_cov[2, 1] / old_cov[1, 1] * old_cov[1, 2]
  exp <- res %*% t(res) + diag(c(0, correction))
  expect_equal(obs, exp)
  
})
