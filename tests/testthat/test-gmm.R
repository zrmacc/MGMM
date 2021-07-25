library(MGMM)

test_that("GMM Complete Data.", {
  skip_on_cran()
  withr::local_seed(102)
  mu1 <- c(-2, -2)
  mu2 <- c(2, 2)
  pi <- c(0.7, 0.3)
  data <- rGMM(
    n = 2e3, 
    d = 2, 
    k = 2, 
    means = list(mu1, mu2),
    miss = 0,
    pi = pi
  )
  fit <- FitMix(data, maxit = 20)
  mu <- mean(fit)
  
  expect_equal(mu[[2]], mu1, tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(mu[[1]], mu2, tolerance = 0.1, ignore_attr = TRUE)
  
  sigma <- vcov(fit)
  expect_equal(sigma[[1]], diag(2), tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(sigma[[2]], diag(2), tolerance = 0.1, ignore_attr = TRUE)
  
  prop <- fit@Proportions
  expect_equal(prop, rev(pi), tolerance = 0.1, ignore_attr = TRUE)
  
})

# -----------------------------------------------------------------------------


test_that("GMM Incomplete Data.", {
  skip_on_cran()
  withr::local_seed(102)
  mu1 <- c(-2, -2)
  mu2 <- c(2, 2)
  pi <- c(0.4, 0.6)
  data <- rGMM(
    n = 2e3, 
    d = 2, 
    k = 2, 
    means = list(mu1, mu2),
    miss = 0.1,
    pi = pi
  )
  fit <- FitMix(data, maxit = 20)
  mu <- mean(fit)
  
  expect_equal(mu[[1]], mu1, tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(mu[[2]], mu2, tolerance = 0.1, ignore_attr = TRUE)
  
  sigma <- vcov(fit)
  expect_equal(sigma[[1]], diag(2), tolerance = 0.2, ignore_attr = TRUE)
  expect_equal(sigma[[2]], diag(2), tolerance = 0.2, ignore_attr = TRUE)
  
  prop <- fit@Proportions
  expect_equal(prop, pi, tolerance = 0.1, ignore_attr = TRUE)
})