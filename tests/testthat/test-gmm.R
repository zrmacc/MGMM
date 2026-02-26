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
  fit <- FitGMM(data, k = 2, maxit = 20)
  mu <- mean(fit)
  # Component order may differ from data generation
  expect_equal(sort(sapply(mu, function(x) x[1])), sort(c(mu1[1], mu2[1])), tolerance = 0.1, ignore_attr = TRUE)
  sigma <- vcov(fit)
  expect_equal(length(sigma), 2)
  prop <- fit@Proportions
  expect_equal(sum(prop), 1, tolerance = 1e-6)
  expect_true(all(prop > 0))
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
  fit <- FitGMM(data, k = 2, maxit = 20)
  mu <- mean(fit)
  expect_equal(sort(sapply(mu, function(x) x[1])), sort(c(mu1[1], mu2[1])), tolerance = 0.1, ignore_attr = TRUE)
  sigma <- vcov(fit)
  expect_equal(length(sigma), 2)
  prop <- fit@Proportions
  expect_equal(sum(prop), 1, tolerance = 1e-6)
  expect_equal(sort(prop), sort(pi), tolerance = 0.15, ignore_attr = TRUE)
})


# -----------------------------------------------------------------------------

test_that("Rank deficient covariance matrix.", {
  skip_on_cran()
  withr::local_seed(101)
  
  d <- 10
  n <- 9
  data <- rGMM(
    n = n, 
    d = d, 
    k = 2, 
    means = list(rep(2, d), rep(-2, d)),
    miss = 0.1
  )
  fit <- expect_error(FitGMM(data, k = 1, lambda = 1e-1, report = FALSE), NA)
  
})

# -----------------------------------------------------------------------------

test_that("FitGMM with k=1 returns mvn object.", {
  withr::local_seed(106)
  data <- rGMM(n = 100, d = 2, k = 1, miss = 0)
  fit <- FitGMM(data, k = 1, report = FALSE)
  expect_s4_class(fit, "mvn")
  expect_true(methods::is(fit, "mvn"))
  expect_equal(length(fit@Mean), 2)
  expect_equal(dim(fit@Covariance), c(2, 2))
  expect_equal(dim(fit@Completed), c(100, 2))
})

test_that("S3 methods work on mix fit.", {
  withr::local_seed(107)
  data <- rGMM(n = 80, d = 2, k = 2, means = list(c(-1, -1), c(1, 1)), miss = 0)
  fit <- FitGMM(data, k = 2, report = FALSE)
  expect_s4_class(fit, "mix")
  expect_equal(length(mean(fit)), 2)
  expect_equal(length(vcov(fit)), 2)
  ll <- suppressWarnings(logLik(fit))
  expect_true(is.numeric(ll) && length(ll) == 1)
  expect_error(capture.output(print(fit)), NA)
})

test_that("FitGMM with fix_means keeps initial means.", {
  skip_on_cran()
  withr::local_seed(108)
  data <- rGMM(n = 200, d = 2, k = 2, means = list(c(-2, 0), c(2, 0)), miss = 0.05)
  init_means <- list(c(-2, 0), c(2, 0))
  fit <- FitGMM(data, k = 2, init_means = init_means, fix_means = TRUE, maxit = 15, report = FALSE)
  expect_equal(fit@Means[[1]], init_means[[1]], tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(fit@Means[[2]], init_means[[2]], tolerance = 1e-6, ignore_attr = TRUE)
})

