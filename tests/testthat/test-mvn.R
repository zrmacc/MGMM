library(MGMM)

test_that("MVN Complete Data.", {
  withr::local_seed(101)
  data <- rGMM(n = 1e3, d = 2, k = 1, miss = 0)
  fit <- FitGMM(data, k = 1)
  expect_equal(mean(fit), c(0, 0), tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(vcov(fit), diag(2), tolerance = 0.1, ignore_attr = TRUE)
})


# -----------------------------------------------------------------------------

test_that("MVN Incomplete Data.", {
  withr::local_seed(101)
  data <- rGMM(n = 1e3, d = 2, k = 1, miss = 0.2)
  fit <- FitGMM(data, k = 1, report = FALSE)
  expect_equal(mean(fit), c(0, 0), tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(vcov(fit), diag(2), tolerance = 0.15, ignore_attr = TRUE)
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
    k = 1, 
    means = rep(2, d),
    miss = 0.1
  )
  fit <- expect_error(FitGMM(data, k = 1, lambda = 1e-1, report = FALSE), NA)
  
})

# -----------------------------------------------------------------------------

test_that("S3 methods work on mvn fit.", {
  withr::local_seed(109)
  data <- rGMM(n = 50, d = 2, k = 1, miss = 0)
  fit <- FitGMM(data, k = 1, report = FALSE)
  expect_s4_class(fit, "mvn")
  expect_equal(length(mean(fit)), 2)
  expect_equal(dim(vcov(fit)), c(2, 2))
  ll <- suppressWarnings(logLik(fit))
  expect_true(is.numeric(ll) && length(ll) == 1)
  expect_error(capture.output(print(fit)), NA)
})

