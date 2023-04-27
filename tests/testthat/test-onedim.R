library(MGMM)
library(testthat)

test_that("Clustering on 1d data, 1 cluster.", {
  
  # No missingness.
  withr::local_seed(101)
  data <- rGMM(n = 1e2, d = 1, k = 1, miss = 0)
  expect_error(FitGMM(data = data, report = FALSE), NA)
  
  # With missingness.
  data <- rGMM(n = 3e2, d = 1, k = 1, miss = 0.2)
  expect_error({fit <- FitGMM(data = data, report = FALSE)}, NA)
  
  # Generate imputation.
  expect_error({imp <- GenImputation(fit)}, NA)
  
})


test_that("Clustering on 1d data, multiple clusters.", {
  
  # No missingness.
  withr::local_seed(101)
  data <- rGMM(n = 2e1, d = 1, k = 2, means = list(-2, 2), miss = 0)
  expect_error(FitGMM(data = data, k = 3, report = FALSE), NA)
  
  # With missingness.
  withr::local_seed(101)
  data <- rGMM(n = 2e1, d = 1, k = 2, means = list(-2, 2), miss = 0.2)
  expect_error({fit <- FitGMM(data = data, k = 3, report = FALSE)}, NA)
  
  # Generate imputation.
  expect_error({imp <- GenImputation(fit)}, NA)
  
})
