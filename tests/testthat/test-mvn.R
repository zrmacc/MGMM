library(MGMM)

test_that("MVN Complete Data.", {
  
  withr::local_seed(101)
  data <- rGMM(n = 1e3, d = 2, k = 1, miss = 0)
  fit <- FitMVN(data)
  expect_equal(mean(fit), c(0, 0), tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(vcov(fit), diag(2), tolerance = 0.1, ignore_attr = TRUE)
  
})

# -----------------------------------------------------------------------------

test_that("MVN Incomplete Data.", {
  
  withr::local_seed(101)
  data <- rGMM(n = 1e3, d = 2, k = 1, miss = 0.2)
  fit <- FitMVN(data, report = FALSE)
  expect_equal(mean(fit), c(0, 0), tolerance = 0.1, ignore_attr = TRUE)
  expect_equal(vcov(fit), diag(2), tolerance = 0.15, ignore_attr = TRUE)
  
})