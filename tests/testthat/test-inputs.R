library(MGMM)

test_that("FitGMM errors on invalid inputs.", {
  withr::local_seed(112)
  data <- rGMM(n = 30, d = 2, k = 2, miss = 0)

  # Data must be matrix
  expect_error(FitGMM(as.data.frame(data), k = 2), "numeric matrix")

  # init_means wrong length
  expect_error(
    FitGMM(data, k = 2, init_means = list(c(0, 0))),
    "one is required for each"
  )

  # init_covs wrong length
  expect_error(
    FitGMM(data, k = 2, init_covs = list(diag(2))),
    "one is required for each"
  )

  # init_means wrong dimension
  expect_error(
    FitGMM(data, k = 2, init_means = list(c(0, 0, 0), c(1, 1, 1))),
    "length of ncol"
  )

  # fix_means without init_means
  expect_error(
    FitGMM(data, k = 2, fix_means = TRUE),
    "initial values are required"
  )
})

test_that("rGMM errors on invalid inputs.", {
  # miss out of range
  expect_error(rGMM(n = 10, d = 2, k = 1, miss = -0.1), "\\[0,1\\)")
  expect_error(rGMM(n = 10, d = 2, k = 1, miss = 1), "\\[0,1\\)")

  # pi length must equal k
  expect_error(rGMM(n = 10, d = 2, k = 2, pi = c(0.5)), "length k")

  # pi must be in (0,1)
  expect_error(rGMM(n = 10, d = 2, k = 2, pi = c(0, 1)), "reside in")
})

test_that("ChooseK errors when k0 < 2.", {
  data <- matrix(rnorm(20), ncol = 2)
  expect_error(ChooseK(data, k0 = 1, k1 = 3, boot = 2, report = FALSE), "At least 2 clusters")
})
