library(MGMM)

test_that("rGMM returns correct dimensions and structure.", {
  withr::local_seed(101)
  n <- 100
  d <- 3
  k <- 2
  data <- rGMM(n = n, d = d, k = k, miss = 0)
  expect_equal(nrow(data), n)
  expect_equal(ncol(data), d)
  expect_equal(colnames(data), paste0("y", seq_len(d)))
  # Row names encode true cluster (1 or 2)
  expect_true(all(rownames(data) %in% c("1", "2")))
})

test_that("rGMM k=1 has single cluster and no NA when miss=0.", {
  withr::local_seed(102)
  data <- rGMM(n = 50, d = 2, k = 1, miss = 0, means = c(1, 2))
  expect_equal(nrow(data), 50)
  expect_equal(ncol(data), 2)
  expect_true(all(rownames(data) == "1"))
  expect_false(any(is.na(data)))
})

test_that("rGMM with miss > 0 introduces approximately correct proportion of NAs.", {
  withr::local_seed(103)
  n <- 500
  d <- 4
  data <- rGMM(n = n, d = d, k = 2, miss = 0.2)
  prop_miss <- mean(is.na(data))
  expect_true(prop_miss >= 0.15 && prop_miss <= 0.25)
})

test_that("rGMM respects custom pi.", {
  withr::local_seed(104)
  n <- 1000
  data <- rGMM(n = n, d = 2, k = 2, pi = c(0.8, 0.2), miss = 0,
               means = list(c(0, 0), c(10, 10)))
  tab <- table(rownames(data))
  expect_true(tab["1"] / n > 0.7 && tab["1"] / n < 0.9)
  expect_true(tab["2"] / n > 0.1 && tab["2"] / n < 0.3)
})

test_that("rGMM with custom means and covs runs and has correct dimension.", {
  withr::local_seed(105)
  mu1 <- c(-1, 1)
  mu2 <- c(1, -1)
  cov1 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  cov2 <- 0.5 * diag(2)
  data <- rGMM(n = 50, d = 2, k = 2, means = list(mu1, mu2), covs = list(cov1, cov2))
  expect_equal(nrow(data), 50)
  expect_equal(ncol(data), 2)
})
