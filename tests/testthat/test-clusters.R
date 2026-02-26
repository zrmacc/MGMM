library(MGMM)

test_that("Clustering on Complete Data.", {
  skip_on_cran()
  withr::local_seed(102)
  mu1 <- c(-4, -4)
  mu2 <- c(0, 0)
  mu3 <- c(4, 4)
  data <- rGMM(
    n = 2e2, 
    d = 2, 
    k = 3, 
    means = list(mu1, mu2, mu3),
    miss = 0
  )
  choose_k <- ChooseK(data, k0 = 2, k1 = 4, boot = 10, report = FALSE)
  expect_true(3 %in% choose_k$Choices$k_opt)
  expect_true(3 %in% choose_k$Choices$k_1se)
})

# -----------------------------------------------------------------------------

test_that("ClustQual returns expected structure.", {
  skip_on_cran()
  withr::local_seed(110)
  data <- rGMM(n = 100, d = 2, k = 2, means = list(c(-2, -2), c(2, 2)), miss = 0)
  fit <- FitGMM(data, k = 2, report = FALSE)
  qual <- ClustQual(fit)
  expect_type(qual, "list")
  expect_named(qual, c("BIC", "CHI", "DBI", "SIL"))
  expect_true(is.numeric(qual$BIC))
  expect_true(is.numeric(qual$CHI))
  expect_true(is.numeric(qual$DBI))
  expect_true(is.numeric(qual$SIL))
  expect_equal(length(qual$BIC), 1)
  expect_equal(length(qual$CHI), 1)
})

test_that("ChooseK returns Choices and Results when fits succeed.", {
  skip_on_cran()
  withr::local_seed(111)
  data <- rGMM(n = 80, d = 2, k = 2, means = list(c(-1, -1), c(1, 1)), miss = 0)
  choose_k <- ChooseK(data, k0 = 2, k1 = 3, boot = 5, report = FALSE)
  expect_named(choose_k, c("Choices", "Results"))
  expect_true(is.data.frame(choose_k$Choices))
  expect_true(all(c("Metric", "k_opt", "Metric_opt", "k_1se", "Metric_1se") %in% names(choose_k$Choices)))
  expect_true(is.data.frame(choose_k$Results))
})
