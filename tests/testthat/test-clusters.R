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
  choose_k <- ChooseK(data, k0 = 2, k1 = 4, report = FALSE)
  expect_true(3 %in% choose_k$Choices$k_opt)
  expect_true(3 %in% choose_k$Choices$k_1se)
})
