library(MGMM)

test_that("Partition data.", {
  
  data <- data.frame(rbind(c(NA, NA), c(NA, 0), c(0, NA), c(0, 0)))
  part_data <- PartitionData(data)
  
  expect_equal(part_data$n_row, 4)
  expect_equal(part_data$n_col, 2)
  
  # Complete cases.
  expect_equal(part_data$n0, 1) 
  
  # Partial cases.
  expect_equal(part_data$n1, 2)  
  
  # Empty cases.
  expect_equal(part_data$n2, 1)  
  
  # Original order.
  expect_equal(part_data$init_order, c(4, 2, 3, 1))
  
})

test_that("Reconstitute data.", {
  
  data <- rbind(c(0, 0), c(1, 1), c(2, 2))
  
  # Case of no missing data.
  split_data <- PartitionData(data)
  recov_data <- ReconstituteData(split_data)
  expect_equal(data, recov_data)
  
  # Case of incomplete data.
  data[2, 2] <- NA
  split_data <- PartitionData(data)
  recov_data <- ReconstituteData(split_data)
  expect_equal(data, recov_data)
  
  # Case of incomplete data and empty data.
  data[2, 2] <- NA
  data[3, 1] <- NA
  data[3, 2] <- NA
  split_data <- PartitionData(data)
  recov_data <- ReconstituteData(split_data)
  expect_equal(data, recov_data)
  
})
