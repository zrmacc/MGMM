# Purpose: Partition data by missingness pattern.
# Updated: 2021-12-09


#' Partition Data by Missingness Pattern
#' 
#' Returns a list with the input data split in separate matrices for complete
#' cases, incomplete cases, and empty cases.
#' 
#' @param data Data.frame.
#' @return List containing:
#' \itemize{
#'   \item The original row and column names: `orig_row_names`, `orig_col_names`.
#'   \item The original row and column numbers: `n_row` and `n_col`.
#'   \item The complete cases `data_comp`.
#'   \item The incomplete cases `data_incomp`.
#'   \item The empty cases `data_empty`.
#'   \item Counts of complete `n0`, incomplete `n1`, and empty `n2` cases.
#'   \item Initial order of the observations `init_order`. 
#' }
#' @export

PartitionData <- function(data) {
  
  d <- ncol(data)
  idx <- seq(1:nrow(data))
  is_comp <- complete.cases(data)
  is_incomp <- !is_comp
  
  # Complete cases
  data_comp <- data[is_comp, , drop = FALSE]
  idx_comp <- idx[is_comp]
  
  # Incomplete cases
  data_incomp <- data[is_incomp, , drop = FALSE]
  idx_incomp <- idx[is_incomp]
  
  # Empty cases
  is_empty <- apply(data_incomp, 1, function(x){
    sum(is.na(x)) == d
  })
  data_empty <- data_incomp[is_empty, , drop = FALSE]
  idx_empty <- idx_incomp[is_empty]
  
  # Remove empty cases
  data_incomp <- data_incomp[!is_empty, , drop = FALSE]
  idx_incomp <- idx_incomp[!is_empty]
  
  # Output
  out <- list()
  out$orig_row_names <- rownames(data)
  out$orig_col_names <- colnames(data)
  
  out$n_row <- nrow(data)
  out$n_col <- ncol(data)
  
  out$n0 <- nrow(data_comp)
  out$n1 <- nrow(data_incomp)
  out$n2 <- nrow(data_empty)
  
  out$data_comp <- data_comp
  out$data_incomp <- data_incomp
  out$data_empty <- data_empty
  
  out$idx_comp <- idx_comp
  out$idx_incomp <- idx_incomp
  out$idx_empty <- idx_empty
  out$init_order <- c(idx_comp, idx_incomp, idx_empty)
  return(out)
}


#' Reconstitute Data
#' 
#' Reassembles a data matrix split by missingness pattern.
#' 
#' @param split_data Split data are returned by \code{\link{PartitionData}}.
#' @return Numeric matrix.
#' @export

ReconstituteData <- function(split_data) {
  
  d <- split_data$n_col 
  out <- rbind(
    split_data$data_comp,
    split_data$data_incomp,
    split_data$data_empty
  )
  
  # Restore initial order. 
  init_order <- split_data$init_order
  out <- out[order(init_order), , drop = FALSE]
  
  # Output
  rownames(out) <- split_data$orig_row_names
  colnames(out) <- split_data$orig_col_names  

  return(out)  
} 
