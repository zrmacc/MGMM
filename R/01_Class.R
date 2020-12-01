# Purpose: Define custom classes.
# Updated: 20/07/16

#' Mixture Model Class
#'
#' Defines the class returned by the fitting functions.
#'
#' @slot Assignments Maximum a posteriori assignment.
#' @slot Completed Completed data, with missing values replaced by posterior expectations.
#' @slot Components Components.
#' @slot Covariances Fitted cluster covariances.
#' @slot Density Component density at observation.
#' @slot Means Fitted cluster means.
#' @slot Objective Final value of the EM objective.
#' @slot Proportions Fitted cluster proportions.
#' @slot Responsibilities Posterior membership probabilities.
#' @name mix-class
#' @rdname mix-class
#' @exportClass mix

setClass(
  Class = "mix", 
  representation = representation(
    Assignments = "matrix", 
    Completed = "matrix",
    Components = "numeric", 
    Covariances = "list", 
    Density = "matrix", 
    Means = "list", 
    Objective = "numeric",
    Proportions = "numeric", 
    Responsibilities = "matrix"
    )
  )


#' Print for Fitted Mixture Model
#'
#' Print method for objects of class \code{mix}.
#'
#' @param x A \code{mix} object.
#' @param ... Unused.
#' @export

print.mix <- function(x, ...) {
  
  # Components
  k <- x@Components
  
  # Parameters
  pi <- signif(x@Proportions, digits = 3)
  q <- signif(x@Objective)
  
  # Display
  cat(paste0("Normal Mixture Model with ", k, " Components."), "\n")
  
  cat("Cluster Proportions:\n")
  print(pi)
  cat("\n")
  
  cat("Final Objective:\n")
  print(q)
  cat("\n")
}


#' Show for Fitted Mixture Models
#' @param object A \code{mix} object.
#' @rdname mix-method
#' @importFrom methods show

setMethod(
  f = "show", 
  signature = c(object = "mix"), 
  definition = function(object) {print.mix(x = object)}
  )

