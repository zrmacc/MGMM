# Purpose: Define custom classes.
# Updated: 2021-07-24

# -----------------------------------------------------------------------------

#' Multivariate Normal Model Class
#'
#' Defines a class to hold multivariate normal models.
#'
#' @slot Completed Completed data, with missing values imputed to their
#'   posterior expectations.
#' @slot Covariance Fitted covariance matrix.
#' @slot Data Original data, with missing values present.
#' @slot Mean Fitted mean vector.
#' @slot Objective Final value of the EM objective.
#' @name mvn-class
#' @rdname mvn-class
#' @exportClass mvn

methods::setClass(
  Class = "mvn", 
  representation = representation(
    Completed = "matrix",
    Covariance = "matrix",
    Data = "matrix",
    Mean = "vector", 
    Objective = "numeric"
  )
)


#' Log likelihood for Fitted MVN Model
#'
#' @param object A \code{mvn} object.
#' @param ... Unused.
#' @export

logLik.mvn <- function(object, ...) {
  warning("Returning EM objective rather than a true log likelihood.\n")
  return(object@Objective)
}


#' Mean for Fitted MVN Model
#'
#' @param x A \code{mvn} object.
#' @param ... Unused.
#' @export

mean.mvn <- function(x, ...) {
  return(x@Mean)
}


#' Print for Fitted MVN Model
#'
#' @param x A \code{mvn} object.
#' @param ... Unused.
#' @export

print.mvn <- function(x, ...) {
  
  # Parameters.
  mu <- signif(x@Mean, digits = 3)
  sigma <- signif(x@Covariance, digits = 3)
  obj <- signif(x@Objective, digits = 3)
  
  # Display.
  cat(paste0("Multivariate Normal Model."), "\n\n")
  
  cat("Estimated mean:\n")
  print(mu)
  cat("\n")
  
  cat("Estimated covariance:\n")
  print(sigma)
  cat("\n")
  
  cat("Final Objective:\n")
  print(obj)
  cat("\n")
}


#' Covariance for Fitted MVN Model
#'
#' @param object A \code{mvn} object.
#' @param ... Unused.
#' @export

vcov.mvn <- function(object, ...) {
  return(object@Covariance)
}


#' Show for Multivariate Normal Models
#' @param object A \code{mvn} object.
#' @rdname mvn-method

methods::setMethod(
  f = "show", 
  signature = c(object = "mvn"), 
  definition = function(object) {print.mvn(x = object)}
)


# -----------------------------------------------------------------------------

#' Mixture Model Class
#'
#' Defines a class to hold Gaussian Mixture Models.
#'
#' @slot Assignments Maximum a posteriori assignments.
#' @slot Completed Completed data, with missing values imputed to their
#'   posterior expectations.
#' @slot Components Number of components.
#' @slot Covariances List of fitted cluster covariance matrices.
#' @slot Data Original data, with missing values present.
#' @slot Density Density of each component at each example.
#' @slot Means List of fitted cluster means.
#' @slot Objective Final value of the EM objective.
#' @slot Proportions Fitted cluster proportions.
#' @slot Responsibilities Posterior membership probabilities for each example.
#' @name mix-class
#' @rdname mix-class
#' @exportClass mix

methods::setClass(
  Class = "mix", 
  representation = representation(
    Assignments = "matrix", 
    Completed = "matrix",
    Components = "numeric", 
    Covariances = "list",
    Data = "matrix",
    Density = "matrix", 
    Means = "list", 
    Objective = "numeric",
    Proportions = "numeric", 
    Responsibilities = "matrix"
    )
  )


#' Log likelihood for Fitted GMM
#'
#' @param object A \code{mix} object.
#' @param ... Unused.
#' @export

logLik.mix <- function(object, ...) {
  warning("Returning EM objective rather than a true log likelihood.\n")
  return(object@Objective)
}


#' Mean for Fitted GMM
#'
#' @param x A \code{mix} object.
#' @param ... Unused.
#' @export

mean.mix <- function(x, ...) {
  return(x@Means)
}

#' Print for Fitted GMM
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
  cat(paste0("Gaussian Mixture Model with ", k, " Components."), "\n\n")
  
  cat("Cluster Proportions:\n")
  print(pi)
  cat("\n")
  
  cat("Final Objective:\n")
  print(q)
  cat("\n")
}


#' Covariance for Fitted GMM
#'
#' @param object A \code{mix} object.
#' @param ... Unused.
#' @export

vcov.mix <- function(object, ...) {
  return(object@Covariances)
}


#' Show for Fitted Mixture Models
#' @param object A \code{mix} object.
#' @rdname mix-method

methods::setMethod(
  f = "show", 
  signature = c(object = "mix"), 
  definition = function(object) {print.mix(x = object)}
)

