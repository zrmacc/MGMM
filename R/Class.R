########################
# For Mixtures
########################

#' Mixture Model Class
#'
#' Defines the class returned by functions that 
#' fit mixture models. 
#'
#' @slot Components Components.
#' @slot Means Fitted cluster means.
#' @slot Covariances Fitted cluster covariances.
#' @slot Proportions Fitted cluster proportions.
#' @slot Objective Final value of the EM objective. 
#' @slot Responsibilities Final cluster responsibilities.
#' @slot Assignments Final cluster assignments.
#' @name mix-class
#' @rdname mix-class
#' @exportClass mix

setClass(Class="mix",representation=representation(Components="numeric",Means="list",Covariances="list",Proportions="numeric",
                                                   Objective="numeric",Responsibilities="data.frame",Assignments="data.frame"));

#' Print for Fitted Mixture Model
#'
#' Print method for objects of class \code{mix}.
#'
#' @param x A \code{mix} object.
#' @param ... Unused.
#' @export

print.mix = function(x,...){
  # Components
  k = x@Components;
  # Parameters
  pi = signif(x@Proportions,digits=3);
  q = signif(x@Objective);
  # Display
  cat(paste0("Normal Mixture Model with ",k," Components."),"\n");
  cat("Cluster Proportions:\n");
  print(pi);
  cat("\n");
  cat("Final Objective:\n");
  print(q);
  cat("\n");
}

#' Show for Fitted Survival Models
#' @param object A \code{mix} object.
#' @rdname mix-method
#' @importFrom methods show

setMethod(f="show",signature=c(object="mix"),definition=function(object){print.mix(x=object)});
