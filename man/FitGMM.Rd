% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_Estimation.R
\name{FitGMM}
\alias{FitGMM}
\title{Estimate Multivariate Normal Mixture}
\usage{
FitGMM(
  data,
  k = 1,
  init_means = NULL,
  fix_means = FALSE,
  init_covs = NULL,
  lambda = 0,
  init_props = NULL,
  maxit = 100,
  eps = 1e-06,
  report = TRUE
)
}
\arguments{
\item{data}{Numeric data matrix.}

\item{k}{Number of mixture components. Defaults to 1.}

\item{init_means}{Optional list of initial mean vectors.}

\item{fix_means}{Fix the means to their starting value? Must provide initial
values.}

\item{init_covs}{Optional list of initial covariance matrices.}

\item{lambda}{Optional ridge term added to covariance matrix to ensure 
positive definiteness.}

\item{init_props}{Optional vector of initial cluster proportions.}

\item{maxit}{Maximum number of EM iterations.}

\item{eps}{Minimum acceptable increment in the EM objective.}

\item{report}{Report fitting progress?}
}
\value{
\itemize{
  \item For a single component, an object of class \code{mvn}, containing
  the estimated mean and covariance, the final objective function, and the
  imputed data. 
  \item For a multicomponent model \eqn{k>1}, an object of class \code{mix}, 
  containing the estimated means, covariances, cluster proportions, cluster
  responsibilities, and observation assignments.
}
}
\description{
Given an \eqn{n \times d} matrix of random vectors, estimates the parameters
of a Gaussian Mixture Model (GMM). Accommodates arbitrary patterns of missingness
at random (MAR) in the input vectors.
}
\details{
Initial values for the cluster means, covariances, and proportions are
specified using \code{M0}, \code{S0}, and \code{pi0}, respectively. If the
data contains complete observations, i.e. observations with no missing
elements, then \code{fit.GMM} will attempt to initialize these parameters
internally using K-means. If the data contains no complete observations, then
initial values are required for \code{M0}, \code{S0}, and \code{pi0}.
}
\examples{
\donttest{
# Single component without missingness
# Bivariate normal observations
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
data <- rGMM(n = 1e3, d = 2, k = 1, means = c(2, 2), covs = sigma)
fit <- FitGMM(data, k = 1)

# Single component with missingness
# Trivariate normal observations
mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = sigma)
fit <- FitGMM(data, k = 2)

# Two components without missingness
# Trivariate normal observations
mean_list <- list(c(-2, -2, -2), c(2, 2, 2))
sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
data <- rGMM(n = 1e3, d = 3, k = 2, means = mean_list, covs = sigma)
fit <- FitGMM(data, k = 2)

# Four components with missingness
# Bivariate normal observations
# Note: Fitting is slow.
mean_list <- list(c(2, 2), c(2, -2), c(-2, 2), c(-2, -2))
sigma <- 0.5 * diag(2)
data <- rGMM(
n = 1000, 
d = 2, 
k = 4, 
pi = c(0.35, 0.15, 0.15, 0.35), 
m = 0.1, 
means = mean_list, 
covs = sigma)
fit <- FitGMM(data, k = 4)
}
}
\seealso{
See \code{\link{rGMM}} for data generation, and \code{\link{ChooseK}} for selecting
the number of clusters.
}
