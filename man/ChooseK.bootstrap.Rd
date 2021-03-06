% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_Clustering.R
\name{ChooseK.bootstrap}
\alias{ChooseK.bootstrap}
\title{Bootstrap Quality Metrics.}
\usage{
ChooseK.bootstrap(
  boot,
  data,
  k,
  init_means,
  fix_means,
  init_covs,
  init_props,
  maxit,
  eps
)
}
\arguments{
\item{boot}{Bootstrap replicates.}

\item{data}{Numeric data matrix.}

\item{k}{Number of clusters.}

\item{init_means}{Optional list of initial mean vectors.}

\item{fix_means}{Fix the means to their starting value? Must initialize.}

\item{init_covs}{Optional list of initial covariance matrices.}

\item{init_props}{Optional vector of initial cluster proportions.}

\item{maxit}{Maximum number of EM iterations.}

\item{eps}{Minimum acceptable increment in the EM objective.}
}
\value{
Numeric matrix of clustering metrics. Returns null if the models
  fails to fit the observed data.
}
\description{
Bootstrap Quality Metrics.
}
