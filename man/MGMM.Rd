% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_About.R
\docType{package}
\name{MGMM}
\alias{MGMM}
\title{MGMM: Missingness Aware Gaussian Mixture Models}
\description{
Parameter estimation and classification for Gaussian Mixture Models (GMMs) in
the presence of missing data. This package complements existing
implementations by allowing for both missing elements in the input vectors
and full (as opposed to strictly diagonal) covariance matrices. Estimation is
performed using an expectation conditional maximization algorithm that
accounts for missingness of both the cluster assignments and the vector
components. The output includes the marginal cluster membership
probabilities; the mean and covariance of each cluster; the posterior
probabilities of cluster membership; and a completed version of the input
data, with missing values imputed to their posterior expectations. For
additional details, please see McCaw ZR, Julienne H, Aschard H. "MGMM: an R
package for fitting Gaussian Mixture Models on Incomplete Data."
<doi:10.1101/2019.12.20.884551>.
}
\author{
Zachary R. McCaw
}
