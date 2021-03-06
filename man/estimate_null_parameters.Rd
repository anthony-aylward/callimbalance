% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callimbalance.R
\name{estimate_null_parameters}
\alias{estimate_null_parameters}
\title{Estimate Null Parameters}
\usage{
estimate_null_parameters(data_frame, minimum_coverage = 10,
  n_breaks = 11, spline_order = 4, cores = 1)
}
\arguments{
\item{data_frame}{a data frame containing allele counts}

\item{minimum_coverage}{filter variants to those with at least this coverage}

\item{n_breaks}{n_breaks for npbin}

\item{spline_order}{spline order for npbin}

\item{cores}{number of cores to use}
}
\value{
list of null parameters
}
\description{
estimate the shape parameters of a beta-binomial null dist
}
