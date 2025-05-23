% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WRMT.R
\name{adaptive.threshold.weights}
\alias{adaptive.threshold.weights}
\title{Obtain adaptive thresholds with universal caliper c and weights w
#' @description
This function calculates adaptive thresholds for WRAT analysis with universal
caliper c and weights w. To be called by function \code{WRAT()}.}
\usage{
adaptive.threshold.weights(
  pair.data,
  num.outcome,
  weight,
  caliper = 0.2,
  clinical.threshold = NULL
)
}
\arguments{
\item{pair.data}{A data frame for pairwise comparison,
returned by \code{data.to.pair()}.}

\item{num.outcome}{A numeric value of the number of outcomes.}

\item{weight}{A numeric vector indicates the weights w for adaptive thresholds.}

\item{caliper}{A numeric value indicates the universal caliper c for adaptive thresholds.}

\item{clinical.threshold}{A numeric vector represents the minimal clinical thresholds.
Default to 0 for all endpoints.}
}
\value{
A list of outputs:
\describe{
\item{threshold}{A numeric vector of obtained adaptive thresholds.}
\item{stage}{A numeric vector indicating outcome for each comparison stage.}
}
}
\description{
Obtain adaptive thresholds with universal caliper c and weights w
#' @description
This function calculates adaptive thresholds for WRAT analysis with universal
caliper c and weights w. To be called by function \code{WRAT()}.
}
