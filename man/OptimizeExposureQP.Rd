% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptimizeExposureQP.R
\name{OptimizeExposureQP}
\alias{OptimizeExposureQP}
\title{Quadratic programming optimization of signature activities}
\usage{
OptimizeExposureQP(spectrum, signatures)
}
\arguments{
\item{spectrum}{Mutational signature spectrum
as a numeric vector or single column data frame or matrix.}

\item{signatures}{Matrix or data frame of signatures from which to
reconstruct \code{spectrum}. Rows are mutation types and columns are
signatures. Should have column names for interpretable results. Cannot be a
vector because the column names are needed.}
}
\value{
A vector of exposures with names being the colnames from
  \code{signatures}.
Code adapted from \code{SignatureEstimation::decomposeQP}.
}
\description{
Quadratic programming optimization of signature activities
}
\keyword{internal}
