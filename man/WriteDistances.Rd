% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShowSigActivity.R
\name{WriteDistances}
\alias{WriteDistances}
\title{Write distances from the output generated by \code{AddSigActivity}}
\usage{
WriteDistances(list.of.sig.activity, output.dir, base.filename, mutation.type)
}
\arguments{
\item{list.of.sig.activity}{A list of contributing signature activity
information for multiple spectra. See the return value of
\code{\link{AddSigActivity}} for more details.}

\item{output.dir}{The directory to save the results. Create
this directory if it does not exist.}

\item{base.filename}{Optional. \code{base.filename} will be appended to the start
of the names of files generated inside \code{output.dir}.}

\item{mutation.type}{The mutation type of the mutational spectrum.}
}
\description{
Write distances from the output generated by \code{AddSigActivity}
}
\keyword{internal}
