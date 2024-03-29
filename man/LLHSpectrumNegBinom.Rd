% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LLHSpectrumNegBinom.R
\name{LLHSpectrumNegBinom}
\alias{LLHSpectrumNegBinom}
\title{Likelihood that 1 observed spectrum was generated from a vector of expected
mutation counts using negative binomial distribution}
\usage{
LLHSpectrumNegBinom(spectrum, expected.counts, nbinom.size, verbose = FALSE)
}
\arguments{
\item{spectrum}{An observed spectrum (a numeric vector)}

\item{expected.counts}{A vector of (integer) expected mutation counts, one
expected count for each mutation type. We want to know the
likelihood that this model generated the observed
spectrum, assuming each mutational types generates counts according to
a negative binomial distribution with
the given \code{expected.counts} (argument \code{mu}
to \code{\link[stats]{NegBinomial}}) and dispersion parameter
\code{nbinom.size}.}

\item{nbinom.size}{The dispersion parameter for the negative
binomial distribution; smaller is more dispersed.
See \code{\link[stats]{NegBinomial}}.}

\item{verbose}{If \code{TRUE} print messages under some circumstances.}
}
\value{
\code{log(likelihood(spectrum | expected.counts))}, or,
in more detail,
the sum of the negative binomial likelihoods
that each element of the
spectrum (i.e., the count for each mutation type e.g. ACT > AAT)
was generated from the expected count for that mutation type.
}
\description{
Likelihood that 1 observed spectrum was generated from a vector of expected
mutation counts using negative binomial distribution
}
\keyword{internal}
