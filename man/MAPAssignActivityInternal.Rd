% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAPAssignActivity1.R
\name{MAPAssignActivityInternal}
\alias{MAPAssignActivityInternal}
\title{Find a Maximum A Posteriori assignment of signature exposures for one spectrum.}
\usage{
MAPAssignActivityInternal(
  spect,
  sigs,
  sigs.presence.prop,
  max.level = 5,
  p.thresh = 0.05,
  m.opts = DefaultManyOpts(),
  max.mc.cores = min(20, 2^max.level),
  max.subsets = 1000,
  max.presence.proportion = 0.99,
  progress.monitor = NULL,
  seed = NULL,
  use.sparse.assign = FALSE
)
}
\arguments{
\item{spect}{A single spectrum.}

\item{sigs}{A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.}

\item{sigs.presence.prop}{The proportions of samples that contain each
signature. A numerical vector (values between 0 and 1), with names
being a subset of \code{colnames(sigs)}. See \code{\link{ExposureProportions}}
for more details.}

\item{max.level}{The maximum number of signatures to try removing.}

\item{p.thresh}{If
the p value for a better reconstruction with a
set of signatures (as opposed to
without that set of signatures)
is > than this argument, then we can use exposures without this set.}

\item{m.opts}{See \code{\link{DefaultManyOpts}}.}

\item{max.mc.cores}{The maximum number of cores to use.
On Microsoft Windows machines it is silently changed to 1.}

\item{max.subsets}{This argument provides a way to 
heuristically limit the
amount of time spent by this function. Larger values of this
argument will tend to allow longer running times.
The algorithm
successively tries to remove all subsets of 1 signature, 2
signatures, 3 signatures, etc., down to \code{max.level}.
(Not every subset is tested at each level; if a subset was
already found to be necessary the algorithm does not test
supersets of that subset.) If at any level the algorithm
needs to test more than \code{max.subsets} this function will
not proceed.}

\item{max.presence.proportion}{The maximum value of the proportion
of tumors that must have a given signature.}

\item{progress.monitor}{Function called at the start of each
new level (number of signatures to try excluding). Must
take named arguments \code{value} and \code{detail}, and
no others. Designed for a \code{\link[ipc]{AsyncProgress}}
progress bar function.}

\item{seed}{Random seed; set this to get reproducible results. (The
numerical optimization is in two phases; the first, global phase
might rarely find different optima depending on the random
seed.)}

\item{use.sparse.assign}{Whether to use sparse assignment. If \code{TRUE},
arguments designed for Maximum A Posteriori assignment such as
\code{sigs.presence.prop} will be ignored.}
}
\description{
This function also can do sparse assignment by specifying \code{use.sparse.assign = TRUE}.
}
\keyword{internal}
