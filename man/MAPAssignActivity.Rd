% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAPAssignActivity.R
\name{MAPAssignActivity}
\alias{MAPAssignActivity}
\title{Find Maximum A Posteriori (MAP) assignment of signature exposures that
explain multiple spectra}
\usage{
MAPAssignActivity(
  spectra,
  sigs,
  sigs.presence.prop,
  output.dir,
  max.level = 5,
  p.thresh = 0.05,
  m.opts = DefaultManyOpts(),
  num.parallel.samples = 5,
  mc.cores.per.sample = min(20, 2^max.level),
  progress.monitor = NULL,
  seed = NULL,
  max.subsets = 1000,
  use.sparse.assign = FALSE,
  drop.low.mut.samples = TRUE
)
}
\arguments{
\item{spectra}{The spectra (multiple spectra) to be reconstructed.}

\item{sigs}{A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.}

\item{sigs.presence.prop}{The proportions of samples that contain each
signature. A numerical vector (values between 0 and 1), with names
being a subset of \code{colnames(sigs)}. See \code{\link{ExposureProportions}}
for more details.}

\item{output.dir}{Directory path to save the output file.}

\item{max.level}{The maximum number of signatures to try removing.}

\item{p.thresh}{If
the p value for a better reconstruction with a
set of signatures (as opposed to
without that set of signatures)
is > than this argument, then we can use exposures without this set.}

\item{m.opts}{See \code{\link{DefaultManyOpts}}.}

\item{num.parallel.samples}{The (maximum) number of samples to run in
parallel. On Microsoft Windows machines it is silently changed to 1. Each
sample in turn can require multiple cores, as governed by
\code{mc.cores.per.sample}.}

\item{mc.cores.per.sample}{The maximum number of cores to use for each
sample. On Microsoft Windows machines it is silently changed to 1.}

\item{progress.monitor}{Function called at the start of each
new level (number of signatures to try excluding). Must
take named arguments \code{value} and \code{detail}, and
no others. Designed for a \code{\link[ipc]{AsyncProgress}}
progress bar function.}

\item{seed}{Random seed; set this to get reproducible results. (The
numerical optimization is in two phases; the first, global phase
might rarely find different optima depending on the random
seed.)}

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

\item{use.sparse.assign}{Whether to use sparse assignment. If \code{TRUE},
arguments designed for Maximum A Posteriori assignment such as
\code{sigs.presence.prop} will be ignored.}

\item{drop.low.mut.samples}{Whether to exclude low mutation samples from
the analysis. If \code{TRUE(default)}, samples with SBS total mutations less
than 100, DBS or ID total mutations less than 25 will be dropped.}
}
\value{
A list with the elements:
\itemize{
\item \code{proposed.assignment}: Proposed signature assignment for \code{spectra}
with the highest MAP found. If \code{use.sparse.assign = TRUE}, this will
be the most sparse set of signatures that can plausibly explain \code{spectra}.
\item \code{proposed.reconstruction}: Proposed reconstruction of \code{spectra} based on \code{MAP}.
If \code{use.sparse.assign = TRUE}, this will be the reconstruction based on
sparse assignment.
\item \code{reconstruction.distances}: Various distances and similarities
between \code{spectra} and \code{proposed.reconstruction}.
\item \code{all.tested}: All tested possible ways to reconstruct each
sample in \code{spectra}.
\item \code{alt.solutions}: A \code{tibble} showing all the alternative solutions
that are statistically as good as the \code{proposed.assignment} that can
plausibly reconstruct \code{spectra}.
\item \code{time.for.assignment}: Value from \code{system.time} for running
\code{MAPAssignActivity1} for each sample in \code{spectra}.
\item \code{error.messages}: Only appearing if there are errors running
\code{MAPAssignActivity}.
}

The elements \code{proposed.assignment}, \code{proposed.reconstruction},
\code{reconstruction.distances}, \code{all.tested},
\code{time.for.assignment} will be \code{NULL} if the algorithm could not
find the optimal reconstruction or there are errors coming out for
\strong{all} samples.
}
\description{
This function also can do sparse assignment by specifying \code{use.sparse.assign = TRUE}.
}
\keyword{internal}
