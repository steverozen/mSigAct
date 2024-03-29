% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AddSigActivity.R
\name{AddSigActivity}
\alias{AddSigActivity}
\title{Add contributing signature activity information for multiple spectra}
\usage{
AddSigActivity(
  spectra,
  exposure,
  sigs,
  sigs.presence.prop,
  nbinom.size = 5,
  likelihood.dist = "multinom",
  use.sparse.assign = FALSE
)
}
\arguments{
\item{spectra}{The spectra (multiple spectra) to be reconstructed.}

\item{exposure}{Exposures as a numerical matrix (or data.frame) with
signatures in rows and samples in columns. Rownames are taken as the
signature names and column names are taken as the sample IDs.}

\item{sigs}{The signatures with which we are trying to reconstruct \code{spectra}.
A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog. The column
names of \code{sigs} should be a superset of row names of \code{exposure}.}

\item{sigs.presence.prop}{The proportions of samples that contain each
signature. A numerical vector (values between 0 and 1), with names
being a subset of \code{colnames(sigs)}. See \code{\link{ExposureProportions}}
for more details.}

\item{nbinom.size}{The dispersion parameter for the negative binomial
distribution; smaller is more dispersed. See
\code{\link[stats]{NegBinomial}}.}

\item{likelihood.dist}{The probability distribution used to calculate the
likelihood, can be either "multinom" (multinomial distribution) or
"neg.binom" (negative binomial distribution).}

\item{use.sparse.assign}{Whether to use sparse assignment. If \code{TRUE},
arguments designed for Maximum A Posteriori assignment such as
\code{sigs.presence.prop} will be ignored.}
}
\value{
A list of lists containing output for each sample in \code{spectra}.

Each sublist has the following elements:
\itemize{
\item \code{original.spect}: The original spectrum with total mutation counts
added to its column name. An additional attribute "exposure" from
\code{exposure} is also added.
\item \code{reconstructed.spect}: The reconstructed spectrum using \code{sigs}
and \code{exposure}. Its column name has the total mutation counts and cosine
similarity with the original spectrum.
\item \code{contributing.sigs}: The contributing signatures to the original
spectrum. The column names of each contributing signature has mutation
counts attributed to this signature, its contribution proportion and proposed
etiology(if the etiology is unknown, then will be blank.)
\item \code{distances}: Various distances and similarities between the original
spectrum and \code{reconstructed.spect}.
}
}
\description{
Add contributing signature activity information for multiple spectra
}
\details{
This function calls \code{\link{ReconstructSpectrum}},
\code{\link{LLHSpectrumNegBinom}} and \code{\link{LLHSpectrumMAP}}.
}
\section{Remark}{
 The column names of \code{spectra} should be the same as the
column name of \code{exposure}.
}

\examples{
\dontrun{
spectra <- PCAWG7::spectra$PCAWG$SBS96[, 1:2, drop = FALSE]
exposure <- PCAWG7::exposure$PCAWG$SBS96[, 1:2, drop = FALSE]
sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                 cancer.type = "Biliary-AdenoCA")
retval <- AddSigActivity(spectra, exposure, sigs, sigs.prop)
}
}
\keyword{internal}
