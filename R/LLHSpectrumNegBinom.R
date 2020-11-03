#' Likelihood that \strong{1} observed spectrum was generated from a vector of expected counts.
#'
#' Returns the sum of the negative binomial likelihoods
#' that each each element of the
#' spectrum (i.e., the count for each mutation type e.g. ACT > AAT)
#' was generated from expected count for that mutation type.
#'
#' @param spectrum An observed spectrum (a numeric vector)
#'
#' @param expected.counts A vector of (integer) expected mutation counts, one
#' expected count for each mutation type. We want to know the
#' likelihood that this model generated the observed
#' spectrum, assuming each mutational types generates counts according to
#' a negative binomial distribution with
#' the given \code{expected.counts} (argument \code{mu}
#' to \code{\link[stats]{dnbinom}}) and dispersion parameter
#' \code{nbinom.size}.
#'
#' @param nbinom.size The \code{size} parameter that
#' governs dispersion. See \code{\link[stats]{dnbinom}}.
#' Smaller values correspond to larger dispersion.
#'
#' @return \code{log(likelihood(spectrum | expected.counts))}
#'
#' @importFrom stats dnbinom
#'
#' @export

LLHSpectrumNegBinom <-function(spectrum, expected.counts, nbinom.size) {

  stopifnot(length(spectrum) == length(expected.counts))
  loglh <- 0
  loglh0 <- sum(stats::dnbinom(
    x = spectrum, mu = expected.counts, size = nbinom.size, log =TRUE))
  for (i in 1:length(spectrum)) { # Iterate over each channel in the
    # spectrum and sum the log
    # likelihoods.

    nbinom <- stats::dnbinom(x=spectrum[i],
                             mu = expected.counts[i],
                             size=nbinom.size,
                             log = TRUE)

    loglh <-loglh + nbinom
  }
  stopifnot(mode(loglh) == 'numeric' )
  if (!isTRUE(all.equal(loglh, loglh0, tolerance = 1e-11))) {
    stop("loglh != loglh0, difference = ", loglh - loglh0)
  }
  return(loglh)
}



# https://cran.r-project.org/web/packages/DirichletReg/DirichletReg.pdf

