#' Likelihood that 1 observed spectrum was generated from a vector of expected
#' mutation counts using negative binomial distribution
#'
#' @param spectrum An observed spectrum (a numeric vector)
#'
#' @param expected.counts A vector of (integer) expected mutation counts, one
#' expected count for each mutation type. We want to know the
#' likelihood that this model generated the observed
#' spectrum, assuming each mutational types generates counts according to
#' a negative binomial distribution with
#' the given \code{expected.counts} (argument \code{mu}
#' to \code{\link[stats]{NegBinomial}}) and dispersion parameter
#' \code{nbinom.size}.
#'
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.
#'
#' @param verbose If \code{TRUE} print messages under some circumstances.
#'
#' @return \code{log(likelihood(spectrum | expected.counts))}, or,
#' in more detail,
#' the sum of the negative binomial likelihoods
#' that each element of the
#' spectrum (i.e., the count for each mutation type e.g. ACT > AAT)
#' was generated from the expected count for that mutation type.
#'
#' @importFrom stats dnbinom
#'
#' @export

LLHSpectrumNegBinom <-
  function(spectrum, expected.counts, nbinom.size, verbose = FALSE) {

    stopifnot(length(spectrum) == length(expected.counts))

    loglh0 <- sum(stats::dnbinom(
      x = spectrum, mu = expected.counts, size = nbinom.size, log =TRUE))
    
    if (is.nan(loglh0)) {
      warning("logl9 is Nan, changing to -Inf")
      loglh0 = -Inf
    }
    if (loglh0 == -Inf && verbose) {
      message("logh0== -Inf for spectrum = ")
      message(paste0(spectrum, collapse = " "))
      message("expected.counts = ")
      message(paste(expected.counts, collapse = " "))
    }

    stopifnot(mode(loglh0) == 'numeric' )
    return(loglh0)
  }


#' A verbose version of \code{\link{LLHSpectrumNegBinom}} for testing
#'
#' We use a separate function so as not to slow down the heavily
#' used \code{\link{LLHSpectrumNegBinom}} and to provide more
#' information in the output
#'
#' @inheritParams LLHSpectrumNegBinom
#'
#' @export
#'
#' @return A \code{\link[tibble]{tibble}} with
#'    self-explanatory columns and rows.

LLHSpectrumNegBinomDebug <-
  function(spectrum, expected.counts, nbinom.size, verbose = FALSE) {

    stopifnot(length(spectrum) == length(expected.counts))

    loglik <- stats::dnbinom(
      x = spectrum, mu = expected.counts, size = nbinom.size, log =TRUE)
    rr <- tibble::tibble(mut.type       = names(spectrum),
                         log.likelihood = loglik,
                         spectrum       = spectrum,
                         expected       = expected.counts,
                         spec.minus.exp = spectrum - expected.counts
    )

    return(rr)
  }


# https://cran.r-project.org/web/packages/DirichletReg/DirichletReg.pdf
