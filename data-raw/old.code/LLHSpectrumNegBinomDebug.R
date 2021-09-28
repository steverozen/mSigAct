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