#' Negative binomial log likelihood that an observed spectrum was generated from an observed vector of expected counts.
#'  
#' @param spectrum An observed spectrum (a numeric vector)
#' 
#' @param expected.counts A vector of expected mutation counts, one
#' expected count for each mutation type. We want to know the 
#' likelihood that this model generated the observed 
#' spectrum, assuming each mutational generates counts according to
#' a negative binomial distribution with 
#' the given \code{expected.counts} and dispersion parameter
#' \code{nbinom.size}. See \code{\link[stats]{nbinom}}.
#' 
#' @param nbinom.size The \code{size} parameter that
#' governs dispersion. See \code{\link[stats]{nbinom}}.
#' Smaller values correspond to larger dispersion.
#'
#' @return  \code{log(likelihood(spectrum | expected.counts))}
#'
#' @export

LogLHNegBinom <-function(spectrum, expected.counts, nbinom.size) {
  
  stopifnot(length(spectrum) == length(expected.counts))
  loglh <- 0
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
  return(loglh)
  # -loglh
}
