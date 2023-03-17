#' Likelihood that 1 observed spectrum was generated from a vector of expected
#' mutation counts using multinomial distribution
#' 
#' @param spectrum An observed spectrum (a numeric vector)
#'
#' @param expected.counts A vector of expected mutation counts, one expected
#'   count for each mutation type. We want to know the likelihood that this
#'   model generated the observed spectrum, assuming each mutational type
#'   generates counts according to a multinomial distribution with the given
#'   \code{expected.counts} (argument \code{prob} to
#'   \code{\link[stats]{Multinom}}).
#'
#' @param verbose If \code{TRUE} print messages under some circumstances.
#'
#' @return \code{log(likelihood(spectrum | expected.counts))}, or, in more
#'   detail, the multinomial likelihood that each element of the spectrum (i.e.,
#'   the count for each mutation type e.g. ACT > AAT) was generated from the
#'   expected count for that mutation type using multinomial distribution.
#'   
#' @importFrom stats dmultinom
#'
#' @keywords internal
LLHSpectrumMultinom <-
  function(spectrum, expected.counts, verbose = FALSE) {
    
    stopifnot(length(spectrum) == length(expected.counts))
    
    if (Sys.getenv("MN_LOGLH_THRESH") == "") {
      thresh_value <- -1
    } else {
      thresh_value <- as.numeric(Sys.getenv("MN_LOGLH_THRESH"))
    }
    indices <- which(spectrum >= thresh_value)
    
    if (length(indices) == 0) {
      indices <- seq_len(length(spectrum))
    }
    
    loglh0 <- stats::dmultinom(x = spectrum[indices],
                               prob = expected.counts[indices],
                               log = TRUE)
    if (is.nan(loglh0)) {
      warning("loglh0 is Nan, changing to -Inf")
      loglh0 = -Inf
    }
    if (loglh0 == -Inf && verbose) {
      message("loglh0== -Inf for spectrum = ")
      message(paste0(spectrum, collapse = " "))
      message("expected.counts = ")
      message(paste(expected.counts, collapse = " "))
    }
    
    stopifnot(mode(loglh0) == 'numeric' )
    return(loglh0)
  }