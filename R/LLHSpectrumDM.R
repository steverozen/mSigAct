#' Likelihood that 1 observed spectrum was generated from a vector of expected
#' mutation counts using Dirichlet-multinomial distribution
#' 
#' @param spectrum An observed spectrum (a numeric vector)
#'
#' @param expected.counts A vector of expected mutation counts, one expected
#'   count for each mutation type. We want to know the likelihood that this
#'   model generated the observed spectrum, assuming each mutational type
#'   generates counts according to a Dirichlet-multinomial distribution with the
#'   given \code{expected.counts} (argument \code{alpha} to
#'   \code{\link[MGLM]{rdirmn}}).
#'   
#' @param cp.factor Concentration parameter factor. When calculating the
#'   dirichlet multinomial likelihood, the concentration parameters \code{alpha}
#'   is calculated as follows: 
#'   \code{cp.factor} * \code{expected.counts} / \code{sum(expected.counts)}
#'
#' @param verbose If \code{TRUE} print messages under some circumstances.
#'
#' @return \code{log(likelihood(spectrum | expected.counts))}, or, in more
#'   detail, the Dirichlet-multinomial likelihood that each element of the
#'   spectrum (i.e., the count for each mutation type e.g. ACT > AAT) was
#'   generated from the expected count for that mutation type using
#'   Dirichlet-multinomial distribution.
#'   
LLHSpectrumDM <-
  function(spectrum, expected.counts, cp.factor = 1000, verbose = FALSE) {
    stopifnot(length(spectrum) == length(expected.counts))
    
    # We need to convert spectrum to a matrix so that we can transpose later
    if (!inherits(spectrum, "matrix")) {
      spectrum <- as.matrix(spectrum)
    }
    
    if (!inherits(expected.counts, "matrix")) {
      expected.counts <- as.matrix(expected.counts)
    }
    
    prob <- expected.counts[, 1] / sum(expected.counts[, 1])
    
    # Calculate the Dirichlet-multinomial likelihood 
    # Need to transpose spectrum first so that it has the correct format to 
    # be used in MGLM: the columns represent number of categories
    # while rows represent number of observations
    loglh0 <- MGLM::ddirmn(Y = t(spectrum), alpha = prob * cp.factor)
    
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