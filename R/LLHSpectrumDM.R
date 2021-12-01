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
#' @param num.replicates Number of bootstrap replicates for \code{expected.counts}.   
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
  function(spectrum, expected.counts, num.replicates = 1000, verbose = FALSE) {
    
    stopifnot(length(spectrum) == length(expected.counts))
    
    # First bootstrap from expected.counts to get reconstructed spectrum
    # replicates This is needed for later fitting Dirichlet distribution to get
    # the maximum likelihood estimates of alpha 
    # (If we do not do bootstrapping, we only have one sample and the MLE
    # estimates of alpha in Dirichlet distribution will not be accurate)
    replicates <- stats::rmultinom(n = num.replicates, 
                                   size = sum(expected.counts), 
                                   prob = expected.counts)
    
    # Transpose replicates so that the columns represent mutation types
    replicates.t <- t(replicates)
    
    # Fit Dirichlet distribution to the reconstructed spectrum replicates to
    # get MLE estimates of alpha in Dirichlet distribution
    fit.dirichlet <- compositions::fitDirichlet(replicates.t)
    
    # Get the MLE estimates of alpha for each mutation type
    fit.dirichlet <- tryCatch(expr = {
      compositions::fitDirichlet(replicates.t)
    },
    error = function(errors) {
      # Error can occur when some of the cell values in the replicates.t matrix
      # are zero
      if (grep("system is computationally singular", errors$message)) {
        return(NA)
      }
    })
    
    if (is.na(fit.dirichlet)) {
      # Replace zeros in replicates.t matrix with the smallest positive number
      replicates.t[replicates.t == 0] <- .Machine$double.xmin
      
      # Refit Dirichlet distribution to replicates.t
      fit.dirichlet <- compositions::fitDirichlet(replicates.t)
    }
    
    alpha.mle <- fit.dirichlet$alpha
    
    # Calculate the Dirichlet-multinomial likelihood 
    loglh0 <- MGLM::ddirmn(Y = spectrum, alpha = alpha.mle)
    
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