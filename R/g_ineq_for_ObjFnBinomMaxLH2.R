#' Function to constrain the sum of estimated exposures to the number of mutations in the spectrum.
#'
#' See \code{\link[nloptr]{nloptr}} to understand how this function is used.
#'
#' @keywords internal
#'
#' @param exp A numeric vector of exposures.
#'
#' @param spectrum The observed spectrum we are trying to reconstruct.
#'
#' @param sigs The signatures with which we are trying to reconstruct
#'  the spectrum. (Ignored in this function but used by \code{\link[nloptr]{nloptr}}.)
#'
#' @param nbinom.size Dispersion parameter.
#' (Ignored in this function but used by \code{\link[nloptr]{nloptr}}.)

g_ineq_for_ObjFnBinomMaxLH2 <- function(exp, # Parameters to optimize
                                        spectrum,
                                        sigs,
                                        nbinom.size
                                        ) {

  retval <- abs(sum(exp) - sum(spectrum))
  # message("gi ", retval)
  return(retval)
}

