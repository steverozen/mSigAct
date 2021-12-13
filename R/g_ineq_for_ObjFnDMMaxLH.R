#' Function to constrain the sum of estimated exposures to the number of
#' mutations in the spectrum.
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
#' @param cp.factor Concentration parameter factor. When calculating the
#'   Dirichlet multinomial likelihood, the concentration parameters \code{alpha}
#'   is calculated as follows: 
#'   \code{cp.factor} * \code{reconstruction} / \code{sum(reconstruction)}  
#' (Ignored in this function but used by \code{\link[nloptr]{nloptr}}.)
#'
g_ineq_for_ObjFnDMMaxLH <- function(exp, # Parameters to optimize
                                    spectrum,
                                    sigs,
                                    cp.factor
) {
  
  retval <- abs(sum(exp) - sum(spectrum))
  # message("gi ", retval)
  return(retval)
}
