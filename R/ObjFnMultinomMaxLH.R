#' The multinomial likelihood objective function
#'
#' Can be used as the
#' objective function for \code{\link{MAPAssignActivity}},
#' \code{\link{MAPAssignActivity1}},
#' \code{\link{SparseAssignActivity}},
# \code{\link{SparseAssignActivity1}},  ##### Function removed
#' \code{\link{SignaturePresenceTest}}
#' and \code{\link{SignaturePresenceTest1}}.
#' (Internally used by by \code{\link[nloptr]{nloptr}}.)
#'
#' @param exp The matrix of exposures ("activities").
#' 
#' @param spectrum The spectrum to assess.
#' 
#' @param sigs The matrix of signatures.
#'
#' @keywords internal
#' 
#' @return
#'
#' -1 * log(likelihood(spectrum | reconstruction))
#'
#' \code{\link[nloptr]{nloptr}} minimizes the objective function, so the
#' lower the objective function, the better.
#'
ObjFnMultinomMaxLH <- function(exp, spectrum, sigs) {
  
  if (any(is.na(exp))) return(Inf)
  
  reconstruction <-  sigs %*% exp
  # Maybe faster than ReconstructSpectrum(sigs = sigs, exp = exp)
  
  ## catch errors with NA in the input or in reconstruction.
  if (any(is.na(reconstruction))) {
    save(reconstruction, spectrum, sigs, exp,
         file = "construction.error.Rdata")
  }
  stopifnot(!any(is.na(reconstruction)))
  
  loglh <- LLHSpectrumMultinom(spectrum = spectrum,
                               expected.counts = reconstruction)
  
  return(-loglh)
}