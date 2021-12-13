#' The Dirichlet-multinomial likelihood objective function
#'
#' Can be used as the
#' objective function for \code{\link{MAPAssignActivity}},
#' \code{\link{MAPAssignActivity1}},
#' \code{\link{SparseAssignActivity}},
#' \code{\link{SparseAssignActivity1}},
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
#' @param cp.factor Concentration parameter factor. When calculating the
#'   Dirichlet multinomial likelihood, the concentration parameters \code{alpha}
#'   is calculated as follows: 
#'   \code{cp.factor} * \code{reconstruction} / \code{sum(reconstruction)}  
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
ObjFnDMMaxLH <- function(exp, spectrum, sigs, cp.factor = 1000) {
  if (any(is.na(exp))) return(Inf)
  
  reconstruction <-  sigs %*% exp
  # Maybe faster than ReconstructSpectrum(sigs = sigs, exp = exp)
  
  ## catch errors with NA in the input or in reconstruction.
  if (any(is.na(reconstruction))) {
    save(reconstruction, spectrum, sigs, exp,
         file = "construction.error.Rdata")
  }
  stopifnot(!any(is.na(reconstruction)))
  
  loglh <- LLHSpectrumDM(spectrum = spectrum,
                         expected.counts = reconstruction,
                         cp.factor = cp.factor)
  
  return(-loglh)
}