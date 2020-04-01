#' The negative binomial maximum likelihood objective function.
#' 
#' For use by \code{\link[nloptr]{nloptr}}
#'
#' @param exp The matrix of exposures ("activities").
#' @param spectrum The spectrum to assess.
#' @param sigs The matrix of signatures.
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.
#'
#' @param no.round.ok If \code{TRUE}, allow use of unrounded
#'        reconstruction if some mutation types would have 0
#'        counts in the reconstructed spectrum.
#' @param show.warning If \code{TRUE} print warning if unrounded
#'        reconstructions were used.
#'
#' @keywords internal
#' 
ObjFnBinomMaxLH2 <- 
  function(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE,
           show.warning = FALSE) {
  
  if (any(is.na(exp))) return(Inf)
  
  reconstruction <-  prop.reconstruct(sigs = sigs, exp = exp)
  
  reconstruction2 <- round(reconstruction)
  # Will cause problems if round of the reconstruction is 0 for
  # any channel even if the reconstruction > 0, because then
  # log likelihood will be -inf. The situation is especial likely
  # to occur if mutation counts in the spectrum are low.
  if (any(reconstruction2 == 0) && no.round.ok) {
    if (show.warning) warning("unrounded reconstruction")
  } else {
    reconstruction <- reconstruction2
  }
  rm(reconstruction2)
  
  ## catch errors with NA in the input or in reconstruction.
  if (any(is.na(reconstruction))) {
    save(reconstruction, spectrum, sigs, exp, nbinom.size,
         file = "reconstruction.error.R")
  }
  stopifnot(!any(is.na(reconstruction)))
  
  loglh <- LLHSpectrumNegBinom(spectrum = spectrum, 
                                expected.counts = reconstruction,
                                nbinom.size = nbinom.size)

  return(-loglh)
}
