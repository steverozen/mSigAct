#' The negative binomial maximum likelihood objective function.
#'
#' For use by \code{\link[nloptr]{nloptr}}.
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
#' @param never.round Ignore \code{no.round.ok} and never
#'        round; for experimenting and may be removed in future versions.
#'
#' @keywords internal
#'
ObjFnBinomMaxLH2 <-
  function(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE,
           show.warning = FALSE, round.ok = FALSE) {

  if (any(is.na(exp))) return(Inf)
   
  exp <- round(exp)  
  reconstruction <-  sigs %*% exp
  # Maybe faster than ReconstructSpectrum(sigs = sigs, exp = exp)

  if (round.ok) {
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
  }

  ## catch errors with NA in the input or in reconstruction.
  if (any(is.na(reconstruction))) {
    save(reconstruction, spectrum, sigs, exp, nbinom.size,
         file = "construction.error.Rdata")
  }
  stopifnot(!any(is.na(reconstruction)))

  # deleted mSigBG:: -- do a diff
  loglh <- LLHSpectrumNegBinom(spectrum = spectrum,
                               expected.counts = reconstruction,
                               nbinom.size = nbinom.size)

  return(-loglh)
  }

#' The preferred negative binomial maximum likelihood objective function.
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
#' @return
#'
#' -1 * log(likelihood(spectrum | reconstruction))
#'
#' \code{\link[nloptr]{nloptr}} minimizes the objective function, so the
#' lower the objective function, the better.
#'
#' @param exp A vector of exposures ("activities").
#' @param spectrum The spectrum to assess.
#' @param sigs The matrix of signatures.
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.
#'
#' @keywords internal
#'

ObjFnBinomMaxLHRound <- function(exp, spectrum, sigs, nbinom.size) {
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = TRUE, round.ok = FALSE)
}

