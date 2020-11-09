#' The negative binomial maximum likelihood objective function.
#'
#' For use by \code{\link[nloptr]{nloptr}}
#'
# @param exp The matrix of exposures ("activities").
# @param spectrum The spectrum to assess.
# @param sigs The matrix of signatures.
# @param nbinom.size The dispersion parameter for the negative
#        binomial distribution; smaller is more dispersed.
#        See \code{\link[stats]{NegBinomial}}.
#'
#' @inheritParams ObjFnBinomMaxLHNoRoundOK
#'
#' @param no.round.ok If \code{TRUE}, allow use of unrounded
#'        reconstruction if some mutation types would have 0
#'        counts in the reconstructed spectrum. Deprecated,
#'        in future will always be \code{TRUE}.
#'
#' @param show.warning Deprecated; ignored.
#'
#' @keywords internal
#'
NEW.VERSION.ObjFnBinomMaxLH2 <-
  function(exp, spectrum, sigs, nbinom.size, no.round.ok = TRUE,
           show.warning = FALSE) {

    if (any(is.na(exp))) return(Inf)

    reconstruction.estimate <-  prop.reconstruct(sigs = sigs, exp = exp)

    # The reconstruction.estimate is an abstract model of the probabilities of
    # the counts of each mutation class given the mixture of exposures in exp
    # (which are also probabilities) and the signatures. The probabilities of
    # the counts do not have to be integers.

    # catch errors with NA in the input or in reconstruction.
    if (any(is.na(reconstruction.estimate))) {
      save(reconstruction.estimate, spectrum, sigs, exp, nbinom.size,
           file = "reconstruction.NA.Rdata")
      stop("There arevNAs in reconstruction; see reconstruction.NA.Rdata")
    }

    if (!no.round.ok) {
      # There can be problems if the rounded reconstruction is 0 for
      # any channel (even if the unrounded reconstruction > 0), because then
      # log likelihood will be -inf. The situation is especial likely
      # to occur if mutation counts in the spectrum are low.
      if (any(round(reconstruction.estimate)[spectrum > 0] == 0)) {
        # warning("Possible problem in rounding reconstucted spectrum; ",
        #        "use eval_f = ObjFnBinomMaxLHNoRoundOK")
      }
      # Use the rounded reconstruction estimate
      reconstruction.estimate <- round(reconstruction.estimate)
    }

    loglh <- LLHSpectrumNegBinom(spectrum = spectrum,
                                 expected.counts = reconstruction.estimate,
                                 nbinom.size = nbinom.size)

    return(-loglh)
  }
