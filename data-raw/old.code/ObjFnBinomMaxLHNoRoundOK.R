#' A deprecated negative binomial maximum likelihood objective function.
#'
#' Use \code{\link{ObjFnBinomMaxLHRound}} instead.
#'
#' This function rounds sometimes, which leads to
#' minor differences in log likelihoods of reconstructed spectra
#' (\code{\link{LLHSpectrumNegBinom}})
#' compared to the value returned by this function.
#'
#' @inheritParams ObjFnBinomMaxLHRound
#'
#' @export
#'
ObjFnBinomMaxLHNoRoundOK <- function(exp, spectrum, sigs, nbinom.size) {
  warning("Called ObjFnBinomMaxLHNoRoundOK")
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = TRUE, round.ok = TRUE)
}