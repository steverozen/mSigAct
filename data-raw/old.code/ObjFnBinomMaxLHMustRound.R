#' A deprecated negative binomial maximum likelihood objective function.
#'
#' Use \code{\link{ObjFnBinomMaxLHRound}} instead.
#'
#' This function will lead to errors in some situations
#' when the rounded reconstructed signature contains 0s for
#' mutations classes for which the target spectrum is > 0.
#'
#' @inheritParams ObjFnBinomMaxLHRound
#'
#' @export
#'
ObjFnBinomMaxLHMustRound <- function(exp, spectrum, sigs, nbinom.size) {
  warning("Called ObjFnBinomMaxLHMustRound")
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE, round.ok = TRUE)
}