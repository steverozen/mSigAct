#' Optimize the reconstruction of a spectrum from a set of signatures.
#'
#' @param spectrum The spectrum to be reconsructed.
#'
#' @param sigs The available signatures.
#'
#' @param m.opts Options that govern the numerical optimizaiton.
#'    For documentation see \code{\link{DefaultManyOpts}}.
#'
#' @param eval_f The objective function for
#'  \code{\link[nloptr]{nloptr}}. We have only tested
#'  \code{\link{ObjFnBinomMaxLHNoRoundOK}} and
#'  \code{\link{ObjFnBinomMaxLHMustRound}}.
#'
#' @param ... Additional arguments for \code{eval_f}.
#'
#' Returns a list with elements \describe{
#' \item{\code{loglh}}{The log likelihood of the best solution (set of exposures) found.
#'       For a more general objective function this might be \code{NULL}.}
#' \item{\code{exposure}}{The vector of exposures that generate \code{loglh}, i.e.
#'    the number of mutations ascribed to each signature.}
#' \item{\code{obj.fn.value}}{The objective function value associated with \code{exposure}.}
#' \item{\code{everything.else}}{Everything returned by the function \code{\link{Nloptr1Tumor}.} }
#' }
#'
#' @export
#'
OptimizeExposure <- function(spectrum,
                             sigs,
                             m.opts,
                             eval_f,
                             ...) {

  stopifnot(mode(spectrum) == "numeric")

  r <- Nloptr1Tumor(spectrum = spectrum,
                    sigs,
                    m.opts = m.opts,
                    eval_f      = eval_f,
                    nbinom.size = m.opts$nbinom.size,
                    ...)
  loglh <- r$objective

  if (loglh == Inf && m.opts$trace > 0)
    message("Got -Inf in OptimizeExposure")

  # sum(r$solution) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spectrum) / sum(r$solution)) # sum(recon))

  # TODO, there is redundant info in everything.else
  return(list(loglh = -loglh, exposure = exp, obj.fn.value = r$objective, everything.else = r))
}

