#' Optimize the reconstruction of a spectrum from a set of signatures.
#'
#' @param spectrum The spectrum to be reconstructed.
#'
#' @param sigs The available signatures.
#'
#' @param m.opts Options that govern the numerical optimization.
#'    For documentation see \code{\link{DefaultManyOpts}}.
#'
#' @param ... Additional arguments for \code{eval_f}.
#'
#' Returns a list with elements \describe{
#' \item{\code{loglh}}{The log likelihood of the best
#'       solution (set of exposures) found.
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
                             ...) {

  stopifnot(mode(spectrum) == "numeric")

  r <- Nloptr1Tumor(spectrum    = spectrum,
                    sigs        = sigs,
                    m.opts      = m.opts,
                    nbinom.size = m.opts$nbinom.size,
                    ...)
  loglh <- r$objective
  if (m.opts$trace > 10) message("loglh from Nlopter1Tumor = ", loglh)

  if (loglh == Inf && m.opts$trace > 0)
    message("Got -Inf in OptimizeExposure")

  # sum(r$solution) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spectrum) / sum(r$solution)) # sum(recon))

  message("TEST TEST sum(spectrum) = ", sum(spectrum), "; sum(r$solution) = ", sum(r$solution))

  # TODO, there is redundant info in everything.else
  return(list(loglh = -loglh, exposure = exp, obj.fn.value = r$objective, everything.else = r))
}
