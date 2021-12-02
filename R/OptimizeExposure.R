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
#' @return
#' A list with elements \describe{
#' \item{\code{loglh}}{The log likelihood of the best
#'       solution (set of exposures) found.}
#' \item{\code{exposure}}{The vector of exposures that generated \code{loglh}, i.e.
#'    the number of mutations ascribed to each signature.}
#' \item{\code{objective}}{The final value of the objective function.}
#' \item{\code{solution}}{The optimum exposures. Deprecated.}
#' \item{\code{warnings}}{A character vector of warnings.}
#' \item{\code{global.search.diagnostics}}{Diagnostics from \code{\link[nloptr]{nloptr}}.}
#' \item{\code{local.search.diagnostics}}{Diagnostics from \code{\link[nloptr]{nloptr}}.}
#' }
#'
#'
#' @keywords internal
#'
OptimizeExposure <- function(spectrum,
                             sigs,
                             m.opts,
                             ...) {

  stopifnot(mode(spectrum) == "numeric")
  
  if (!is.null(m.opts$nbinom.size)) {
    r <- Nloptr1Tumor(spectrum    = spectrum,
                      sigs        = sigs,
                      m.opts      = m.opts,
                      nbinom.size = m.opts$nbinom.size,
                      ...)
    
  } else if (!is.null(m.opts$num.replicates)) {
    r <- Nloptr1Tumor(spectrum    = spectrum,
                      sigs        = sigs,
                      m.opts      = m.opts,
                      num.replicates = m.opts$num.replicates,
                      ...)
  } else {
    r <- Nloptr1Tumor(spectrum    = spectrum,
                      sigs        = sigs,
                      m.opts      = m.opts,
                      ...)
  }

  loglh <- r$objective
  if (m.opts$trace > 10) message("Negative loglh from Nlopter1Tumor = ", loglh)

  if (loglh == Inf && m.opts$trace > 0)
    message("Got -Inf in OptimizeExposure")

  # sum(r$solution) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spectrum) / sum(r$solution)) # sum(recon))

  # message("TEST TEST sum(spectrum) = ", sum(spectrum), "; sum(r$solution) = ", sum(r$solution))

  r$solution <- NULL
  r$objective <- NULL

  return(c(list(loglh = -loglh,
              exposure = exp),
              r))
}
