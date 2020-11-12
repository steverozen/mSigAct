#' Optimize the reconstruction of a spectrum from a set of signatures.
#'
#' @param spectrum The spectrum to be reconstructed.
#'
#' @param sigs The available signatures.
#'
#' @param m.opts Options that govern the numerical optimization.
#'    For documentation see \code{\link{DefaultManyOpts}}.
#'
#' @param eval_f The objective function for
#'  \code{\link[nloptr]{nloptr}}. We have only tested
#'  \code{\link{ObjFnBinomMaxLHNoRoundOK}}.
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
  if (m.opts$trace > 10) message("loglh from Nlopter1Tumor = ", loglh)

  if (loglh == Inf && m.opts$trace > 0)
    message("Got -Inf in OptimizeExposure")

  # sum(r$solution) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spectrum) / sum(r$solution)) # sum(recon))

  # TODO, there is redundant info in everything.else
  return(list(loglh = -loglh, exposure = exp, obj.fn.value = r$objective, everything.else = r))
}


#' Optimize assignment of a fixed set of signature activities for a \strong{single} tumor.
#'
#' @inheritParams OptimizeExposure
#'
#' @keywords internal
Nloptr1Tumor <- function(spectrum,
                         sigs,
                         m.opts = NULL,
                         eval_f,
                         ... ) {
  if (!"matrix" %in% class(sigs)) {
    if (!"numeric" %in% class(sigs)) {
      stop("Unexpected class argument: ", class(sigs))
    }
    # Otherwise sigs is a numeric vector, not a matrix. We assume the caller
    # intended it as a single column matrix.
    sigs <- matrix(sigs, ncol = 1)
  }

  if (nrow(sigs) != length(spectrum)) {
    # If we get here there is an error. We save spectrum and sigs, which seems
    # useful for debugging in a call to mclapply (multi-core lapply).
    save(spectrum, sigs, file = 'spectrum.sigs.debug.Rdata')
    stop("nrow(sigs) != length(spectrum), look in file spectrum.sigs.debug.Rdata")
  }

  stopifnot(mode(spectrum) == "numeric")
  foo <- spectrum
  storage.mode(spectrum) <- "double"
  stopifnot(foo == spectrum)

  if (is.null(m.opts)) m.opts <- DefaultManyOpts()

  use.old <- TRUE
  if (ncol(sigs) == 1 || use.old) {
    # x0 is uniform distribution of mutations across signatures
    # (Each signature gets the same number of mutations)
    x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))
  } else {
    stop("This branch does not work")
    # We get x0 that minimizes the Frobenius norm
    # x0 <- sum(spectrum) * as.vector(ICAMS.shiny:::findSigExposures(M = spectrum, P = sigs)$exposures)
    message("QP estimate is ", paste(x0, collapse = " "))
    m.opts$global.opts$maxeval <- 1
    m.opts$local.opts$maxeval <- 30000  }

  if (!is.null(m.opts$global.opts)) { # WARNING, ADDITIONAL CODE WOULD NEED TO BE CHANGED DISABLE THIS BRANCH

    global.res <- nloptr::nloptr(
      x0       = x0,
      eval_f   = eval_f,
      lb       = rep(0, ncol(sigs)),
      ub       = rep(sum(spectrum), ncol(sigs)),
      opts     = m.opts$global.opts,
      spectrum = spectrum,
      sigs     = sigs,
      ...)
    if (m.opts$trace > 0) {
      message("global.res$objective = ", global.res$objective)
    }
    if (global.res$iterations == m.opts$global.opts[["maxeval"]]) {
      # warning("reached maxeval on global optimization: ", global.res$iterations)
    }
  }

  if (use.old) {
    my.x0 <- global.res$solution
  } else {
    my.x0 <- x0
  }
  local.res <- nloptr::nloptr(
    x0=my.x0,
    # x0       = global.res$solution,
    eval_f   = eval_f,
    lb       = rep(0, ncol(sigs)),
    ub       = rep(sum(spectrum) + 1e-2, ncol(sigs)),
    opts     = m.opts$local.opts,
    spectrum = spectrum,
    sigs     = sigs,
    ...)
  if (m.opts$trace > 0)
    message("local.res$objective = ", local.res$objective)
  message("local.res$iterations = ", local.res$iterations)

  if (local.res$iterations == m.opts$local.opts[["maxeval"]]) {
    warning("reached maxeval on local optimization: ", local.res$iterations)
  }

  names(local.res$solution) <- colnames(sigs)
  return(list(objective =  local.res$objective,
              solution   = local.res$solution,
              global.res = global.res,
              local.res  = local.res))
}

