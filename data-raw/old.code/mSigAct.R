DefaultGlobalOpts <- function() {
  return(
    list(algorithm    = "NLOPT_GN_DIRECT",
         xtol_rel      = 1e-9,
         # print_level = print_level,
         print_level   = 0,
         maxeval       = 10000))
}

DefaultLocalOpts <- function() {
  return(list(algorithm   = "NLOPT_LN_COBYLA",
              xtol_rel    = 1e-15,
              print_level = 0,
              maxeval     = 10000))
}

#' Set default options for many functions, especially \code{\link[nloptr]{nloptr}}.
#'
#' @export
#'
#' @return A list with the following elements
#' \describe{
#'   \item{global.opts}{Options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the global optimization phase.}
#'
#'   \item{local.opts}{Options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the local optimization phase.}
#'
#'   \item{nbinom.size}{The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.}
#'
#'   \item{trace}{If > 0 print progress messages.}
#' }
DefaultManyOpts <- function() {
  return(list(
    global.opts = DefaultGlobalOpts(),
    local.opts  = DefaultLocalOpts(),
    nbinom.size = 5,
    trace       = 0
  ))
}

#' Optimize assignment of a fixed set of signature activities for a \strong{single} tumor.
#'
#' @inheritParams OptimizeExposure
#'
#' @keywords internal
old.Nloptr1Tumor <- function(spectrum,
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
    # If we get here there is an error. We save specturm and sigs, which seems
    # useful for debugging in a call to mclapply (multi-core lapply).
    save(spectrum, sigs, file = 'spectrum.sigs.debug.Rdata')
    stop("nrow(sigs) != length(spectrum), look in file spectrum.sigs.debug.Rdata")
  }

  stopifnot(mode(spectrum) == "numeric")
  foo <- spectrum
  storage.mode(spectrum) <- "double"
  stopifnot(foo == spectrum)

  if (is.null(m.opts)) m.opts <- DefaultManyOpts()

  # x0 is uniform distribution of mutations across signatures
  # (Each signature gets the same number of mutations)
  x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))

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
      message("globa.res$objective = ", global.res$objective)
    }
    if (global.res$iterations == m.opts$global.opts[["maxeval"]]) {
      # warning("reached maxeval on global optimization: ", global.res$iterations)
    }
  }

  local.res <- nloptr::nloptr(
    # x0=x0,
    x0       = global.res$solution,
    eval_f   = eval_f,
    lb       = rep(0, ncol(sigs)),
    ub       = rep(sum(spectrum) + 1e-2, ncol(sigs)),
    opts     = m.opts$local.opts,
    spectrum = spectrum,
    sigs     = sigs,
    ...)
  if (m.opts$trace > 0)
    message("local.res$objective = ", local.res$objective)

  if (local.res$iterations == m.opts$local.opts[["maxeval"]]) {
    warning("reached maxeval on local optimization: ", local.res$iterations)
  }

  names(local.res$solution) <- colnames(sigs)
  return(list(objective =  local.res$objective,
              solution   = local.res$solution,
              global.res = global.res,
              local.res  = local.res))
}


# Helper function, given signatures (sigs) and exposures (exp), return a
# *proportional* reconstruction; in general, it is *not necessarily* scaled to
# the actual spectrum counts.
#' error checked (?) function to get reconstructed something?
#'
#' @keywords internal
prop.reconstruct <- function(sigs, exp) {
  stopifnot(length(exp) == ncol(sigs))
  return(as.matrix(sigs) %*% exp)
}

#' A deprecated negative binomial maximum likelihood objective function.
#'
#' Use \code{\link{ObjFnBinomMaxLHNoRoundOK}} instead.
#'
#' This function will lead to errors in some situations
#' when the rounded reconstructed signature contains 0s for
#' mutations classes for which the target spectrum is > 0.
#'
#' @inheritParams ObjFnBinomMaxLHNoRoundOK
#'
#' @export
#'
ObjFnBinomMaxLHMustRound <- function(exp, spectrum, sigs, nbinom.size) {
  message("ObjFnBinomMaxLHMustRound is deprecated, use ObjFnBinomMaxLHNoRoundOK")
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE)
}

#' The preferred negative binomial maximum likelihood objective function.
#'
#' Can be used as the
#' objective function for \code{\link{SparseAssignActivity}},
#' \code{\link{SparseAssignActivity1}},
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
#' @param exp The matrix of exposures ("activities").
#' @param spectrum The spectrum to assess.
#' @param sigs The matrix of signatures.
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.
#'
#' @export
#'
ObjFnBinomMaxLHNoRoundOK <- function(exp, spectrum, sigs, nbinom.size) {
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = TRUE)
}

#' Euclidean reconstruction error.
#'
#' @keywords internal
EDist2SpectRounded <- function(exp, sig.names, spect) {
  reconstruction <-
    prop.reconstruct(
      sigs = PCAWG7::signature$genome$SBS96[ , sig.names],
      exp = round(exp))
  reconstruction <- round(reconstruction)
  class(spect) <- "matrix"
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


#' Euclidean reconstruction error.
#'
#' @keywords internal
EDist2Spect <- function(exp, sig.names, spect) {
  reconstruction <-
    prop.reconstruct(
      sigs =
        PCAWG7::signature$genome$SBS96[ , sig.names], exp = exp)
  class(spect) <- "matrix"
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


#' "Polish" a solution by minimizing Euclidean distance to the input spectrum.
#'
#' This is experimental for testing.
#'
#' @keywords internal
Polish <- function(exp, sig.names, spect) {
  class(spect) <- "matrix" # Otherwise cbind will use the ICAMS catalog
                           # method, which may complain that the reconstructed
                           # output is not a catalog.

  retval <- nloptr::nloptr(x0 = exp,
                            eval_f = EDist2Spect,
                            lb = rep(0, length(exp)),
                            ub = rep(sum(exp), length(exp)),
                            opts = list(algorithm   = "NLOPT_LN_COBYLA",
                                        maxeval     = 1000,
                                        print_level = 0,
                                        xtol_rel    = 0.001,
                                        xtol_abs    = 0.0001),
                            sig.names = sig.names,
                            spect     = spect)

  names(retval$solution) <- names(exp)
  return(retval$solution)

}


Adj.mc.cores <- function(mc.cores) {
  if (Sys.info()["sysname"] == "Windows" && mc.cores > 1) {
    message("On Windows, changing mc.cores from ", mc.cores, " to 1")
    return(1)
  }
  return(mc.cores)
}

#' Test whether a given signature is plausibly present in a catalog
#'
#' @param spectra The catalog (matrix) to analyze. This could be
#'   an \code{\link[ICAMS]{ICAMS}} catalog or a numerical matrix.
#'
#' @param sigs A catalog of signatures from which to choose.
#' This could be
#'   and \code{\link[ICAMS]{ICAMS}} catalog or a numerical matrix.
#'
#' @param target.sig.index The index of the signature the presence
#' of which we want to test.
#'
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#'
#' @param m.opts If \code{NULL} use the return from
#'    calling \code{\link{DefaultManyOpts}}. For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'
#' @param mc.cores Number of cores to use. Always silently
#'  changed to 1 on Microsoft Windows.
#'
#' @export

SignaturePresenceTest <- function(
  spectra, sigs, target.sig.index, m.opts = NULL, eval_f, mc.cores = 10) {

  # check if signatures sum to 1
  all.col.sums <- colSums(sigs)

  if (!all.equal(all.col.sums,
                 rep(1, ncol(sigs)),
                 tolerance = 1e-3,
                 check.names = FALSE)) {
    stop(
      "Argument sigs does not seem to be a signature ",
      "catalog because some columns do not sum to 1")
  }

  # Need to match exactly one signature name
  stopifnot(length(target.sig.index) == 1)

  spectra.as.list <- split(t(spectra), 1:ncol(spectra))

  names(spectra.as.list) <- colnames(spectra)

  if (is.null(m.opts)) m.opts <- DefaultManyOpts()

  mc.cores <- Adj.mc.cores(mc.cores)

  out.res <-
    parallel::mclapply(
      X                = spectra.as.list,
      FUN              = SignaturePresenceTest1,
      mc.cores         = mc.cores,
      sigs             = sigs,
      target.sig.index = target.sig.index,
      m.opts           = m.opts,
      eval_f           = eval_f)

  return(out.res)
}
