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
#' @param m.opts If \code{NULL} use the return from
#'    calling \code{\link{DefaultManyOpts}}. For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'
#' @param mc.cores Number of cores to use. Always silently
#'  changed to 1 on Microsoft Windows.
#'
#' @export

SignaturePresenceTest <- function(
  spectra, sigs, target.sig.index, m.opts = NULL, mc.cores = 10) {

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
      m.opts           = m.opts)

  return(out.res)
}
