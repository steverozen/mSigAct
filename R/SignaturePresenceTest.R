#' Test whether a given signature is plausibly present in a catalog.
#'
#' @param spectra The catalog (matrix) to analyze. This could be
#'   an \code{\link[ICAMS]{ICAMS}} catalog or a numerical matrix.
#'
#' @param sigs A catalog of signatures from which to choose.
#' This could be
#'   and \code{\link[ICAMS]{ICAMS}} catalog or a numerical matrix.
#'
#' @param target.sig.index The index of the signature the presence
#' of which we want to test. It can also be the signature id (e.g. "SBS22").
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'    
#' @param seed Random seed; set this to get reproducible results. (The
#'   numerical optimization is in two phases; the first, global phase
#'   might rarely find different optima depending on the random
#'   seed.)
#'
#' @param mc.cores Number of cores to use. Always silently
#'  changed to 1 on Microsoft Windows.
#'
#' @export
#' 
#' @return A list of test results for each sample in \code{spectra}.
#' Each sublist contains the following elements:
#' 
#' * loglh.with: The maximum log likelihood of the reconstructed spectrum using
#' all the signatures.
#'
#' * loglh.without: The maximum log likelihood of the reconstructed spectrum
#' without the target signature.
#'
#' * statistic: Likelihood ratio test statistic.
#'
#' * chisq.p: P-value of the likelihood ratio test. The null hypothesis is we
#' can plausibly reconstruct the spectrum without the target signature.
#'
#' * exp.with: The exposure using all the signatures which generates the maximum
#' log likelihood \code{loglh.with}.
#'
#' * exp.without: The exposure not using the target signature which generates
#' the maximum log likelihood \code{loglh.without}.
#' 
#' @md
#' 
#' @examples 
#' indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
#' spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
#' sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96",
#'                                  cancer.type = "Lung-AdenoCA")
#' sigs.to.use <- sigs[, names(sigs.prop), drop = FALSE]
#' # Test whether SBS17a is plausibly present in the spectra
#' sig.presence.test.out <- SignaturePresenceTest(spectra = spectra,
#'                                                sigs = sigs.to.use, 
#'                                                target.sig.index = "SBS17a",
#'                                                seed = 2581,
#'                                                mc.cores = 2)

SignaturePresenceTest <- function(
  spectra, sigs, target.sig.index, m.opts = DefaultManyOpts(), seed = NULL, mc.cores = 2) {

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

  mc.cores <- Adj.mc.cores(mc.cores)

  out.res <-
    parallel::mclapply(
      X                = spectra.as.list,
      FUN              = SignaturePresenceTest1,
      mc.cores         = mc.cores,
      sigs             = sigs,
      target.sig.index = target.sig.index,
      m.opts           = m.opts,
      seed             = seed)

  return(out.res)
}
