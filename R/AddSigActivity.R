#' Add contributing signature activity information for one sample
#'
#' @details This function calls \code{\link{ReconstructSpectrum}} and
#'   \code{\link{LLHSpectrumNegBinom}}.
#'
#' @param spect A single spectrum.
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param sigs The signatures with which we are trying to reconstruct \code{spect}.
#' A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog. The column
#' names of \code{sigs} should be a superset of row names of \code{exposure}.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being a subset of \code{colnames(sigs)}. See \code{\link{ExposureProportions}}
#'    for more details.
#'
#' @param nbinom.size The dispersion parameter for the negative binomial
#'   distribution; smaller is more dispersed. See
#'   \code{\link[stats]{NegBinomial}}.
#'
#' @return A list of elements:
#' * \code{original.spect}: The original \code{spect} with total mutation counts
#' added to its column name. An additional attribute "exposure" from
#' \code{exposure} is also added.
#'
#' * \code{reconstructed.spect}: The reconstructed spectrum using \code{sigs}
#' and \code{exposure}. Its column name has the total mutation counts and cosine
#' similarity with the original \code{spect}.
#'
#' * \code{contributing.sigs}: The contributing signatures to the original
#' \code{spect}. The column names of each contributing signature has mutation
#' counts attributed to this signature, its contribution proportion and proposed
#' etiology. (If the etiology is unknown, then will be blank.)
#'
#' * \code{distances}: Various distances and similarities between the original
#' spectrum and \code{reconstructed.spect}.
#'
#' @section Remark: The column names of \code{spect} should be the same as the
#'   column name of \code{exposure}.
#'
#' @md
#'
#' @keywords internal
AddSigActivity1 <- function(spect, exposure, sigs,
                            sigs.presence.prop, nbinom.size = 5) {
  if (ncol(spect) != ncol(exposure)) {
    stop("The number of samples in spectrum is not equal to the number of ",
         "samples in exposure")
  }

  if (colnames(spect) != colnames(exposure)) {
    stop("The sample name in spectrum ", colnames(spect), " is not the same ",
         "as the sample name in exposure ", colnames(exposure))
  }

  attr(spect, "exposure") <- exposure

  exposure <- as.matrix(exposure) # In case it is a data frame

  # Remove rows with zero exposure
  exposure <- exposure[exposure > 0, , drop = FALSE]

  sigs.not.available <- setdiff(rownames(exposure), colnames(sigs))
  if (length(sigs.not.available) > 0) {
    stop("Some signatures used in exposure is not available in the signatures ",
         "provided ", paste(sigs.not.available, collapse = " "))
  }

  # Sort exposure according to mutation counts
  exposure <- exposure[order(exposure[, 1], decreasing = TRUE), , drop = FALSE]

  sigs.names <- rownames(exposure)
  sigs1 <- sigs[, sigs.names, drop = FALSE]

  # Get the mutation type of the spectrum
  mut.type <- GetMutationType(spect)

  if (!is.null(mut.type)) {
    ets <- PCAWG7::GetEtiology(mut.type, colnames(sigs1))
    colnames(sigs1) <-
      paste0(colnames(sigs1), " (exposure = ", round(exposure[, 1]),
             ", contribution = ",
             round(exposure[, 1]/sum(exposure[, 1]), 2), ") ",
             ets)
    # etiologies <- sigs.etiologies[[mut.type]]
    # colnames(sigs1) <-
    #  paste0(colnames(sigs1), " (exposure = ", round(exposure[, 1]),
    #         ", contribution = ",
    #         round(exposure[, 1]/sum(exposure[, 1]), 2), ") ",
    #         etiologies[colnames(sigs1), ])
  } else {
    colnames(sigs1) <-
      paste0(colnames(sigs1), " (exposure = ", round(exposure[, 1]),
             ", contribution = ",
             round(exposure[, 1]/sum(exposure[, 1]), 2), ")")
  }

  recon.spect <- ReconstructSpectrum(sigs = sigs1, exp = exposure)
  distances <-
    DistanceMeasures(spect = spect, recon = recon.spect, nbinom.size = nbinom.size,
                     model = sigs.names, sigs.presence.prop = sigs.presence.prop)
  reconstructed.spectrum <- round(recon.spect)

  colnames(reconstructed.spectrum) <-
    paste0("reconstructed (count = ", round(colSums(reconstructed.spectrum)),
           ", cosine similarity = ",
           round(distances$value["cosine"], 5), ")")
  colnames(spect) <- paste0(colnames(spect), " (count = ",colSums(spect), ")")
  sig.activity <- list(original.spect = spect,
                       reconstructed.spect = reconstructed.spectrum,
                       contributing.sigs = sigs1,
                       distances = distances)
  return(sig.activity)
}

#' Add contributing signature activity information for multiple spectra
#'
#' @details This function calls \code{\link{ReconstructSpectrum}} and
#'   \code{\link{LLHSpectrumNegBinom}}.
#'
#' @param spectra The spectra (multiple spectra) to be reconstructed.
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param sigs The signatures with which we are trying to reconstruct \code{spectra}.
#' A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog. The column
#' names of \code{sigs} should be a superset of row names of \code{exposure}.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being a subset of \code{colnames(sigs)}. See \code{\link{ExposureProportions}}
#'    for more details.
#'
#' @param nbinom.size The dispersion parameter for the negative binomial
#'   distribution; smaller is more dispersed. See
#'   \code{\link[stats]{NegBinomial}}.
#'
#' @return A list of lists containing output for each sample in \code{spectra}.
#'
#' Each sublist has the following elements:
#' * \code{original.spect}: The original spectrum with total mutation counts
#' added to its column name. An additional attribute "exposure" from
#' \code{exposure} is also added.
#'
#' * \code{reconstructed.spect}: The reconstructed spectrum using \code{sigs}
#' and \code{exposure}. Its column name has the total mutation counts and cosine
#' similarity with the original spectrum.
#'
#' * \code{contributing.sigs}: The contributing signatures to the original
#' spectrum. The column names of each contributing signature has mutation
#' counts attributed to this signature, its contribution proportion and proposed
#' etiology(if the etiology is unknown, then will be blank.)
#'
#' * \code{distances}: Various distances and similarities between the original
#' spectrum and \code{reconstructed.spect}.
#'
#' @section Remark: The column names of \code{spectra} should be the same as the
#'   column name of \code{exposure}.
#'
#' @md
#'
#' @export
#'
#' @examples
#' spectra <- PCAWG7::spectra$PCAWG$SBS96[, 1:2, drop = FALSE]
#' exposure <- PCAWG7::exposure$PCAWG$SBS96[, 1:2, drop = FALSE]
#' sigs <- PCAWG7::signature$genome$SBS96
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96",
#'                                  cancer.type = "Biliary-AdenoCA")
#' retval <- AddSigActivity(spectra, exposure, sigs, sigs.prop)
AddSigActivity <-
  function(spectra, exposure, sigs, sigs.presence.prop, nbinom.size = 5) {
  if (ncol(spectra) != ncol(exposure)) {
    stop("The number of samples in spectrum is not equal to the number of ",
         "samples in exposure")
  }

  if (!all(colnames(spectra) == colnames(exposure))) {
    stop("The sample name in spectrum ", paste(colnames(spectra), collapse = " "),
         " is not the same as the sample name in exposure ",
         paste(colnames(exposure), collapse = ""))
  }

  ret <- lapply(1:ncol(spectra), FUN = function(x) {
    spect <- spectra[, x, drop = FALSE]
    expo <- exposure[, x, drop = FALSE]
    out <- AddSigActivity1(spect = spect, exposure = expo, sigs = sigs,
                           sigs.presence.prop = sigs.presence.prop,
                           nbinom.size = nbinom.size)
  })

  names(ret) <- colnames(spectra)
  return(ret)
}
