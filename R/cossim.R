#' Cosine similarity with useful argument types
#'
#' @param v1 A vector or single-column matrix
#' @param v2 A vector or single-column matrix
#'
#' @export
#' @examples 
#' spectrum <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
#' SBS96.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
#' exposure <- PCAWG7::exposure$PCAWG$SBS96[, 1, drop = FALSE]
#' reconstructed.spectrum <- ReconstructSpectrum(sigs = SBS96.sigs,
#'                                               exp = exposure,
#'                                               use.sig.names = TRUE)
#' cosine <- cossim(spectrum, reconstructed.spectrum)
cossim <- function(v1, v2) {
  df <- rbind(as.vector(v1),
              as.vector(v2))
  return(suppressMessages(philentropy::distance(x = df, 
                                                method = "cosine", 
                                                test.na = FALSE)))
}
