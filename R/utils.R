#' @keywords internal
GetListOfCombinedSpectra <- function(sample.ids, SBS.spectra, ID.spectra) {
  retval <- lapply(sample.ids, FUN = function(x) {
    SBS.spect <- SBS.spectra[, x, drop = FALSE]
    ID.spect <- ID.spectra[, x, drop = FALSE]
    list(SBS.spect, ID.spect)
  })
  retval1 <- do.call("c", retval)
  return(retval1)
}