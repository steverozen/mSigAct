#' @keywords internal
GetSigActivity <- function(spectra, exposure, sigs, sig.id, output.dir, 
                           cancer.type = NULL) {
  sample.names1 <- colnames(spectra)
  sample.names2 <- colnames(exposure)
  sample.names.diff <- setdiff(sample.names1, sample.names2)
  if (length(sample.names.diff) > 0) {
    stop("Some samples in spectra do not have corresponding exposure ",
         paste(sample.names.diff, collapse = " "))
  }
  
  if (!is.null(cancer.type)) {
    indices <- grep(pattern = cancer.type, x = colnames(exposure))
    exposure <- exposure[, indices, drop = FALSE]
    spectra <- spectra[, colnames(exposure), drop = FALSE]
  }
  
  exposure.one.sig <- exposure[sig.id, ]
  non.zero.exposure.samples <- names(exposure.one.sig[exposure.one.sig > 0])
  
  if (length(non.zero.exposure.samples) > 0) {
    spectra.to.use <- spectra[, non.zero.exposure.samples, drop = FALSE]
    exposure.to.use <- exposure[, non.zero.exposure.samples, drop = FALSE]
    
    sig.activity <- AddSigActivity(spectra = spectra.to.use, 
                                   exposure = exposure.to.use, 
                                   sigs = sigs, 
                                   use.sparse.assign = TRUE)
    ShowSigActivity(list.of.sig.activity = sig.activity, 
                    output.dir = output.dir)
    return(RemoveZeroActivitySig(exposure.to.use))
  } else {
    msg <- paste0("There are no non-zero exposure for ", sig.id)
    if (is.null(cancer.type)) {
      warning(msg)
    } else {
      warning(msg, " in cancer type ", cancer.type)
    }
    invisible(NULL)
  }
}
