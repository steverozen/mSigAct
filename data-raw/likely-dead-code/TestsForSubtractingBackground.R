PlotFactorizations <- function(out.dir,
                               spectra,
                               bg.sig.info,
                               solution,
                               sig.number = 96,
                               ref.genome = NULL,
                               region = "genome")
{
  if (!dir.exists(out.dir)) {
    if (!dir.create(out.dir, recursive = TRUE)) {
      stop("Cannot create ", out.dir)
    }
  }
  sig <- Solution2Signature(solution,
                            sig.number,
                            ref.genome,
                            region)
  
  b <- solution[(sig.number + 1):length(solution)]
  if (length(b) != ncol(spectra)) {
    stop("The number of estimates of the contribution of the target sequence (",
         length(b), ")\n",
         "does not match the number of input spectra (", ncol(spectra), ")")
  }
  total.counts <- colSums(spectra)
  for (i in 1:ncol(spectra)) {
    bg.counts <- round(b[i] * mSigAct::HepG2.background.info$background.sig)
    attr(bg.counts, "catalog.type") <- "counts"
    sig.counts <- round((total.counts[i] - b[i]) * sig)
    attr(sig, "catalog.type") <- "counts"
    tmp <- cbind(spectra[ , i, drop = FALSE],
                 bg.counts,
                 sig.counts,
                 spectra[ , i, drop = FALSE] - bg.counts)
    
    # TODO(Steve) average the spectra minus the counts and see
    # what they look like
    # TODO(get the pcawg signatures and add them in at different
    # concentrations, with and without noise)
    name <- colnames(spectra)[i]
    colnames(tmp) <- c("Orig", "BG", "Exp*Sig", "Orig-BG")
    ICAMS::PlotCatalogToPdf(tmp, paste0(out.dir, "/", name, ".pdf"))
  }
  return(data.frame(sample = colnames(spectra),
                    spectrum.count = total.counts,
                    bg.count  = b,
                    target.sig.count = total.counts - b))
}



# Calls to SeparateSignatureFromBackground are out of date
# 
#' Test \code{SeparateSignatureFromBackground} on background-only spectra.
#' @keywords internal
Test0 <- function(start.b.fraction = 0.9, maxeval = 10000) {
  
  spectra <- ICAMS::ReadCatalog(
    system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct"))
  
  res <- SeparateSignatureFromBackground(
    spectra = spectra,
    bg.sig.info = mSigAct::HepG2.background.info,
    maxeval=maxeval, 
    print_level=1,
    start.b.fraction = start.b.fraction)
  
  return(res)
  
}

PlotTest0 <- function(out.dir, retval)  {
  PlotFactorizations(out.dir,
                     spectra = mSigAct::HepG2.background.spectra,
                     bg.sig.info = mSigAct::HepG2.background.info,
                     solution = retval$solution)
}

Test1 <- function(start.b.fraction = 0.1, maxeval = 10000) {
  res <- SeparateSignatureFromBackground(
    spectra = mSigAct::cisplatin.exposed.HepG2.96,
    bg.sig.info = mSigAct::HepG2.background.info,
    maxeval=maxeval, 
    print_level=1,
    start.b.fraction = start.b.fraction)
  
  return(res)
}

PlotTest1 <- function(out.dir, test1.retval)  {
  PlotFactorizations(out.dir,
                     spectra = mSigAct::cisplatin.exposed.HepG2.96,
                     bg.sig.info = mSigAct::HepG2.background.info,
                     solution = test1.retval$solution)
}

