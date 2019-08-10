
#' Test \code{FindSignatureMinusBackground} on background-only spectra.
#' @keywords internal
Test0 <- function(start.b.fraction = 0.9, maxeval = 10000) {
  
  spectra <- ICAMS::ReadCatalog(
    system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct"))
  
  res <- FindSignatureMinusBackground(
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
  res <- FindSignatureMinusBackground(
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

