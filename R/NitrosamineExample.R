RenderNitroExample <- function(whichnitro, bg.inflate.factor = 1) {
  rmarkdown::render(
    input = devtools::package_file(
      file.path("data-raw/NitrosamineExample.Rmd")),
    output_file = file.path(
      devtools::package_file("data-raw"),
      paste0("Nitrosamine-", whichnitro, "-x", bg.inflate.factor, ".html")),
    output_options = list(title      = paste("Background subtraction for", 
                                              whichnitro,
                                              "with bg inflation",
                                              bg.inflate.factor)),
    params = list(
      whichnitro = whichnitro,
      bgfactor   = bg.inflate.factor))
}

if (FALSE) {
  for (mynitro in c("NDEA",
                    "NDMA",
                    "NPIP",
                    "NPYR")) {
    RenderNitroExample(mynitro, bg.inflate.factor = 8)
  }
  RenderNitroExample("NDMA", bg.inflate.factor = 2)
}


if (FALSE) {
  # Possibly turn this into a wrapper function that generates multiple plots?
  # 
  # Make sure ICAMS can generate the correct plots?
  
  # Find the signature minus background and generate plots and a report
  # Two parameters, the background info and the spectra
  # Also location to put the output.
  
  
  NDEA <- mSigAct::nitrosamine.examples$catSBS96[ , 1:2]
  ret <- 
    FindSignatureMinusBackground(
      spectra     = NDEA,
      bg.sig.info = HepG2.background.info,
      m.opts      = NULL,
      start.b.fraction = 0.5)
  
  inferred.sig <- ICAMS::as.catalog(ret$inferred.target.sig, 
                                    catalog.type = "counts.signature",
                                    region = "genome",
                                    infer.rownames = TRUE)
  ICAMS::PlotCatalogToPdf(inferred.sig, "NDEA.inferred.sig2.pdf")
  
  inferred.target.spectra <- inferred.sig %*% matrix(ret$exposures.to.target.sig, nrow = 1)
  inferred.target.spectra <- ICAMS::as.catalog(inferred.target.spectra,
                                               catalog.type   = "counts",
                                               region         = "genome",
                                               infer.rownames = TRUE)
  ICAMS::PlotCatalogToPdf(inferred.target.spectra, "NDEA.inferred.target.spectra2.pdf")
  total.counts <- apply(NDEA, MARGIN = 2, sum)
  bg.counts <- total.counts - ret$exposures.to.target.sig
  inferred.bg.spectra <- HepG2.background.info$background.sig %*% matrix(bg.counts, nrow = 1)
  inferred.bg.spectra <- ICAMS::as.catalog(inferred.bg.spectra,
                                  catalog.type   = "counts",
                                  region         = "genome",
                                  infer.rownames = TRUE)
  ICAMS::PlotCatalogToPdf(inferred.bg.spectra, "NDEA.inferred.bg.spectra2.pdf")
  
  total.spectra <- inferred.target.spectra + inferred.bg.spectra
  ICAMS::PlotCatalogToPdf(total.spectra, "NDEA.reconstructed2.pdf")
  mSigAct:::cossim(total.spectra[ ,1], NDEA[ ,1])
  mSigAct:::cossim(total.spectra[ ,2], NDEA[ ,2])
  dist(rbind(total.spectra[ ,1], NDEA[ ,1]), method = "euclidean")
  dist(rbind(total.spectra[ ,2], NDEA[ ,2]), method = "euclidean")
}
