library(ICAMS) # Needed temporarily pending ICAMS update
MakeHepG2BackgroundPart1()
HepG2.background.info <- MakeHepG2BackgroundPart2(maxeval = 50000)
tmp.sigs <- cbind(HepG2.background.info$background.sig, HepG2.background.info$mean.sig)
colnames(tmp.sigs) <- c("Max.Likelihood.Sig", "Mean.Sig")
attr(tmp.sigs, "catalog.type") <- "counts.signature"
PlotCatalogToPdf(tmp.sigs, "data-raw/spectra.for.background.signatures/2estimates.of.HepG2.background.pdf")
