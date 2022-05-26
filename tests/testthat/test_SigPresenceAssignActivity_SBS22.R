context("SigPresenceAssignActivity for SBS22")

test_that("SigPresenceAssignActivity for SBS Liver tumor", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  library(ICAMS)
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[5], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  
  sig.names <- setdiff(names(sigs.prop), cosmicsig::possible_artifacts())
  sig.names <- setdiff(sig.names, MSISigs()$SBS)
  
  
  output.dir <- "./tests/testthat/output/"
  ICAMS::PlotCatalogToPdf(catalog = spectra,
                          file = file.path(output.dir, "test.spectra.pdf"))
  
  
  
  sigs.to.use <- sigs[, sig.names, drop = FALSE]
  sbs.h8.file <-
    "~/packages/attribution_benchmarking/real_data/source_file/SBS_H8_SBS96.csv"
  sbs.catalog <- 
    ICAMS::ReadCatalog(file = sbs.h8.file, catalog.type = "counts.signature")
  sigs.to.use <- cbind(sigs.to.use, sbs.catalog)
  
  ICAMS::PlotCatalogToPdf(catalog = sigs.to.use,
                          file = file.path(output.dir, "test.signatures.pdf"))
  
  
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 5
  my.opts$trace <- 1e15
  retval.sbs22 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir, "with.sbs22"),
                              max.level = ncol(sigs.to.use) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  pcawg.exp <- PCAWG7::exposure$PCAWG$SBS96
  pcawg.exp2 <- pcawg.exp[, colnames(spectra), drop = FALSE]
  pcawg.assign <- GetSampleSigActivity(spectra = spectra,
                                  exposure = pcawg.exp2,
                                  sigs = sigs, 
                                  sample.names = colnames(spectra), 
                                  output.dir = file.path(output.dir, "pcawg.assign"))
  
  sigs.to.use2 <- sigs.to.use[, colnames(sigs.to.use) != "SBS22", drop = FALSE]
  retval.no.sbs22 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use2,
                              m.opts = my.opts,
                              output.dir = file.path(output.dir, "without.sbs22"),
                              max.level = ncol(sigs.to.use2) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use2),
                              num.parallel.samples = 1,
                              mc.cores.per.sample = 60,
                              seed = 2561, 
                              sig.pres.test.nbinom.size = 5)
  
  
  expect_equal(nrow(retval1$proposed.assignment), 8)
  expect_equal(nrow(retval2$proposed.assignment), 10)
  unlink(file.path(tempdir(), "Liver-HCC-1"), recursive = TRUE)
  unlink(file.path(tempdir(), "Liver-HCC-2"), recursive = TRUE)
})