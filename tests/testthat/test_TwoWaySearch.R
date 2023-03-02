test_that("Use forward search for MAPAssignActivity", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.use <- sigs[, names(sigs.prop), drop = FALSE]
  
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 8
  my.opts$trace <- 100
  
  retval1 <- 
    PresenceAssignActivity(spectra = spectra,
                           sigs = sigs.to.use,
                           output.dir = file.path(tempdir(), "Liver-HCC-1"),
                           m.opts = my.opts,
                           max.level = ncol(sigs.to.use) - 1,
                           p.thresh = 0.001 / (4 * ncol(sigs.to.use)),
                           num.parallel.samples = 1,
                           mc.cores.per.sample = 10,
                           seed = 2561, 
                           use.two.way.search = TRUE)
  
})
  