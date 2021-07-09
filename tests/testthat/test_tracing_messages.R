context("Test tracing messages")

test_that("Test tracing messages", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- PCAWG7::signature$genome$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  my.opts <- DefaultManyOpts()
  MAP.out <- MAPAssignActivity(spectra = spectra,
                               sigs = sigs,
                               sigs.presence.prop = sigs.prop,
                               output.dir = file.path(tempdir(), "Liver-HCC.1"),
                               max.level = length(sigs.prop) - 1,
                               p.thresh = 0.05 / ncol(spectra),
                               m.opts = my.opts,
                               num.parallel.samples = 1,
                               mc.cores.per.sample = 1)
  
  my.opts$trace <- 1
  MAP.out <- MAPAssignActivity(spectra = spectra,
                               sigs = sigs,
                               sigs.presence.prop = sigs.prop,
                               output.dir = file.path(tempdir(), "Liver-HCC.2"),
                               max.level = length(sigs.prop) - 1,
                               p.thresh = 0.05 / ncol(spectra),
                               m.opts = my.opts,
                               num.parallel.samples = 1,
                               mc.cores.per.sample = 1)
  
  unlink(file.path(tempdir(), "Liver-HCC.1"), recursive = TRUE)
  unlink(file.path(tempdir(), "Liver-HCC.2"), recursive = TRUE)
})