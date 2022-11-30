test_that("SigPresenceAssignActivity for SBS Liver tumor", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.use <- sigs[, names(sigs.prop), drop = FALSE]
  retval1 <- 
    SigPresenceAssignActivity(spectra = spectra,
                              sigs = sigs.to.use,
                              output.dir = file.path(tempdir(), "Liver-HCC-1"),
                              max.level = ncol(sigs.to.use) - 1,
                              p.thresh = 0.05 / ncol(sigs.to.use),
                              num.parallel.samples = 2,
                              mc.cores.per.sample = 30,
                              seed = 2561)
  
  retval2 <- 
    SparseAssignActivity(spectra = spectra,
                         sigs = sigs.to.use,
                         output.dir = file.path(tempdir(), "Liver-HCC-2"),
                         max.level = ncol(sigs.to.use) - 1,
                         p.thresh = 0.05 / ncol(sigs.to.use),
                         num.parallel.samples = 2,
                         mc.cores.per.sample = 30,
                         seed = 2561, 
                         max.subsets = 1e15)
  expect_lt(retval1$time.for.assignment$elapsed, 
            retval2$time.for.assignment$elapsed)
  expect_equal(nrow(retval1$proposed.assignment), 8)
  expect_equal(nrow(retval2$proposed.assignment), 10)
  unlink(file.path(tempdir(), "Liver-HCC-1"), recursive = TRUE)
  unlink(file.path(tempdir(), "Liver-HCC-2"), recursive = TRUE)
})