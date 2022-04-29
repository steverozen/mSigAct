context("SPAA_MAP")

test_that("SPAA_MAP for SBS Liver tumor", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[92], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.use <- sigs[, names(sigs.prop), drop = FALSE]
  retval1 <- 
    SPAA_MAP(spectra = spectra,
             sigs = sigs.to.use,
             sigs.presence.prop = sigs.prop,
             output.dir = file.path(tempdir(), "Liver-HCC-1"),
             max.level = ncol(sigs.to.use) - 1,
             p.thresh = 0.05 / ncol(sigs.to.use),
             num.parallel.samples = 1,
             mc.cores.per.sample = 100,
             seed = 2561)
  
  retval2 <- 
    MAPAssignActivity(spectra = spectra,
                      sigs = sigs.to.use,
                      sigs.presence.prop = sigs.prop,
                      output.dir = file.path(tempdir(), "Liver-HCC-2"),
                      max.level = ncol(sigs.to.use) - 1,
                      p.thresh = 0.05 / ncol(sigs.to.use),
                      num.parallel.samples = 1,
                      mc.cores.per.sample = 100,
                      seed = 2561, 
                      max.subsets = 1e15)
  expect_lt(retval1$time.for.assignment$elapsed, 
            retval2$time.for.assignment$elapsed)
  expect_equal(nrow(retval1$proposed.assignment), 11)
  expect_equal(nrow(retval2$proposed.assignment), 11)
  unlink(file.path(tempdir(), "Liver-HCC-1"), recursive = TRUE)
  unlink(file.path(tempdir(), "Liver-HCC-2"), recursive = TRUE)
})