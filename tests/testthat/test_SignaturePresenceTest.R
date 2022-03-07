context("SignaturePresenceTest")

test_that("Use signature id in SignaturePresenceTest", {
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.test <- sigs[, names(sigs.prop), drop = FALSE]
  test.out <- SignaturePresenceTest(spectra = spectra,
                                    sigs = sigs.to.test,
                                    target.sig.index = "SBS22",
                                    seed = 2892,
                                    mc.cores = 2)
  
  test.out1 <- SignaturePresenceTest(spectra = spectra,
                                     sigs = sigs.to.test,
                                     target.sig.index = 13,
                                     seed = 2892,
                                     mc.cores = 2)
  
  test.out2 <- 
    expect_error(SignaturePresenceTest(spectra = spectra,
                                       sigs = sigs.to.test,
                                       target.sig.index = "SBS10a",
                                       seed = 2892,
                                       mc.cores = 1))
  expect_equal(test.out, test.out1)
})

test_that("Test all signatures in one spectrum", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.test <- sigs[, names(sigs.prop), drop = FALSE]
  test.out <- TestAllSigs(spectrum = spectra, sigs = sigs.to.test, 
                          seed = 2892, mc.cores = 30)
  expect_equal(test.out[1], 3.313776e-05, check.attributes = FALSE)
})



