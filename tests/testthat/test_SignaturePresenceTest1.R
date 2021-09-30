context("SignaturePresenceTest1")

test_that("Use signature id in SignaturePresenceTest1", {
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- PCAWG7::signature$genome$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.test <- sigs[, names(sigs.prop), drop = FALSE]
  test.out <- SignaturePresenceTest1(spectrum = spectra,
                                     sigs = sigs.to.test,
                                     target.sig.index = "SBS22",
                                     seed = 2892)
  
  test.out1 <- SignaturePresenceTest1(spectrum = spectra,
                                     sigs = sigs.to.test,
                                     target.sig.index = 13,
                                     seed = 2892)
  
  test.out2 <- 
    expect_error(SignaturePresenceTest1(spectrum = spectra,
                                        sigs = sigs.to.test,
                                        target.sig.index = "SBS10a",
                                        seed = 2892))
  expect_equal(test.out$statistic, 155.937, tolerance = 1e-5)
  expect_equal(test.out, test.out1)
})