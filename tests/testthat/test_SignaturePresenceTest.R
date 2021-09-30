context("SignaturePresenceTest")

test_that("Use signature id in SignaturePresenceTest", {
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
  sigs <- PCAWG7::signature$genome$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs.to.test <- sigs[, names(sigs.prop), drop = FALSE]
  test.out <- SignaturePresenceTest(spectra = spectra,
                                    sigs = sigs.to.test,
                                    target.sig.index = "SBS22")
  
  test.out1 <- SignaturePresenceTest(spectra = spectra,
                                     sigs = sigs.to.test,
                                     target.sig.index = 13)
  
  test.out2 <- 
    expect_error(SignaturePresenceTest(spectra = spectra,
                                       sigs = sigs.to.test,
                                       target.sig.index = "SBS10a"))
  expect_equal(test.out[[1]]$statistic, 156.1548, tolerance = 1e-5)
  expect_equal(test.out, test.out1)
})