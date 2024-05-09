test_that("OptimizeExposureQPBootstrap", {
  rr <- OptimizeExposureQP(spectrum = PCAWG7::spectra$PCAWG$SBS96[ , 1],
                            signatures = cosmicsig::COSMIC_v3.0$signature$GRCh37$SBS96[ , 1:5])
  expect_equal(rr,
               c(SBS1 = 1383.896, SBS2 = 1684.294,
                 SBS3 = 6958.131,  SBS4 = 1196.792,
                 SBS5 = 3676.888), tolerance = 1e-5)

})




