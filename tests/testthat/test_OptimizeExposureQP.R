test_that("OptimizeExposureQPBootstrap", {
  rr <- OptimizeExposureQP(spectrum = PCAWG7::spectra$PCAWG$SBS96[ , 1],
                            signatures = PCAWG7::COSMIC.v3.0$signature$genome$SBS96[ , 1:5])
  expect_equal(rr,
               c(SBS1 = 1384.08972804323, SBS2 = 1684.34516525904,
                 SBS3 = 6960.32503149651,  SBS4 = 1196.35199280255,
                 SBS5 = 3674.88808239867))

})




