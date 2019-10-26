context("Test Background Subtraction on Methyeugenol data from Kucab 2019")

test_that("Meg, 0.25, Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.Test$x0.25.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.control.bg,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])),
    0.9457784, tolerance = 1e-5)
})

test_that("Meg, 0.5, Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.Test$x0.5.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.control.bg,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])),
    # 0.9761126, 
    0.9760834,
    tolerance = 1e-5)
})

test_that("Meg, 0.5, No Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.Test$x0.5.no.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.control.bg,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])),
    0.9805251, tolerance = 1e-5)
})



