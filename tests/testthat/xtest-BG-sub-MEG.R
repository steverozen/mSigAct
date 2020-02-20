context("Test Background Subtraction on Methyeugenol data from Kucab 2019")

test_that("Meg, 0.25, Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.kucab.bg.test$x0.25.noise$test.spectra)
  set.seed(1099,"L'Ecuyer-CMRG")
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.background.info,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(
      # lsa::cosine
      cossim(ma.test$target.sig,
                  mSigAct::BG.MEG.kucab.bg.test$sig.ICAMS[ , 1])),
    0.9457784, 
    tolerance = 2e-2)
  testthat::expect_equal(
    ma.test$exposures.to.target.sig, 
    c(215.5293, 241.6403, 195.5088),
    tolerance = 0.5)
})

test_that("Meg, 0.5, Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.kucab.bg.test$x0.5.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.background.info,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(
      # lsa::cosine
      cossim(ma.test$target.sig,
                  mSigAct::BG.MEG.kucab.bg.test$sig.ICAMS[ , 1])),
    # 0.9761126, 
    0.9760834,
    tolerance = 1e-2)
})

test_that("Meg, 0.5, No Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.kucab.bg.test$x0.5.no.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.background.info,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(
      # lsa::cosine
      cossim(ma.test$target.sig,
                  mSigAct::BG.MEG.kucab.bg.test$sig.ICAMS[ , 1])),
    0.9805251,
    tolerance = 1e-2)
})
