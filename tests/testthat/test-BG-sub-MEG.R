context("Test Background Subtraction on Methyeugenol data from Kucab 2019")

test_that("Meg, 0.25, Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.Test$x0.25.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.control.bg,
    # algorithm = 'NLOPT_LN_COBYLA',
    # maxeval = 10000, 
    # print_level = 0,
    # xtol_rel = 0.001,  # 0.0001,)
    # xtol_abs = 0.0001,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])),
    0.945781, tolerance = 1e-5)
})

test_that("Meg, 0.5, Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.Test$x0.5.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.control.bg,
    # algorithm = 'NLOPT_LN_COBYLA',
    # maxeval = 10000, 
    # print_level = 0,
    # xtol_rel = 0.001,  # 0.0001,)
    # xtol_abs = 0.0001,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])),
    0.9761126, tolerance = 1e-5)
})

test_that("Meg, 0.5, No Noise", {
  test.spectra <- mSigAct::KucabToICAMSSpectra(
    mSigAct::BG.MEG.Test$x0.5.no.noise$test.spectra)
  ma.test <- FindSignatureMinusBackground(
    test.spectra, 
    mSigAct::kucab.control.bg,
    # algorithm = 'NLOPT_LN_COBYLA',
    # maxeval = 10000, 
    # print_level = 0,
    # xtol_rel = 0.001,  # 0.0001,)
    # xtol_abs = 0.0001,
    m.opts = FindSigMinusBGOpt(),
    start.b.fraction = 0.67)
  testthat::expect_equal(
    as.vector(lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])),
    0.980528, tolerance = 1e-5)
})


# lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])
# [1,] 0.980528 # add.noise = FALSE, ratio = 0.5
# [1,] 0.9761126 # add.noise = TRUE, ratio = 0.5
# [1,] 0.945781 # add.noise = TRUE, ratio = 0.25



