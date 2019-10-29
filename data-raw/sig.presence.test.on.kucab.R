# Test methyleugenol



SPIKE.TEST <- BG.MEG.kucab.bg.test$x0.5.noise

# mSigAct tests below

test.spectra <- mSigAct::KucabToICAMSSpectra(SPIKE.TEST$test.spectra)
options(warn = 2)
ma.test <- FindSignatureMinusBackground(
  test.spectra, 
  mSigAct::kucab.background.info,
  algorithm = 'NLOPT_LN_COBYLA',
  # algorithm = "NLOPT_GN_DIRECT", # does not work well
  maxeval = 10000, 
  print_level = 0,
  xtol_rel = 0.001,  # 0.0001,)
  xtol_abs = 0.0001,
  start.b.fraction = 0.67)

lsa::cosine(ma.test$target.sig, BG.MEG.Test$sig.ICAMS[ , 1])
# [1,] 0.980528 # add.noise = FALSE, ratio = 0.5
# [1,] 0.9761126 # add.noise = TRUE, ratio = 0.5
# [1,] 0.945781 # add.noise = TRUE, ratio = 0.25

m.opts <- DefaultManyOpts()
ma.pre <- mSigAct::SignaturePresenceTest1(
     spectrum = test.spectra[ , 1, drop = FALSE],
     sigs    = as.matrix(cbind(mSigAct::kucab.background.info$background.sig, 
                               ma.test$target.sig)),
     target.sig.index = 1,
     m.opts = m.opts,
     eval_f = mSigAct:::ObjFnBinomMaxLH)

