test_that("Compare and plot log likelihood for one SBS kidney sample", {
  spectra <- PCAWG7::spectra$PCAWG$SBS96
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  kidney_prop <-
    ExposureProportions(
      mutation.type = "SBS96",
      cancer.type = "Kidney-RCC"
    )
  kidney_sigs <- sigs[, names(kidney_prop), drop = FALSE]

  kidney_spectrum <- spectra[, "Kidney-RCC::SP102755", drop = FALSE]

  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 1

  retval <-
    SignaturePresenceTest(
      spectra = kidney_spectrum,
      sigs = kidney_sigs,
      target.sig.index = "SBS22",
      m.opts = my.opts,
      seed = 3819,
      mc.cores = 1
    )
  p.value <- retval[[1]][4]

  loglh.file <- file.path(tempdir(), "sbs22.loglh.pdf")

  xx <- CompareAndPlotLoglh(
    spectrum = kidney_spectrum,
    exposure = retval[[1]]$exp.with,
    sigs = kidney_sigs,
    sig.to.test = "SBS22",
    neg.binom.size = 1,
    file = loglh.file
  )
  expect_true(xx)

  unlink(x = loglh.file, recursive = TRUE)
})
