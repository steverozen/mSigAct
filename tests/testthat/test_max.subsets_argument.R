context("Testing max.subsets argument")

test_that("Testing max.subsets argument for MAPAssignActivity", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
  sigs <- PCAWG7::signature$genome$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Lung-AdenoCA")
  sigs <- sigs[, names(sigs.prop), drop = FALSE]
  output.dir <- file.path(tempdir(), "SBS96")
  
  retval <- MAPAssignActivity(
    spectra                 = spectra,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = DefaultManyOpts(likelihood.dist = "neg.binom"),
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 10,
    seed                    = 8787,
    max.subsets             = 1)
  expect_null(retval$proposed.assignment)
  expect_equal(length(retval$error.messages), 2)
  
  # One sample will have NULL return for proposed.assignment
  retval2 <- MAPAssignActivity(
    spectra                 = spectra,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = DefaultManyOpts(likelihood.dist = "neg.binom"),
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 10,
    seed                    = 8787,
    max.subsets             = 10)
  
  expect_true(!is.null(retval2$proposed.assignment))
  expect_equal(length(retval2$error.messages), 1)
  
})

test_that("Testing max.subsets argument for MAPAssignActivity1", {
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spect <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- PCAWG7::signature$genome$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Lung-AdenoCA")
  sigs <- sigs[, names(sigs.prop), drop = FALSE]
  
  retval <- MAPAssignActivity1(
    spect                   = spect,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = DefaultManyOpts(likelihood.dist = "neg.binom"),
    max.mc.cores            = 10,
    seed                    = 8787,
    max.subsets             = 10)
   
  expect_null(retval$proposed.assignment)
})