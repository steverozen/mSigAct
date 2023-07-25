context("Testing max.subsets argument")

test_that("Testing max.subsets argument for MAPAssignActivity", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Lung-AdenoCA")
  sigs <- sigs[, names(sigs.prop), drop = FALSE]
  output.dir <- file.path(tempdir(), "SBS96")
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 5
  
  retval <- MAPAssignActivity(
    spectra                 = spectra,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 30,
    seed                    = 8787,
    max.subsets             = 1)
  expect_equal(sum(retval$proposed.assignment), 0)
  expect_equal(length(retval$error.messages), 2)
  
  # One sample will have NULL return for proposed.assignment
  retval2 <- MAPAssignActivity(
    spectra                 = spectra,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 10,
    seed                    = 8787,
    max.subsets             = 10)
  
  expect_equal(colSums(retval2$proposed.assignment), c(0, 40159), check.attributes = FALSE)
  expect_equal(nchar(retval2$error.messages), c(231, 0), check.attributes = FALSE)
  
})

test_that("Testing max.subsets argument for MAPAssignActivity1", {
  indices <- grep("Lung-SCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spect <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Lung-SCC")
  sigs <- sigs[, names(sigs.prop), drop = FALSE]
  my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts$nbinom.size <- 5
  
  retval <- MAPAssignActivity1(
    spect                   = spect,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    max.mc.cores            = 10,
    seed                    = 8787,
    max.subsets             = 2)
   
  expect_equal(sum(retval$proposed.assignment), 0)
  expect_equal(nchar(retval$error.messages), 227)
})