context("Testing CustomizeManyOpts function")

test_that("Testing CustomizeManyOpts function for SBS", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
  sigs <- sigs[, names(sigs.prop), drop = FALSE]
  output.dir <- file.path(tempdir(), "SBS96")
  
  TestFunction <- function(my.opts) {
    retval <- PresenceAssignActivity(
      spectra                 = spectra,
      sigs                    = sigs,
      output.dir              = output.dir,
      p.thresh                = 0.001 / (4 * ncol(sigs)),
      m.opts                  = my.opts,
      num.parallel.samples    = 1,
      mc.cores.per.sample     = 10,
      seed                    = 8257,
      save.files              = FALSE,
      use.forward.search      = TRUE)
    
    return(retval)
  }
  
  my.opts1 <- DefaultManyOpts()
  my.opts1$trace <- 1e10
  retval1 <- TestFunction(my.opts1)
  
  my.loglh.fn1 <- function(spectrum, expected.counts) {
    loglh0 <- 
      stats::dmultinom(x = spectrum, prob = expected.counts, log = TRUE)
    return(loglh0)
  }
  my.opts2 <- CustomizeManyOpts(loglh.fn = my.loglh.fn1)
  my.opts2$trace <- 1e10
  retval2 <- TestFunction(my.opts2)
  
  
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
    m.opts                  = DefaultManyOpts(likelihood.dist = "neg.binom"),
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 10,
    seed                    = 8787,
    max.subsets             = 10)
  
  expect_equal(colSums(retval2$proposed.assignment), c(0, 40159), check.attributes = FALSE)
  expect_equal(nchar(retval2$error.messages), c(231, 0), check.attributes = FALSE)
  
})

test_that("Testing max.subsets argument for MAPAssignActivity1", {
  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spect <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Liver-HCC")
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
  
  expect_equal(sum(retval$proposed.assignment), 0)
  expect_equal(nchar(retval$error.messages), 231)
})