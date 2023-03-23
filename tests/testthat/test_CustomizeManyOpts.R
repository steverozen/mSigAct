context("Testing CustomizeManyOpts function")

test_that("Testing CustomizeManyOpts function for SBS", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")

  indices <- grep("Liver-HCC", colnames(PCAWG7::spectra$PCAWG$SBS96))
  spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(
    mutation.type = "SBS96",
    cancer.type = "Liver-HCC"
  )
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
      use.forward.search      = TRUE
    )

    return(retval)
  }

  my.opts1 <- DefaultManyOpts()
  retval1 <- TestFunction(my.opts1)

  my.loglh.fn1 <- function(spectrum, expected.counts) {
    loglh0 <-
      stats::dmultinom(x = spectrum, prob = expected.counts, log = TRUE)
    return(loglh0)
  }
  my.opts2 <- CustomizeManyOpts(loglh.fn = my.loglh.fn1)
  retval2 <- TestFunction(my.opts2)

  expect_equal(retval1$proposed.assignment, retval2$proposed.assignment)

  my.opts3 <- DefaultManyOpts(likelihood.dist = "neg.binom")
  my.opts3$nbinom.size <- 8
  retval3 <- TestFunction(my.opts3)

  my.loglh.fn3 <- function(spectrum, expected.counts) {
    loglh0 <-
      sum(stats::dnbinom(
        x = spectrum, mu = expected.counts,
        size = 8, log = TRUE
      ))
    return(loglh0)
  }
  my.opts4 <- CustomizeManyOpts(loglh.fn = my.loglh.fn3)
  retval4 <- TestFunction(my.opts4)
  expect_equal(retval3$proposed.assignment, retval4$proposed.assignment)
})
