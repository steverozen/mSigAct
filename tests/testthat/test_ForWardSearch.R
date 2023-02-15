test_that("Use forward search for MAPAssignActivity", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  spectra <- PCAWG7::spectra$PCAWG$SBS96
  catalog.list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  biliary.catalogs <- catalog.list$`Biliary-AdenoCA`
  
  sample.index <- 1
  catalog <- biliary.catalogs[, sample.index, drop = FALSE]
  SBS.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  mutation.type <- "SBS96"
  cancer.type <- "Biliary-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- SBS.sigs[, names(sigs.prop), drop = FALSE]
  my.opts <- DefaultManyOpts()
  
  retval1 <- MAPAssignActivity1(
    spect                   = catalog,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.001 / 4 * ncol(sigs),
    m.opts                  = my.opts,
    max.mc.cores            = 60,
    seed                    = 2351,
    use.sparse.assign       = TRUE, 
    use.sig.presence.test   = TRUE,
    use.forward.search      = TRUE, 
    sig.pres.test.nbinom.size = 3)
  
  expect_equal(nrow(retval1$proposed.assignment), 5)
  
  
  catalogs <- biliary.catalogs[, 1:3, drop = FALSE]
  
  output.dir <- file.path(tempdir(), "test")
  retval2 <- MAPAssignActivity(
    spectra                 = catalogs,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.001 / 4 * ncol(sigs),
    m.opts                  = my.opts,
    num.parallel.samples    = 3,
    mc.cores.per.sample     = 20,
    seed                    = 2351,
    use.sparse.assign       = TRUE, 
    use.sig.presence.test   = TRUE,
    use.forward.search      = TRUE, 
    sig.pres.test.nbinom.size = 3)
  expect_equal(retval1$proposed.assignment,
               RemoveZeroActivitySig(retval2$proposed.assignment[, 1, drop = FALSE]))
  
  sapply(output.dir, FUN = unlink, recursive = TRUE)
})
