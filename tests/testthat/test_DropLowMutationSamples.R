test_that("DropLowMutationSamples for ID spectra", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  spectra <- PCAWG7::spectra$PCAWG$ID
  catalog.list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  lung.catalogs <- catalog.list$`Lung-AdenoCA`
  
  sample.index <- 6
  catalog <- lung.catalogs[, sample.index, drop = FALSE]
  ID.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  mutation.type <- "ID"
  cancer.type <- "Lung-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- ID.sigs[, names(sigs.prop), drop = FALSE]
  my.opts <- DefaultManyOpts()
  
  output.dir <- file.path(tempdir(), paste0("test", 1:3))
  retval1 <- expect_message(
    MAPAssignActivity(
    spectra                 = catalog,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir[1],
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 1,
    mc.cores.per.sample     = 30,
    seed                    = 2351),
    "Samples with total mutations less than 25 were excluded in the analysis for ID")
  
  retval2 <- expect_message(
    MAPAssignActivity1(
    spect                   = catalog,
    sigs                    = sigs,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    max.mc.cores            = 30,
    seed                    = 2351,
    use.sparse.assign       = TRUE),
    "Samples with total mutations less than 25 were excluded in the analysis for ID")
  
  
  retval3 <- MAPAssignActivity(
    spectra                 = catalog,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir[2],
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 1,
    mc.cores.per.sample     = 30,
    seed                    = 2351,
    drop.low.mut.samples    = FALSE)
  
  expect_equal(nrow(retval3$proposed.assignment), 3)
  
  retval4 <- MAPAssignActivity(
    spectra                 = catalog,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir[3],
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    m.opts                  = my.opts,
    num.parallel.samples    = 1,
    mc.cores.per.sample     = 30,
    seed                    = 2351,
    use.sparse.assign       = TRUE,
    drop.low.mut.samples    = FALSE)
  
  expect_equal(nrow(retval4$proposed.assignment), 2)
  
  sapply(output.dir, FUN = unlink, recursive = TRUE)
})