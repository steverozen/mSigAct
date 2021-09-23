test_that("DropLowMutationSamples for ID spectra", {
  spectra <- PCAWG7::spectra$PCAWG$ID
  catalog.list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  lung.catalogs <- catalog.list$`Lung-AdenoCA`
  
  sample.index <- 6
  catalog <- lung.catalogs[, sample.index, drop = FALSE]
  ID.sigs <- PCAWG7::signature$genome$ID
  mutation.type <- "ID"
  cancer.type <- "Lung-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- ID.sigs[, names(sigs.prop), drop = FALSE]
  my.opts <- DefaultManyOpts()
  
  output.dir <- file.path(tempdir(), paste0("test", 1))
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
  
  sapply(output.dir, FUN = unlink, recursive = TRUE)
})