test_that("Edge case when spectrum have very low mutation counts", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")

  spectra <- PCAWG7::spectra$PCAWG$DBS78
  breast_samples <-
    c(
      "Breast-AdenoCA::SP117244", "Breast-AdenoCA::SP118073",
      "Breast-AdenoCA::SP124195", "Breast-AdenoCA::SP117850", 
      "Breast-AdenoCA::SP117593", "Breast-AdenoCA::SP118077",
      "Breast-AdenoCA::SP117728", "Breast-AdenoCA::SP4535"
    )
  catalog <- spectra[, breast_samples, drop = FALSE]
  colSums(catalog)
  
  DBS.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
  mutation.type <- "DBS78"
  cancer.type <- "Breast-AdenoCA"
  sigs.prop <- ExposureProportions(
    mutation.type = mutation.type,
    cancer.type = cancer.type
  )
  sigs <- DBS.sigs[, names(sigs.prop), drop = FALSE]

  my.opts <- DefaultManyOpts()

  retval1 <- PresenceAssignActivity(
    spectra                 = catalog,
    sigs                    = sigs,
    output.dir              = tempdir(),
    p.thresh                = 0.05 / ncol(sigs),
    m.opts                  = my.opts,
    num.parallel.samples    = 8,
    mc.cores.per.sample     = 7,
    seed                    = 2351,
    drop.low.mut.samples    = FALSE
  )

  expect_true(all(retval1$error.messages == ""))
})
