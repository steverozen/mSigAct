context("Use multinomial distribution to calculate likelihood")

test_that("Use multinomial distribution to calculate likelihood", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  prior.prop <- ExposureProportions("SBS96", cancer.type = "Skin-Melanoma")
  prior.prop1 <- 
    prior.prop[!names(prior.prop) %in% cosmicsig::possible_artifacts()]
  
  skin.spectra <- 
    PCAWG7::spectra$other.genome$SBS96[, "Skin-Melanoma::ML_31_T_01", drop = FALSE]
  
  # Use negative binomial distribution to calculate log likelihood
  MAP.out1 <- 
    MAPAssignActivity(spectra = skin.spectra,
                      sigs = cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96,
                      sigs.presence.prop = prior.prop1,
                      max.level = length(prior.prop1) - 1,
                      m.opts = DefaultManyOpts(likelihood.dist = "neg.binom"),
                      mc.cores.per.sample = 30,
                      max.subsets = 1e10,
                      seed = 8386,
                      output.dir = file.path(tempdir(), "skin.mapout1"))
  
  expect_equal(MAP.out1$reconstruction.distances$proposed.assignment$cosine,
               0.9517956, tolerance = 1e-2)
  
  # Use multinomial distribution to calculate log likelihood
  my.opts <- DefaultManyOpts(likelihood.dist = "multinom")
  MAP.out2 <- 
    MAPAssignActivity(spectra = skin.spectra,
                      sigs = cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96,
                      sigs.presence.prop = prior.prop1,
                      m.opts = my.opts,
                      max.level = length(prior.prop1) - 1,
                      mc.cores.per.sample = 30,
                      max.subsets = 1e10,
                      seed = 8386,
                      output.dir = file.path(tempdir(), "skin.mapout2"))
  expect_equal(MAP.out2$reconstruction.distances$proposed.assignment$cosine,
               0.9947313, tolerance = 1e-5)
  
  unlink(file.path(tempdir(), "skin.mapout1"), recursive = TRUE)
  unlink(file.path(tempdir(), "skin.mapout2"), recursive = TRUE)
})