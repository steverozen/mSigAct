context("Use multinomial distribution to calculate likelihood")

test_that("Use multinomial distribution to calculate likelihood", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  prior.prop <- ExposureProportions("SBS96", cancer.type = "Skin-Melanoma")
  prior.prop1 <- prior.prop[!names(prior.prop) %in% PossibleArtifacts()]
  
  skin.spectra <- 
    PCAWG7::spectra$other.genome$SBS96[, "Skin-Melanoma::ML_31_T_01", drop = FALSE]
  
  # Use negative binomial distribution to calculate log likelihood
  my.opts1 <- DefaultManyOpts(likelihood.dist = "neg.binom")
  MAP.out1 <- 
    MAPAssignActivity(spectra = skin.spectra,
                      sigs = PCAWG7::signature$genome$SBS96,
                      sigs.presence.prop = prior.prop1,
                      max.level = length(prior.prop1) - 1,
                      m.opts = my.opts1,
                      mc.cores.per.sample = 30,
                      max.subsets = 1e10,
                      output.dir = file.path(tempdir(), "skin.mapout1"))
  
  expect_equal(MAP.out1$reconstruction.distances$MAP.distances$cosine,
               0.9517956)
  
  # Use multinomial distribution to calculate log likelihood
  my.opts2 <- DefaultManyOpts()
  my.opts2$trace <- 100
  MAP.out2 <- 
    MAPAssignActivity(spectra = skin.spectra,
                      sigs = PCAWG7::signature$genome$SBS96,
                      sigs.presence.prop = prior.prop1,
                      max.level = length(prior.prop1) - 1,
                      m.opts = my.opts2,
                      mc.cores.per.sample = 30,
                      max.subsets = 1000000,
                      output.dir = file.path(tempdir(), "skin.mapout2"))
  unlink(file.path(tempdir(), "skin.mapout1"))
  
})