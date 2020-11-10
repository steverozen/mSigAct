catSBS96.1 <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
my.sig.SBS96 <-
  ICAMS.shiny::CancerTypeToSigSubset(cancer.type = "Biliary-AdenoCA", tumor.cohort = "PCAWG",
                                     sig.type = "SBS96", region = "genome")
library(mSigAct)

mm <- DefaultManyOpts()
mm$trace <- 100

foo <- SparseAssignActivity(spectra = catSBS96.1, sigs = my.sig.SBS96,
                            max.level = 16,
                            p.thresh = 0.01,
                            m.opts = mm,
                            mc.cores.per.sample = 100)

