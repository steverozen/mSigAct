
MakeAndSaveHepG2Tests <- function(num.replicates = 20) {
  
  sig.names.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
  contribution <- c(0.1, 0.5, 0.9)
  
  HepG2.bg.tests.no.noise <-
    MakeSynTestGrid(sig.names.to.test = sig.names.to.test,
                    contribution      = contribution,
                    bg.sig.info       = mSigAct::HepG2.background.info,
                    num.replicates    = num.replicates)
  usethis::use_data(HepG2.bg.tests.no.noise, overwrite = TRUE)
}
