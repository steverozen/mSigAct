setwd("tests.id.2/")     # <--------

library(mSigAct)

# mutation.type = "SBS96"
mutation.type = "ID"      # <------
# mutation.type = "SBS192"
# mutation.type = "DBS78"

# devtools::load_all(".")

p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
  PCAWG7::spectra$PCAWG[[mutation.type]])

cancer.types <- names(p7)

mm <- mSigAct::DefaultManyOpts()
mm$trace <- 100

# debug(OneMAPAssignTest)
# debug(MAPAssignActivity1)
for (tt in cancer.types) {
  message("cancer type = ", tt)
  mSigAct:::PCAWGMAPTest(cancer.type             = tt,
                         sample.index            = 1,
                         mutation.type           = mutation.type,
                         max.mc.cores            = 50,
                         out.dir                 = TRUE,
                         m.opts                  = mm,
                         max.level               = 100,
                         max.presence.proportion = 0.99) }

# nice R --vanilla < background.test.MAP.R &> tests.sbs192.2/bg.SBS192.x &
