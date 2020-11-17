setwd("tests.sbs.gobal.eval.1000")     # <--------

library(mSigAct)

mutation.type = "SBS96"
# mutation.type = "ID"      # <------
# mutation.type = "SBS192"
# mutation.type = "DBS78"

# devtools::load_all(".")

p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
  PCAWG7::spectra$PCAWG[[mutation.type]])

cancer.types <- names(p7)

mm <- mSigAct::DefaultManyOpts()
mm$trace <- 100
mm$global.opts$maxeval <- 1000

# debug(OneMAPAssignTest)
# debug(MAPAssignActivity1)

total.time <- system.time(
for (tt in cancer.types) {
  message("cancer type = ", tt)
  mSigAct:::PCAWGMAPTest(cancer.type             = tt,
                         sample.index            = 1,
                         mutation.type           = mutation.type,
                         max.mc.cores            = 50,
                         out.dir                 = TRUE,
                         m.opts                  = mm,
                         max.level               = 100,
                         max.presence.proportion = 0.99) })

print(total.time)

# nice R --vanilla < background.test.MAP.R &> tests.sbs.gobal.eval.1000/log.txt &

