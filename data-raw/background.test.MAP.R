
library(mSigAct)

# mutation.type = "SBS96"
# mutation.type = "ID"      # <------
mutation.type = "SBS192"
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

for (tt in cancer.types[2]) {
  set.seed(101010+1, kind = "L'Ecuyer-CMRG")
  message("cancer type = ", tt)
  mSigAct:::PCAWGMAPTest(cancer.type             = tt,
                         sample.index            = 1,
                         mutation.type           = mutation.type,
                         max.mc.cores            = 50,
                         out.dir                 = "foo.long", # TRUE,
                         m.opts                  = mm,
                         max.level               = 100,
                         max.presence.proportion = 0.99) }

)

print(total.time)

# cd tests.id.gobal.eval.1000     # <--------
# nice R --vanilla < ../background.test.MAP.R &> log.txt &

