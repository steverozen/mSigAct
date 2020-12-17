
library(mSigAct)
mutation.type <- "SBS96"
# mutation.type <- "ID"
# mutation.type <- "SBS192"
# mutation.type <- "DBS78"

#  ls */*/no.results.txt
#
# ls */DIFF* | wc
# devtools::load_all(".")

p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
  PCAWG7::spectra$PCAWG[[mutation.type]])

cancer.types <- names(p7)
# cancer.types <- cancer.types[-(1:7)] # starting where we left off

# Liver-HCC has 24 signatures, index 20
# Stomach-AdenoCA has 20 signatures, index 35
# Biliary-AdenoCA has 17


mm <- mSigAct::DefaultManyOpts()
# mm$trace <- 100
mm$global.opts$maxeval <- 10000

# debug(OneMAPAssignTest)
# debug(MAPAssignActivity1)

total.time <- system.time(
  for (tt in cancer.types[1]) {
    # for (ii in 1:min(ncol(p7[[tt]]),5)) {
    for (ii in 2) {
      message("sample index = ", ii)
      set.seed(101010+1, kind = "L'Ecuyer-CMRG")
      message("cancer type = ", tt)
      xx <- mSigAct::YPCAWGMAPTest(
        cancer.type             = tt,
        sample.index            = ii,
        mutation.type           = mutation.type,
        max.mc.cores            = 50,
        out.dir                 = paste(tt, mutation.type, ii, sep = "-"),
        m.opts                  = mm,
        max.level               = 100,
        max.presence.proportion = 0.99) }}

)

message("TOTAL TIME")
print(total.time)

# cd tests.id.gobal.eval.1000     # <--------
# nice R --vanilla < ~/mSigAct/data-raw/Ybackground.test.MAP.R &> log.txt &

