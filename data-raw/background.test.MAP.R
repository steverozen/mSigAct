library(mSigAct)

mutation.type = "SBS96"
# mutation.type = "ID"
# mutation.type = "SBS192"
# mutation.type = "DBS78"

# devtools::load_all(".")

p7 <- PCAWG7::SplitPCAWGMatrixByTumorType(
  PCAWG7::spectra$PCAWG[[mutation.type]])

cancer.types <- names(p7)
# cancer.types <- "Cervix-AdenoCA"

# debug(OneMAPAssignTest)
# debug(MAPAssignActivity1)
for (tt in cancer.types) {
  message("cancer type = ", tt)
  mSigAct:::PCAWGMAPTest(cancer.type = tt,
                         sample.index = 1,
                         mutation.type = mutation.type,
                         max.mc.cores = 100,
                         out.dir = TRUE,
                         max.level = 100,
                         reduce.must.have = 0.99) }

