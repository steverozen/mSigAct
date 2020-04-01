# Next variable not used.
pre.sa <- c("SBS1","SBS5", "SBS13", 
            "SBS15", "SBS17a", "SBS17b", 
            "SBS18", "SBS28", "SBS30",
            "SBS37", "SBS44", "SBS45",
            "SBS51", "SBS54", "SBS58")

e <- local({load("data-raw/arnoud.error.2020.03.31.RData"); environment()})
ct <- e$catTmp

library(BSgenome.Hsapiens.1000genomes.hs37d5)
# The next line silences a warning due to attempted lookup of a file 
# path on Arnoud's computer
attr(ct, "ref.genome") <- BSgenome.Hsapiens.1000genomes.hs37d5

# These are the H0 signature names.
post.sa <- c("SBS1", "SBS5", "SBS13", "SBS18", "SBS45")

retval <- mSigAct::AnySigSubsetPresent(
  spect           = ct,
  all.sigs        = PCAWG7::signature$genome$SBS96[ , c("SBS30", post.sa)],
  Ha.sigs.indices = 1,
  
  # The eval_f is the key difference between this test and
  # the next.
  eval_f          = mSigAct::ObjFnBinomMaxLHNoRoundOK,
  m.opts          = mSigAct::DefaultManyOpts(),
  max.mc.cores    = 1)

retval$all.Ha.info[[1]][c("p", "sigs.added", "base.sigs")]
# $p
# [1] 0.0401562
# 
# $sigs.added
# [1] "SBS30"
#
# $base.sigs
# [1] "SBS1,SBS5,SBS13,SBS18,SBS45"


retval$all.Ha.info[[1]]$Ha.info[c("loglh", "exposure")]
# $loglh
# [1] -140.6682
# 
# $exposure
# SBS30      SBS1      SBS5     SBS13     SBS18     SBS45 
# 16.169331 50.917434 79.629568  5.092431 13.116205  7.075032 


retval.nr <- mSigAct::AnySigSubsetPresent(
  spect           = ct,
  all.sigs        = PCAWG7::signature$genome$SBS96[ , c("SBS30", post.sa)],
  Ha.sigs.indices = 1,
  eval_f          = mSigAct::ObjFnBinomMaxLHMustRound,
  m.opts          = mSigAct::DefaultManyOpts(),
  max.mc.cores    = 1)

retval.nr$all.Ha.info[[1]][c("p", "sigs.added", "base.sigs")]
# $p
# [1] NaN
# 
# $sigs.added
# [1] "SBS30"
#
# $base.sigs
# [1] "SBS1,SBS5,SBS13,SBS18,SBS45"

retval.nr$all.Ha.info[[1]]$Ha.info[c("loglh", "exposure")]
# $loglh
# [1] -Inf
#
# $exposure
# SBS30     SBS1     SBS5    SBS13    SBS18    SBS45 
# 28.66667 28.66667 28.66667 28.66667 28.66667 28.66667 

