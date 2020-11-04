set.seed(101010, kind = "L'Ecuyer-CMRG")


some.sigs  <-
  PCAWG7::signature$genome$SBS96[ , 1:2, drop = FALSE]
ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
region     <- attr(some.sigs, "region", exact = TRUE)
if (is.null(region)) {
  message("Null region, why?")
  region <- "genome"
}

spect <- round(some.sigs %*% c(100,100))
spect <-
  ICAMS::as.catalog(
    spect,
    ref.genome   = ref.genome,
    region       = region,
    catalog.type = "counts")

library(mSigAct)
sessionInfo()

m.opts <- mSigAct::DefaultManyOpts()

library(profvis)

ww <- profvis::profvis(
   mSigAct:::SparseAssignActivity1(spect       = spect,
                                  sigs         = some.sigs,
                                  eval_f       = ObjFnBinomMaxLHMustRound,
                                  m.opts       = m.opts,
                                  max.level            = 1,
                                  max.mc.cores         = 1),
    prof_output = "out.Rprof"
 )


htmlwidgets::saveWidget(ww, "00profile.ww.html")



