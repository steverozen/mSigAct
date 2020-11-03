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

m.opts <- DefaultManyOpts()

library(profvis)

ww <- profvis(

   SA.out <- SparseAssignActivity1(spect       = spect,
                                  sigs         = some.sigs,
                                  eval_f       = ObjFnBinomMaxLHMustRound,
                                  m.opts       = m.opts,
                                  max.level            = 1,
                                  max.mc.cores         = 1,
                                  num.parallel.samples = 1
  ))


htmlwidgets::saveWidget(ww, "0profile.ww.html")



