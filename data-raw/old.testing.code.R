# OLD testing code, not up to date
# 
TestAnySubsetPresent <- function(extra.sig.indices, eso.indices) {
  
  eso.spectra <- TestEsoSpectra(eso.indices)
  
  m.opts <- DefaultManyOpts()
  
  sigs.plus <- TestEsoSigs(extra.sig.indices)
  
  sigs.plus <- TestEsoSigs(extra.sig)
  
  # Set up variables so result can be a data frame.
  sample.id.v <- c()
  df.v <- c()
  statistic.v <- c()
  p.v <- c()
  sigs.added.v <- c()
  base.sigs.v  <- c()
  for (sample.id in colnames(eso.spectra)[idx.to.test]) {
    out <- AnySigSubsetPresent(eso.spectra[ , sample.id, drop = FALSE],
                               all.sigs = sp96.sig,
                               H0.sigs   = eso.min.sigs, 
                               more.sigs = c("SBS17a", "SBS17b"),
                               m.opts = DefaultManyOpts())
    sample.id.v  <- c(sample.id.v,  rep(sample.id, length(out)))
    sigs.added.v <- c(sigs.added.v, sapply(out, function(x) x$sigs.added))
    df.v         <- c(df.v,         sapply(out, function(x) x$df))
    statistic.v  <- c(statistic.v,  sapply(out, function(x) x$statistic))
    p.v          <- c(p.v,          sapply(out, function(x) x$p))
    base.sigs.v  <- c(base.sigs.v,  sapply(out, function(x) x$base.sigs))
  }
  
  retval <-
    data.frame(sample.id = sample.id.v,
               sigs.added = sigs.added.v,
               df         = df.v,
               statistic  = statistic.v,
               p          = p.v,
               base.sigs  = base.sigs.v)
  return(retval)
}

if (FALSE) {
  any.test <- TestAnySubsetPresent(c(1,2,6))
}
