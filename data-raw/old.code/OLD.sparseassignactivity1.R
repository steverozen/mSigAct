#' Component of \code{\link{SparseAssignActivity}} for one spectrum.
#' @keywords internal

SparseAssignActivity1OLD <- function(
  spect, sigs, max.level = 5, p.thresh = 0.05, eval_f, m.opts) {
  mode(spect) <-  'numeric'
  start <- one.lh.and.exp(spect, 
                          sigs,
                          eval_f  = eval_f,
                          m.opts  = m.opts)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp > 0.5)
  if (m.opts$trace > 0) {
    message('Starting with',
            paste(names(start.exp)[non.0.exp.index], collapse = ","),
            '\n')
    message('max.level =', max.level, '\n')
  }
  if (length(non.0.exp.index) < 2) {
    if (m.opts$trace > 0)
      message("returning, only ", length(non.0.exp.index),
              " non-0 exposures")
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)
  
  c.r.s <- sets::set() # subsets of the signature indices that cannot be removed
  
  max.p <- numeric(max.level)
  best.subset <- non.0.exp.index
  best.try <- start
  best.try.is.start <- T
  
  for (df in 1:max.level) {
    if (m.opts$trace > 1) message('df =', df, '\n')
    subsets <- sets::set_combn(non.0.exp.index, df)
    for (subset in subsets) {
      # subset is of class set (package sets)
      if (is.superset.of.any(subset, c.r.s)) next
      subset.to.remove.v <- as.numeric(subset)
      subset.name <- paste(names(start.exp)[subset.to.remove.v], collapse = ",")
      tmp.set <- setdiff(non.0.exp.index, subset.to.remove.v)
      try.sigs <- sigs[ , tmp.set]
      
      if (length(tmp.set) == 1) {
        try.sigs <- as.matrix(try.sigs)
        rownames(try.sigs) <- rownames(sigs)
        if (m.opts$trace > 0) message("Running new code\n")
      }
      
      # Get the max lh for try.sig / Get the maximum likelihood 
      # reconstruction using try.sig
      try <- one.lh.and.exp(spect, 
                            try.sigs,
                            eval_f    = eval_f,
                            m.opts = m.opts)
      # try contains maxlh and exposure
      statistic <- 2 * (lh.w.all - try$loglh)
      chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
      if (m.opts$trace > 0) {
        message("Trying" , subset.name, ", p = ", chisq.p)
        if (m.opts$trace > 1) 
          message("loglh = ", try$loglh, "; statistic  =", statistic)
      }
      if (chisq.p > p.thresh) {
        # This an acceptable set of exposures
        if (m.opts$trace > 0) message("acceptable")
        if (chisq.p > max.p[df]) {
          # This is the best so far
          max.p[df] <- chisq.p
          best.subset <- tmp.set
          best.try <- try
          best.try.is.start <- F
          if (m.opts$trace > 0) message('best\n')
        } else {
          if (m.opts$trace > 0) message('not best\n')
        }
      } else {
        c.r.s <- sets::set_union(c.r.s, sets::set(subset))
        if (m.opts$trace > 0)  {
          message('cannot remove\n')
        }
      }
    } # end for (subset in subsets)
    if (max.p[df] == 0) break
  } # end for df in 1:max.level
  
  # Need to return the exposures in the context of the 
  # orginal signatures matrix
  out.exp <- numeric(ncol(sigs)) #all zeros
  names(out.exp) <- colnames(sigs)
  if (best.try.is.start) {
    out.exp <- start.exp
  } else {
    out.exp[best.subset] <- best.try$exposure
  }
  if (m.opts$trace > 0) {
    message('max.p =', paste(max.p, collapse = ', '), '\n')
    message('Ending with',
            paste(names(start.exp)[best.subset], collapse = ","),
            '\n')
  }
  stopifnot(abs(sum(out.exp) - sum(spect)) < 1)
  return(out.exp)
}
