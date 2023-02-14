spect <<- spect
#sigs <<- remained.sigs
m.opts <<- m.opts
max.mc.cores <<- max.mc.cores
p.thresh <<- 
step <- 1


ForwardSearch <- function(spect, sigs, m.opts, max.mc.cores, p.thresh) {
  start <- OptimizeExposure(spectrum = spect, sigs = sigs, m.opts  = m.opts)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  num_sig <- ncol(sigs)
  
  if (m.opts$trace >= 0) {
    message("Performing forward search: ")
  }
  
  sigs.to.choose <- sigs
  optimal.sigs <- sigs[, 0, drop = FALSE]

  for (step in seq_len(ncol(sigs))) {
    browser()
    if (m.opts$trace >= 0) {
      message("Forward search step ", step, ": ")
    }
    
    retval <- 
      parallel::mclapply(1:ncol(sigs.to.choose), FUN = function(i) {
        one.sig <- sigs.to.choose[, i, drop = FALSE]
        optimal.sigs <- cbind(optimal.sigs, one.sig)
        
        try.exp <- OptimizeExposure(spectrum = spect, 
                                    sigs = optimal.sigs, m.opts = m.opts)
        if (is.infinite(try.exp$loglh)) {
          # It was not possible to generate the spectrum from the signatures, so
          # try.exp$loglh should be -Inf and the p value the test that the spectrum
          # can be better reconstructed with sigs than with optimal.sigs is 0.
          if (try.exp$loglh < 0) {
            chisq.p <- 0
          } else {
            stop("Probable coding error")
          }
        } else {
          statistic <- 2 * (lh.w.all - try.exp$loglh)
          chisq.p <- stats::pchisq(q = statistic, 
                                   df = num_sig - step, lower.tail = FALSE)
        }
        
        return(list(p.value = chisq.p, exposure = try.exp$exposure))
        
      }, mc.cores = max.mc.cores)
    names(retval) <- colnames(sigs.to.choose)
    
    p.values <- sapply(retval, FUN = "[[", "p.value")
    
    # Find out the largest p-value
    index <- which.max(p.values)
    largest.p.value <- p.values[index]
    selected.sig <- sigs.to.choose[, index, drop = FALSE]
    optimal.sigs <- cbind(optimal.sigs, selected.sig)
    optimal.exposure <- as.matrix(retval[[index]][["exposure"]])
    colnames(optimal.exposure) <- colnames(spect)
    
    if (largest.p.value > p.thresh) {
      return(optimal.exposure = optimal.exposure)
    }
    
    sigs.to.choose <- sigs.to.choose[, -index, drop = FALSE]
    
  }
}

ForwardSearch(spect = spect, sigs = sigs, m.opts = m.opts,
              max.mc.cores = max.mc.cores, p.thresh = p.thresh)