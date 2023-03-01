BackwardSearch <- 
  function(spect, sigs, m.opts, max.mc.cores, p.thresh, optimal.exposure) {
    optimal.sigs <- sigs[, names(optimal.exposure), drop = FALSE]
    num_sig <- ncol(optimal.sigs)
    
    if (num_sig == 1) {
      return(optimal.exposure = optimal.exposure)
    }
    
    if (m.opts$trace >= 0) {
      message("\nStarting backward search for ", colnames(spect))
    }
    start <- 
      OptimizeExposure(spectrum = spect, sigs = optimal.sigs, m.opts = m.opts)
    lh.w.all <- start$loglh # The likelihood with all signatures
    
    if (m.opts$trace >= 1) {
      message("log likelihood using all signatures = ", lh.w.all)
    }
    

    
    sigs.to.test <- optimal.sigs
    
    
}


TwoWaySearch <- function(spect, sigs, m.opts, max.mc.cores, p.thresh) {
  optimal.exposure <-
    ForwardSearch(spect = spect, sigs = sigs,
                  m.opts = m.opts, max.mc.cores = max.mc.cores,
                  p.thresh = p.thresh)
  
  
  
}