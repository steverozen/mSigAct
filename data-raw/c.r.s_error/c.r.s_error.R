
load("data-raw/c.r.s_error/liver_sample_test.Rdata")

calculate_p_thresh <- function(p_thresh_level, sig) {
  if (p_thresh_level == "mid") {
    p_thresh <- 0.05 / ncol(sig)
  } else if (p_thresh_level == "veryhigh") {
    p_thresh <- 0.1
  } else if (p_thresh_level == "low") {
    p_thresh <- 0.01 / (4 * ncol(sig))
  } else if (p_thresh_level == "vlow") {
    p_thresh <- 0.001 / (4 * ncol(sig))
  } else {
    stop("unknown p_thresh_level, ", p_thresh_level)
  }
  return(p_thresh)
}

my_opts$trace <- 1e10

test <-
  mSigAct::PresenceAssignActivity(spectra = liver_sample_to_test, 
                                  sigs = liver_sig_non_msi, 
                                  output.dir = "data-raw/c.r.s_error/test", 
                                  max.level = ncol(liver_sig_non_msi) - 1, 
                                  p.thresh = calculate_p_thresh(p_thresh_level = "vlow", sig = liver_sig_non_msi), 
                                  m.opts = my_opts, 
                                  num.parallel.samples = 1, 
                                  mc.cores.per.sample = 50, 
                                  seed = 145879,
                                  save.files = TRUE)
