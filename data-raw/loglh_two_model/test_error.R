load("data-raw/loglh_two_model/paa_error.Rdata")

calculate_p_thresh <- function(p_thresh_level, sig) {
  if (p_thresh_level == "mid") {
    p_thresh <- 0.05 / ncol(sig)
  } else if (p_thresh_level == "veryhigh") {
    p_thresh <- 0.1
  } else if (p_thresh_level == "low") {
    p_thresh <- 0.01 / (4 * ncol(sig))
  } else if (p_thresh_level == "vlow") {
    p_thresh <- 0.001 / (4 * ncol(sig))
  } else if (grepl(pattern = "denom", x = p_thresh_level)){
    all_numbers <- 
      stringi::stri_extract_all(p_thresh_level, regex="\\d+")[[1]]
    base_num <- as.numeric(all_numbers[1])
    pow <- as.numeric(all_numbers[2])
    p_thresh <- 0.05 / ncol(sig) / (base_num ^ pow)
  } else {
    stop("unknown p_thresh_level, ", p_thresh_level)
  }
  return(p_thresh)
}
Sys.setenv(NB_LOGLH_THRESH = 30)
my_opts$trace <- 1e10
retval <- PresenceAssignActivity(
  spectra                   = sample_test,
  sigs                      = skin_sigs,
  output.dir                = "data-raw/loglh_two_model/test_error",
  max.level                 = ncol(skin_sigs) - 1,
  p.thresh                  = calculate_p_thresh(p_thresh_level, skin_sigs),
  m.opts                    = my_opts,
  num.parallel.samples      = 1,
  mc.cores.per.sample       = 5,
  seed                      = 145879,
  drop.low.mut.samples      = FALSE,
  save.files                = TRUE,
  use.forward.search        = TRUE
)
