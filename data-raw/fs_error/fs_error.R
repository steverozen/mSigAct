load("data-raw/fs_error/fs_error.Rdata")

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

p_thresh_level <- "denom_10_pw_4"

my_opts$nbinom.size <- 8
size <- 8

test2 <-
  mSigAct::PresenceAssignActivity(
    spectra = kidney_sample_to_test,
    sigs = kidney_sig_non_msi,
    output.dir = file.path(
      "data-raw/fs_error/deep_dive",
      paste0("size_", size)
    ),
    max.level = ncol(kidney_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = kidney_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 1,
    mc.cores.per.sample = 50,
    seed = 145879,
    save.files = TRUE,
    use.forward.search = TRUE
  )
