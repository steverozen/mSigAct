dbs_spectra <- PCAWG7::spectra$PCAWG$DBS78
indices <- grep(pattern = "Skin", x = colnames(dbs_spectra))
skin_spectra <- dbs_spectra[, indices]

load("data-raw/loglh_two_model/test_loglh_two_model.Rdata")
library(ICAMS)

skin_sig_prop <-
  mSigAct::ExposureProportions(
    mutation.type = "DBS78",
    cancer.type = "Skin-Melanoma"
  )
skin_sig_names <- names(skin_sig_prop)

skin_sig_non_msi <- sigs[, skin_sig_names]

p_thresh_level <- "vlow"

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

my_opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
my_opts$nbinom.size <- 8

nb_loglh_thresh <- 11
Sys.setenv(NB_LOGLH_THRESH = nb_loglh_thresh)

paa_retval <-
  PresenceAssignActivity(
    spectra = skin_spectra,
    sigs = skin_sig_non_msi,
    output.dir = "data-raw/loglh_two_model/paa/skin_inconsistency",
    max.level = ncol(skin_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = skin_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 30,
    mc.cores.per.sample = 5,
    seed = 145879,
    drop.low.mut.samples = FALSE,
    save.files = FALSE,
    use.forward.search = TRUE
  )
saveRDS(paa_retval, 
        file = "data-raw/loglh_two_model/paa/skin_inconsistency/nb_loglh_thresh_11.Rds")
paa_exp <- paa_retval$proposed.assignment
dbs2_exp <- paa_exp["DBS2", ]
dbs2_exp_non_zero <- dbs2_exp[dbs2_exp > 0]
length(dbs2_exp_non_zero)
