load("data-raw/loglh_two_model/test_loglh_two_model.Rdata")
library(ICAMS)

retval <- CompareAndPlotLoglh(spectrum = spectrum,
                              exposure = exposure,
                              sigs = sigs,
                              sig.to.test = "DBS2",
                              neg.binom.size = 8, 
                              file = "data-raw/loglh_two_model/DBS2_loglh.pdf")


new_sigs <- sigs[, rownames(exposure)]
my_opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
my_opts$nbinom.size <- 8

test <- SignaturePresenceTest1(spectrum = spectrum,
                               sigs = new_sigs,
                               target.sig.index = "DBS2",
                               m.opts = my_opts,
                               seed = 8157)



skin_samples <- 
  readRDS("data-raw/loglh_two_model/dbs2_skin_samples_paa_only_size_8.Rds")
dbs_spectra <- PCAWG7::spectra$PCAWG$DBS78
indices <- grep(pattern = "Skin", x = colnames(dbs_spectra))
skin_spectra <- dbs_spectra[, indices]

View(skin_spectra)
peak_counts <- skin_spectra["CCAA", ]

load("data-raw/loglh_two_model/test_loglh_two_model.Rdata")
library(ICAMS)


skin_sig_prop <-
  mSigAct::ExposureProportions(
    mutation.type = "DBS78",
    cancer.type = "Skin-Melanoma"
  )
skin_sig_names <- names(skin_sig_prop)

skin_sample_to_test <- spectrum
skin_sig_non_msi <- sigs[, skin_sig_names]

p_thresh_level <- "vlow"
#p_thresh_level <- "denom_10_pw_20"
#p_thresh_level <- "denom_10_pw_4"


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
my_opts$trace <- 1e13

nb_loglh_thresh <- 11
Sys.setenv(NB_LOGLH_THRESH = nb_loglh_thresh)

paa_retval <-
  PresenceAssignActivity(
    spectra = skin_spectra,
    sigs = skin_sig_non_msi,
    output.dir = file.path(
      "data-raw/loglh_two_model/paa",
      paste0("p_level_", p_thresh_level, "_nb_loglh_thresh_", nb_loglh_thresh)
    ),
    max.level = ncol(skin_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = skin_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 30,
    mc.cores.per.sample = 5,
    seed = 145879,
    save.files = FALSE,
    drop.low.mut.samples = FALSE,
    use.forward.search = TRUE
  )
saveRDS(paa_retval, file = "data-raw/loglh_two_model/paa/nb_loglh_thresh_11.Rds")
paa_exp <- paa_retval$proposed.assignment
dbs2_exp <- paa_exp["DBS2", ]
dbs2_exp_non_zero <- dbs2_exp[dbs2_exp > 0]
length(dbs2_exp_non_zero)

nb_loglh_thresh <- 20
Sys.setenv(NB_LOGLH_THRESH = nb_loglh_thresh)

paa_retval_20 <-
  PresenceAssignActivity(
    spectra = skin_spectra,
    sigs = skin_sig_non_msi,
    output.dir = file.path(
      "data-raw/loglh_two_model/paa",
      paste0("p_level_", p_thresh_level, "_nb_loglh_thresh_", nb_loglh_thresh)
    ),
    max.level = ncol(skin_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = skin_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 30,
    mc.cores.per.sample = 5,
    seed = 145879,
    save.files = FALSE,
    use.forward.search = TRUE
  )
paa_exp_20 <- paa_retval_20$proposed.assignment
dbs2_exp_20 <- paa_exp_20["DBS2", ]
dbs2_exp_20_non_zero <- dbs2_exp_20[dbs2_exp_20 > 0]
length(dbs2_exp_20_non_zero)


nb_loglh_thresh <- 30
Sys.setenv(NB_LOGLH_THRESH = nb_loglh_thresh)

paa_retval_30 <-
  PresenceAssignActivity(
    spectra = skin_spectra,
    sigs = skin_sig_non_msi,
    output.dir = file.path(
      "data-raw/loglh_two_model/paa",
      paste0("p_level_", p_thresh_level, "_nb_loglh_thresh_", nb_loglh_thresh)
    ),
    max.level = ncol(skin_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = skin_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 30,
    mc.cores.per.sample = 5,
    seed = 145879,
    save.files = FALSE,
    use.forward.search = TRUE
  )
paa_exp_30 <- paa_retval_30$proposed.assignment
dbs2_exp_30 <- paa_exp_30["DBS2", ]
dbs2_exp_30_non_zero <- dbs2_exp_30[dbs2_exp_30 > 0]
length(dbs2_exp_30_non_zero)


test1 <-
  PresenceAssignActivity(
    spectra = skin_sample_to_test,
    sigs = skin_sig_non_msi,
    output.dir = file.path(
      "data-raw/loglh_two_model/deep_dive",
      paste0("p_level_", p_thresh_level, "_nb_loglh_thresh_", nb_loglh_thresh)
    ),
    max.level = ncol(skin_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = skin_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 1,
    mc.cores.per.sample = 10,
    seed = 145879,
    save.files = FALSE,
    use.forward.search = TRUE
  )

