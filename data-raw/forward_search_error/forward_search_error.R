load("data-raw/forward_search_error/forward_search_error.Rdata")

calculate_p_thresh <- function(p_thresh_level, sig) {
  if (p_thresh_level == "mid") {
    p_thresh <- 0.05 / ncol(sig)
  } else if (p_thresh_level == "veryhigh") {
    p_thresh <- 0.1
  } else if (p_thresh_level == "low") {
    p_thresh <- 0.01 / (4 * ncol(sig))
  } else if (p_thresh_level == "vlow") {
    p_thresh <- 0.001 / (4 * ncol(sig))
  } else if (grepl(pattern = "denom", x = p_thresh_level)) {
    all_numbers <-
      stringi::stri_extract_all(p_thresh_level, regex = "\\d+")[[1]]
    base_num <- as.numeric(all_numbers[1])
    pow <- as.numeric(all_numbers[2])
    p_thresh <- 0.05 / ncol(sig) / (base_num^pow)
  } else {
    stop("unknown p_thresh_level, ", p_thresh_level)
  }
  return(p_thresh)
}

compare_and_plot_loglh <-
  function(spectrum, exposure, sigs, sig_to_test, file, size) {
    my_opts <- mSigAct::DefaultManyOpts(likelihood.dist = "neg.binom")
    my_opts$nbinom.size <- size

    sigs <- sigs[, rownames(exposure), drop = FALSE]
    test_sig_spectra <- sigs[, sig_to_test, drop = FALSE]

    original_recon <- mSigAct::ReconstructSpectrum(
      sigs = sigs,
      exp = exposure,
      use.sig.names = TRUE
    )
    original_loglh <- stats::dnbinom(
      x = spectrum, mu = original_recon, size = size, log = TRUE
    )


    new_sigs <- sigs[, colnames(sigs) != sig_to_test, drop = FALSE]
    optimized_exposure <-
      mSigAct:::OptimizeExposure(
        spectrum = spectrum,
        sigs = new_sigs,
        m.opts = my_opts
      )
    new_exposure <- optimized_exposure$exposure

    new_recon <- mSigAct::ReconstructSpectrum(
      sigs = sigs,
      exp = new_exposure,
      use.sig.names = TRUE
    )

    new_loglh <- stats::dnbinom(
      x = spectrum, mu = new_recon, size = size, log = TRUE
    )

    original_minus_loglh <- original_loglh * -1
    new_minus_loglh <- new_loglh * -1
    substracted_loglh <- new_minus_loglh - original_minus_loglh

    to_plot <- cbind(
      new_minus_loglh, original_minus_loglh,
      substracted_loglh
    )

    sample_name <- colnames(spectrum)
    colnames(to_plot) <-
      c(
        paste0(
          sample_name,
          " (negative loglh without ", sig_to_test, ": ",
          round(colSums(new_minus_loglh), 2), ")"
        ),
        paste0(
          sample_name, " (negative loglh with ", sig_to_test, ": ",
          round(colSums(original_minus_loglh), 2), ")"
        ),
        paste0(
          sample_name, " (subtracted negative loglh: ",
          round(colSums(substracted_loglh), 2), ")"
        )
      )

    list_of_catalogs <- list(to_plot, test_sig_spectra)

    mSigAct:::PlotListOfCatalogsToPdf(
      list.of.catalogs = list_of_catalogs,
      file = file
    )
    return(to_plot)
  }

p_thresh_level <- "vlow"

retval <-
  mSigAct::PresenceAssignActivity(
    spectra = liver_sample_to_test,
    sigs = liver_sig_non_msi,
    output.dir = file.path(
      "data-raw/forward_search_error"
    ),
    max.level = ncol(liver_sig_non_msi) - 1,
    p.thresh = calculate_p_thresh(
      p_thresh_level = p_thresh_level,
      sig = liver_sig_non_msi
    ),
    m.opts = my_opts,
    num.parallel.samples = 1,
    mc.cores.per.sample = 50,
    seed = 145879,
    save.files = TRUE,
    use.forward.search = TRUE
  )

retval2 <-
  compare_and_plot_loglh(
    spectrum = liver_sample_to_test,
    exposure = retval$proposed.assignment,
    sigs = liver_sig_non_msi,
    sig_to_test = "SBS8",
    file = "data-raw/forward_search_error/sbs8_loglh.pdf",
    size = 8
  )
