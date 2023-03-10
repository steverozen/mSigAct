#' @keywords internal
CompareAndPlotLoglh <-
  function(spectrum, exposure, sigs, sig.to.test, neg.binom.size, file) {
    if (!is.matrix(exposure) && is.numeric(exposure) && !is.null(names(exposure))) {
      exposure <- as.matrix(exposure)
    } else {
      stop("exposure should be a matrix with rownames or a numeric vector with names")
    }

    my.opts <- DefaultManyOpts(likelihood.dist = "neg.binom")
    my.opts$nbinom.size <- neg.binom.size

    sigs <- sigs[, rownames(exposure), drop = FALSE]

    original.recon <- ReconstructSpectrum(
      sigs = sigs,
      exp = exposure,
      use.sig.names = TRUE
    )
    original.loglh <- stats::dnbinom(
      x = spectrum, mu = original.recon, size = neg.binom.size, log = TRUE
    )

    new.sigs <- sigs[, colnames(sigs) != sig.to.test, drop = FALSE]
    optimized.exposure <-
      OptimizeExposure(
        spectrum = spectrum,
        sigs = new.sigs,
        m.opts = my.opts
      )
    new.exposure <- optimized.exposure$exposure

    new.recon <- ReconstructSpectrum(
      sigs = sigs,
      exp = new.exposure,
      use.sig.names = TRUE
    )

    new.loglh <- stats::dnbinom(
      x = spectrum, mu = new.recon, size = neg.binom.size, log = TRUE
    )

    original.neg.loglh <- original.loglh * -1
    new.neg.loglh <- new.loglh * -1
    substracted.loglh <- new.neg.loglh - original.neg.loglh

    to.plot <- cbind(
      new.neg.loglh, original.neg.loglh,
      substracted.loglh
    )

    sample.name <- colnames(spectrum)
    colnames(to.plot) <-
      c(
        paste0(
          sample.name,
          " (negative loglh without ", sig.to.test, ": ",
          round(colSums(new.neg.loglh), 2), ")"
        ),
        paste0(
          sample.name, " (negative loglh with ", sig.to.test, ": ",
          round(colSums(original.neg.loglh), 2), ")"
        ),
        paste0(
          sample.name, " (subtracted negative loglh: ",
          round(colSums(substracted.loglh), 2), ")"
        )
      )

    test.sig.spectra <- sigs[, sig.to.test, drop = FALSE]
    list.of.catalogs <- list(to.plot, test.sig.spectra)

    PlotListOfCatalogsToPdf(
      list.of.catalogs = list.of.catalogs,
      file = file
    )
    invisible(TRUE)
  }
