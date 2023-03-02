ForwardSearch <- function(spect, sigs, m.opts, max.mc.cores, p.thresh) {
  start <- OptimizeExposure(spectrum = spect, sigs = sigs, m.opts = m.opts)
  lh.w.all <- start$loglh # The likelihood with all signatures
  num_sig <- ncol(sigs)
  
  if (num_sig == 1) {
    return(optimal.exposure = start$exposure)
  }

  if (m.opts$trace >= 0) {
    message("\nStarting forward search for ", colnames(spect))
  }
  
  if (m.opts$trace >= 1) {
    message("log likelihood using all signatures = ", lh.w.all)
  }

  sigs.to.choose <- sigs
  optimal.sigs <- sigs[, 0, drop = FALSE]

  for (step in seq_len(num_sig - 1)) {
    if (m.opts$trace >= 0) {
      message(
        "\nForward search step ", step, ": ",
        ncol(sigs.to.choose), " signatures to test"
      )
    }

    if (m.opts$trace >= 10) {
      message(
        "Signatures to test at this step: ",
        paste(colnames(sigs.to.choose), collapse = " ")
      )
    }

    retval <-
      parallel::mclapply(1:ncol(sigs.to.choose), FUN = function(i) {
        one.sig <- sigs.to.choose[, i, drop = FALSE]
        optimal.sigs <- cbind(optimal.sigs, one.sig)

        try.exp <- OptimizeExposure(
          spectrum = spect,
          sigs = optimal.sigs, m.opts = m.opts
        )
        if (is.infinite(try.exp$loglh)) {
          # It was not possible to generate the spectrum from the signatures, so
          # try.exp$loglh should be -Inf and the p value the test that the spectrum
          # can be better reconstructed with optimal.sigs than with sigs is 0.
          if (try.exp$loglh < 0) {
            chisq.p <- 0
          } else {
            stop("Probable coding error")
          }
        } else {
          statistic <- 2 * (lh.w.all - try.exp$loglh)
          chisq.p <- stats::pchisq(
            q = statistic,
            df = num_sig - step, lower.tail = FALSE
          )
        }

        return(list(p.value = chisq.p, exposure = try.exp$exposure))
      }, mc.cores = max.mc.cores)
    names(retval) <- colnames(sigs.to.choose)

    p.values <- sapply(retval, FUN = "[[", "p.value")

    if (m.opts$trace >= 10) {
      message("p-values for each signature tested: ")
      print(p.values)
      message("p.thresh: ", p.thresh)
    }

    # Find out the largest p-value
    index <- which.max(p.values)
    largest.p.value <- p.values[index]
    selected.sig <- sigs.to.choose[, index, drop = FALSE]
    optimal.sigs <- cbind(optimal.sigs, selected.sig)

    if (m.opts$trace >= 1) {
      message("Signature selected at this step: ", colnames(selected.sig))
      message(
        "Total signatures selected: ",
        paste(colnames(optimal.sigs), collapse = " ")
      )
    }

    if (largest.p.value > p.thresh) {
      if (m.opts$trace >= 0) {
        message("\nFinished forward search for ", colnames(spect))
      }

      return(optimal.exposure = retval[[index]][["exposure"]])
    }

    sigs.to.choose <- sigs.to.choose[, -index, drop = FALSE]

    # At the last step, if the largest.p.value is still <= p.thresh,
    # then we will use all signatures
    if (step == num_sig - 1) {
      if (m.opts$trace >= 0) {
        message("\nNeed to use all signatures to reconstruct the spectrum: ", 
                paste(colnames(sigs), collapse = " "))
      }
      return(optimal.exposure = start$exposure)
    }
  }
}
