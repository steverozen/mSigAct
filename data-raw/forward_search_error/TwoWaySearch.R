save(spect, sigs, m.opts, max.mc.cores, p.thresh, optimal.exposure,
  file = "data-raw/forward_search_error/two_way_search_testdata.Rdata"
)

load("data-raw/forward_search_error/two_way_search_testdata.Rdata")

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

    for (step in seq_len(num_sig - 1)) {
      if (m.opts$trace >= 0) {
        message(
          "\nBackward search step ", step, ": ",
          ncol(sigs.to.test), " signatures to test"
        )
      }

      if (m.opts$trace >= 10) {
        message(
          "Signatures to test at this step: ",
          paste(colnames(sigs.to.test), collapse = " ")
        )
      }

      retval <-
        parallel::mclapply(1:ncol(sigs.to.test), FUN = function(i) {
          optimal.sigs <- sigs.to.test[, -i, drop = FALSE]

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
              df = step, lower.tail = FALSE
            )
          }

          return(list(p.value = chisq.p, exposure = try.exp$exposure))
        }, mc.cores = max.mc.cores)
      names(retval) <- colnames(sigs.to.test)

      p.values <- sapply(retval, FUN = "[[", "p.value")

      if (m.opts$trace >= 10) {
        message("p-values for each signature tested: ")
        print(p.values)
        message("p.thresh: ", p.thresh)
      }

      # Find out the largest p-value
      index <- which.max(p.values)
      largest.p.value <- p.values[index]

      if (largest.p.value < p.thresh) {
        if (m.opts$trace >= 0) {
          message("\nCannot remove any signature at this step")
          message("\nFinished backward search for ", colnames(spect))
        }

        return(optimal.exposure = optimal.exposure)
      }

      removed.sig <- sigs.to.test[, index, drop = FALSE]
      optimal.sigs <- sigs.to.test[, -index, drop = FALSE]

      if (m.opts$trace >= 1) {
        message("Signature removed at this step: ", colnames(removed.sig))
        message(
          "Total signatures remained: ",
          paste(colnames(optimal.sigs), collapse = " ")
        )
      }

      sigs.to.test <- optimal.sigs
      end <-
        OptimizeExposure(spectrum = spect, sigs = optimal.sigs, m.opts = m.opts)
      optimal.exposure <- end$exposure

      if (step == num_sig - 1) {
        if (m.opts$trace >= 0) {
          message("\nFinished backward search for ", colnames(spect))
        }
        return(optimal.exposure = optimal.exposure)
      }
    }
  }


TwoWaySearch <- function(spect, sigs, m.opts, max.mc.cores, p.thresh) {
  optimal.exposure <-
    ForwardSearch(
      spect = spect, sigs = sigs,
      m.opts = m.opts, max.mc.cores = max.mc.cores,
      p.thresh = p.thresh
    )
  
  optimal.exposure <-
    BackwardSearch(
      spect = spect, sigs = sigs,
      m.opts = m.opts, max.mc.cores = max.mc.cores,
      p.thresh = p.thresh, optimal.exposure = optimal.exposure
    )
  
  return(optimal.exposure = optimal.exposure)
    
}

retval <-
  TwoWaySearch(
    spect = spect, sigs = sigs,
    m.opts = m.opts, max.mc.cores = max.mc.cores,
    p.thresh = p.thresh
  )
