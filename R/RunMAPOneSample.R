#' Run \code{MAPAssignActivity1} on one sample and save and plot the results
#'
#' @inheritParams MAPAssignActivityInternal
#'
#' @param output.dir Directory path to save the output file.
#' 
#' @param drop.low.mut.samples Whether to exclude low mutation samples from
#' the analysis. If \code{TRUE(default)}, samples with SBS total mutations less
#' than 100, DBS or ID total mutations less than 25 will be dropped.
#'
#' @importFrom utils write.csv
#'
#' @keywords internal
RunMAPOnOneSample <-
  function(spect,
           sigs,
           sigs.presence.prop,
           output.dir,
           max.level               = 5,
           p.thresh                = 0.05,
           m.opts                  = DefaultManyOpts(),
           max.mc.cores            = min(20, 2^max.level),
           progress.monitor        = NULL,
           seed                    = NULL,
           max.subsets             = 1000,
           use.sparse.assign       = FALSE,
           drop.low.mut.samples    = TRUE) {
    
    if (drop.low.mut.samples) {
      spect <- DropLowMutationSamples(spect)
    } else {
      spect <- spect
    }
    
    if (!dir.exists(output.dir)) {
      dir.create(output.dir, recursive = TRUE)
    }
    
    spect.name <- colnames(spect)
    retval <- MAPAssignActivity1(
      spect                   = spect,
      sigs                    = sigs,
      sigs.presence.prop      = sigs.presence.prop,
      max.level               = max.level,
      p.thresh                = p.thresh,
      m.opts                  = m.opts,
      max.mc.cores            = max.mc.cores,
      progress.monitor        = progress.monitor,
      seed                    = seed,
      max.subsets             = max.subsets,
      use.sparse.assign       = use.sparse.assign,
      drop.low.mut.samples    = drop.low.mut.samples)
    
    if (!is.null(retval$error.messages)) {
      return(retval)
    }
    
    # Get the mutation type of the spectrum
    mut.type <- GetMutationType(spect)
    
    # We cannot use "::" in the file path, otherwise error will occur on Windows
    spect.name <- gsub(pattern = "::", replacement = ".", spect.name)
    
    output.path <- file.path(output.dir, spect.name)
    dir.create(path = output.path, showWarnings = FALSE)
    
    if (use.sparse.assign) {
      algorithm <- "sparse"
    } else {
      algorithm <- "MAP"
    }
    
    save(retval, file = file.path(output.path,
                                  paste0(spect.name, ".", mut.type, ".", algorithm, ".Rdata")))
    
    distance.info <- retval$reconstruction.distances
    
    write.csv(distance.info,
              file = file.path(output.path,
                               paste0(spect.name, ".", mut.type, ".distances.csv")),
              row.names = TRUE)
    
    ICAMS::WriteCatalog(catalog = ICAMS::as.catalog(spect),
                        file = file.path(output.path,
                                         paste0(spect.name, ".", mut.type,
                                                ".catalog.csv")))
    inferred.exposure <- retval$proposed.assignment
    
    # Order inferred.exposure by mutation counts
    inferred.exposure <-
      inferred.exposure[order(inferred.exposure[, 1], decreasing = TRUE), ,
                        drop = FALSE]
    colnames(inferred.exposure) <- colnames(spect)
    
    WriteExposure(exposure = inferred.exposure,
                             file = file.path(output.path,
                                              paste0(spect.name, ".", mut.type,
                                                     ".inferred.exposure.csv")))
    PlotExposureToPdf(inferred.exposure,
                                 file = file.path(output.path,
                                                  paste0(spect.name , ".", mut.type,
                                                         ".inferred.exposure.pdf")))
    
    sigs.names <- rownames(inferred.exposure)
    sigs1 <- sigs[, sigs.names, drop = FALSE]
    
    
    if (!is.null(mut.type)) {
      ets <- cosmicsig::get_etiology(mut.type, colnames(sigs1))
      colnames(sigs1) <-
        paste0(colnames(sigs1), " (exposure = ", round(inferred.exposure[, 1]),
               ", contribution = ",
               round(inferred.exposure[, 1]/sum(inferred.exposure[, 1]), 2), ") ",
               ets)
    } else {
      colnames(sigs1) <-
        paste0(colnames(sigs1), " (exposure = ", round(inferred.exposure[, 1]),
               ", contribution = ",
               round(inferred.exposure[, 1]/sum(inferred.exposure[, 1]), 2), ")")
    }
    
    reconstructed.spectrum <- retval$proposed.reconstruction
    colnames(reconstructed.spectrum) <-
      paste0("Reconstructed spectrum (count = ", round(colSums(reconstructed.spectrum)),
             ", cosine similarity = ", 
             round(distance.info$proposed.assignment["cosine"], 5), ")")
    subtracted.spect <- spect - reconstructed.spectrum
    colnames(subtracted.spect) <- 
      paste0(colnames(spect), " (substracted spectrum)")
    
    colnames(spect) <- paste0(colnames(spect), " (count = ",colSums(spect), ")")
    
    list.of.catalogs <- list(spect, reconstructed.spectrum, subtracted.spect,
                             sigs1)
    PlotListOfCatalogsToPdf(list.of.catalogs,
                            file = file.path(output.path,
                                             paste0(spect.name, ".", mut.type,
                                                    ".proposed.reconstruction.pdf")))
    return(retval)
  }
