#' Find Maximum A Posteriori (MAP) assignment of signature exposures that
#' explain multiple spectra
#'
#' @inheritParams MAPAssignActivityInternal
#'
#' @param spectra The spectra (multiple spectra) to be reconstructed.
#'
#' @param output.dir Directory path to save the output file.
#'
#' @param num.parallel.samples The (maximum) number of samples to run in
#'   parallel. On Microsoft Windows machines it is silently changed to 1. Each
#'   sample in turn can require multiple cores, as governed by
#'   \code{mc.cores.per.sample}.
#'
#' @param mc.cores.per.sample The maximum number of cores to use for each
#'   sample. On Microsoft Windows machines it is silently changed to 1.
#'
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: Proposed signature assignment for \code{spectra}
#' with the highest MAP found.
#'    
#' * \code{proposed.reconstruction}: Proposed reconstruction of \code{spectra} based on \code{MAP}.
#' 
#' * \code{reconstruction.distances}: Various distances and similarities
#' between \code{spectra} and \code{proposed.reconstruction}.
#' 
#' * \code{error.messages}: Only appearing if there are errors running
#' \code{MAPAssignActivity}.
#'
#' * \code{results.details}: Detailed results for each sample in \code{spectra}.
#' 
#' These elements will be \code{NULL} if the algorithm could not find the
#' optimal reconstruction or there are errors coming out.
#' 
#' @md
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' # This is a long running example unless multiple CPU cores are available
#' indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
#' spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
#' sigs <- PCAWG7::COSMIC.v3.1$signature$genome$SBS96
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96", 
#'                                  cancer.type = "Lung-AdenoCA")
#' MAP.out <- MAPAssignActivity(spectra = spectra, 
#'                              sigs = sigs, 
#'                              sigs.presence.prop = sigs.prop, 
#'                              output.dir = file.path(tempdir(), "Lung-AdenoCA"),
#'                              max.level = length(sigs.prop) - 1,
#'                              p.thresh = 0.05 / ncol(spectra),
#'                              num.parallel.samples = 2,
#'                              mc.cores.per.sample = 10)
#'}
MAPAssignActivity <-
  function(spectra,
           sigs,
           sigs.presence.prop,
           output.dir,
           max.level               = 5,
           p.thresh                = 0.05,
           m.opts                  = DefaultManyOpts(),
           num.parallel.samples    = 5,
           mc.cores.per.sample     = min(20, 2^max.level),
           progress.monitor        = NULL,
           seed                    = NULL) {
    f1 <- function(i) {
      retval1 <- RunMAPOnOneSample(
        spect                   = spectra[ , i, drop = FALSE],
        sigs                    = sigs,
        sigs.presence.prop      = sigs.presence.prop,
        output.dir              = output.dir,
        max.level               = max.level,
        p.thresh                = p.thresh,
        m.opts                  = m.opts,
        max.mc.cores            = mc.cores.per.sample,
        progress.monitor        = progress.monitor,
        seed                    = seed)
      
      return(retval1)
    }
    
    num.parallel.samples <- Adj.mc.cores(num.parallel.samples)
    
    retval <- parallel::mclapply(1:ncol(spectra),
                                 f1,
                                 mc.cores = num.parallel.samples)
    names(retval) <- colnames(spectra)
    check.mclapply.result(
      retval, "MAPAssignActivity", colnames(spectra))
    
    proposed.assignment <- GetExposureInfo(list.of.MAP.out = retval)
    # Replace NA to 0 in proposed.assignment
    proposed.assignment[is.na(proposed.assignment)] <- 0
    
    proposed.reconstruction <- GetReconstructionInfo(list.of.MAP.out = retval)
    # Add attributes to proposed.reconstruction to be same as spectra
    proposed.reconstruction <- ICAMS::as.catalog(proposed.reconstruction)
    attr(proposed.reconstruction, "ref.genome") <- attr(spectra, "ref.genome")
    attr(proposed.reconstruction, "region") <- attr(spectra, "region")
    attr(proposed.reconstruction, "abundance") <- attr(spectra, "abundance")
    
    reconstruction.distances <- GetDistanceInfo(list.of.MAP.out = retval)
    error.messages <- lapply(retval, FUN = function(x) {
      return(x$error.messages)
    })
    # Remove NULL elements from error.messages
    error.messages <- Filter(Negate(is.null), error.messages)
    
    if (length(error.messages) == 0) {
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances,
                  results.details              = retval))
    } else {
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances,
                  error.messages               = error.messages,
                  results.details              = retval))
    }
  }

#' Retrieve exposure information from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#' 
#' @keywords internal
GetExposureInfo <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      exposures <- x$proposed.assignment
      count <- matrix(exposures$count, ncol = nrow(exposures))
      colnames(count) <- exposures$sig.id
      return(as.data.frame(count))
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  tmp1 <- tmp[index.of.non.null]
  retval <- do.call(dplyr::bind_rows, tmp1)
  retval1 <- t(retval)
  colnames(retval1) <- names(list.of.MAP.out)[index.of.non.null]
  
  retval2 <- retval1[SortSigId(rownames(retval1)), ]
  return(retval2)
}

#' Retrieve reconstruction information from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#' 
#' @keywords internal
GetReconstructionInfo <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      reconstruction <- x$proposed.reconstruction
      return(reconstruction)
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  tmp1 <- tmp[index.of.non.null]
  retval <- do.call(cbind, tmp1)
  colnames(retval) <- names(list.of.MAP.out)[index.of.non.null]
  
  return(retval)
}

#' Retrieve distance information from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#' 
#' @keywords internal
GetDistanceInfo <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      distances <- x$reconstruction.distances
      values <- matrix(distances$value, ncol = nrow(distances))
      colnames(values) <- distances$method
      return(as.data.frame(values))
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  
  tmp1 <- tmp[index.of.non.null]
  retval <- do.call(dplyr::bind_rows, tmp1)
  rownames(retval) <- names(list.of.MAP.out)[index.of.non.null]
  return(retval)
}

#' Run \code{MAPAssignActivity1} on one sample and save and plot the results
#' 
#' @inheritParams MAPAssignActivityInternal
#' 
#' @param output.dir Directory path to save the output file.
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
           seed                    = NULL) {
    
    if (!dir.exists(output.dir)) {
      dir.create(output.dir)
    }
    
    spect.name <- colnames(spect)
    retval <- mSigAct::MAPAssignActivity1(
      spect                   = spect,
      sigs                    = sigs,
      sigs.presence.prop      = sigs.presence.prop,
      max.level               = max.level,
      p.thresh                = p.thresh,
      m.opts                  = m.opts,
      max.mc.cores            = max.mc.cores,
      progress.monitor        = progress.monitor,
      seed                    = seed)
    
    if (!is.null(retval$error.messages)) {
      return(retval)
    }
    
    # Get the mutation type of the spectrum
    mut.type <- GetMutationType(spect)
    
    # We cannot use "::" in the file path, otherwise error will occur on Windows
    spect.name <- gsub(pattern = "::", replacement = ".", spect.name)
    
    output.path <- file.path(output.dir, spect.name)
    dir.create(path = output.path, showWarnings = FALSE)
    save(retval, file = file.path(output.path, 
                                  paste0(spect.name, ".", mut.type, ".MAP.Rdata")))
    
    distance.info <- retval$reconstruction.distances
    
    write.csv(distance.info, 
              file = file.path(output.path, 
                               paste0(spect.name, ".", mut.type, ".distances.csv")), 
              row.names = TRUE)
    
    ICAMS::WriteCatalog(catalog = ICAMS::as.catalog(spect), 
                        file = file.path(output.path, 
                                         paste0(spect.name, ".", mut.type, 
                                                ".catalog.csv")))
    
    tmp <- dplyr::arrange(retval$proposed.assignment, dplyr::desc(.data$count))
    sig.names <- tmp$sig.id
    inferred.exposure <- matrix(tmp$count)
    rownames(inferred.exposure) <- sig.names
    colnames(inferred.exposure) <- colnames(spect)
    
    ICAMSxtra::WriteExposure(exposure = inferred.exposure, 
                             file = file.path(output.path, 
                                              paste0(spect.name, ".", mut.type, 
                                                     ".inferred.exposure.csv")))
    ICAMSxtra::PlotExposureToPdf(inferred.exposure, 
                                 file = file.path(output.path, 
                                                  paste0(spect.name , ".", mut.type, 
                                                         ".inferred.exposure.pdf")))
    
    sigs.names <- rownames(inferred.exposure)
    sigs1 <- sigs[, sigs.names, drop = FALSE]
    
    
    if (!is.null(mut.type)) {
      etiologies <- sigs.etiologies[[mut.type]]
      colnames(sigs1) <- 
        paste0(colnames(sigs1), " (exposure = ", round(inferred.exposure[, 1]),
               ", contribution = ", 
               round(inferred.exposure[, 1]/sum(inferred.exposure[, 1]), 2), ") ",
               etiologies[colnames(sigs1), ])
    } else {
      colnames(sigs1) <- 
        paste0(colnames(sigs1), " (exposure = ", round(inferred.exposure[, 1]),
               ", contribution = ", 
               round(inferred.exposure[, 1]/sum(inferred.exposure[, 1]), 2), ")")
    }
    
    reconstructed.spectrum <- retval$proposed.reconstruction
    colnames(reconstructed.spectrum) <- 
      paste0("Reconstructed spectrum (count = ", round(colSums(reconstructed.spectrum)),
             ", cosine similarity = ", round(distance.info$value["cosine"], 5), ")")
    colnames(spect) <- paste0(colnames(spect), " (count = ",colSums(spect), ")")
    list.of.catalogs <- list(spect, reconstructed.spectrum, sigs1)
    PlotListOfCatalogsToPdf(list.of.catalogs,
                            file = file.path(output.path, 
                                             paste0(spect.name, ".", mut.type, 
                                                    ".proposed.reconstruction.pdf")))
    return(retval)
  }

#' @keywords internal
SortSigId <- function(sig.id) {
  num <- ICAMSxtra::NumFromId(sig.id)
  sig.id2 <- sig.id[order(num)]
  return(sig.id2)
}

#' @keywords internal
GetMutationType <- function(spect) {
  if (nrow(spect) == 96) {
    return("SBS96")
  } else if (nrow(spect) == 192) {
    return("SBS192")
  } else if (nrow(spect) == 78) {
    return("DBS78")
  } else if (nrow(spect) == 83) {
    return("ID")
  } else {
    return(NULL)
  }
}

<<<<<<< HEAD
=======
#' Retrieve distance information from the output generated by \code{MAPAssignActivity}
#' 
#' @keywords internal
GetDistanceInfo <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (isTRUE(x$success)) {
      distances <- x$MAP.distances
      values <- matrix(distances$value, ncol = nrow(distances))
      colnames(values) <- distances$method
      return(as.data.frame(values))
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  
  tmp1 <- tmp[index.of.non.null]
  retval <- do.call(dplyr::bind_rows, tmp1)
  rownames(retval) <- names(list.of.MAP.out)[index.of.non.null]
  return(retval)
}

#' Retrieve exposure information from the output generated by \code{MAPAssignActivity}
#' 
#' @param list.of.MAP.out Output generated by \code{MAPAssignActivity}.
#'
#' @param type.of.exposure Type of exposure information to retrieve, either
#'   "MAP" or "best.sparse".
#'
#' @keywords internal
GetExposureInfo <- function(list.of.MAP.out, type.of.exposure = "MAP") {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (isTRUE(x$success)) {
    exposures <- x[[type.of.exposure]]
    count <- matrix(exposures$count, ncol = nrow(exposures))
    colnames(count) <- exposures$sig.id
    return(as.data.frame(count))
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  tmp1 <- tmp[index.of.non.null]
  retval <- do.call(dplyr::bind_rows, tmp1)
  retval1 <- t(retval)
  colnames(retval1) <- names(list.of.MAP.out)[index.of.non.null]
  
  retval2 <- retval1[SortSigId(rownames(retval1)), ]
  
  # Change NA to 0
  retval2[is.na(retval2)] <- 0
  return(retval2)
}

>>>>>>> c54bc7811efec976bea963ae3fabecdd8d80b701
#' Plot List of catalogs to Pdf
#' 
#' @param list.of.catalogs List of catalogs in \code{\link{ICAMS}} format.
#'
#' @inheritParams ICAMS::PlotCatalogToPdf
#' 
#' @keywords internal
PlotListOfCatalogsToPdf <- function(list.of.catalogs, 
                                    file, 
                                    plot.SBS12 = FALSE, 
                                    cex     = 0.8,
                                    grid    = TRUE, 
                                    upper   = TRUE, 
                                    xlabels = TRUE,
                                    ylim    = NULL) {
  old.par.tck.value <- graphics::par("tck")
  # Setting the width and length for A4 size plotting
  grDevices::pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  graphics::par(tck = old.par.tck.value)
  
  num.of.catalogs <- length(list.of.catalogs)
  if (nrow(list.of.catalogs[[1]]) == 96) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  } else if (nrow(list.of.catalogs[[1]]) == 192) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 2, 1, 1))
  } else if (nrow(list.of.catalogs[[1]]) == 78) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 3, 2, 2))
  } else if (nrow(list.of.catalogs[[1]]) == 83) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))
  } 
  
  on.exit(graphics::par(opar))
  
  for (i in 1:num.of.catalogs) {
    catalog <- list.of.catalogs[[i]]
    num.of.samples <- ncol(catalog)
    
    for (j in 1:num.of.samples) {
      cat <- catalog[, j, drop = FALSE]
      ICAMS::PlotCatalog(cat, plot.SBS12 = plot.SBS12, cex = cex, grid = grid, 
                         upper = upper, xlabels = xlabels, ylim = ylim)
    }
    
  }
  
  grDevices::dev.off()
  invisible(list(plot.success = TRUE))
}
