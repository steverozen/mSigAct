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
#' @return A list of lists containing output for each sample in \code{spectra}. 
#' Each sublist has the following elements \describe{
#'
#' \item{MAP}{A 2-column \code{tibble} with the attributions with the highest MAP found.
#'    Column 1 contains signature ids; column 2 contains the associated counts. }
#'
#' \item{MAP.row}{A 1-row \code{tibble} with various information on the selected exposure.}
#'
#' \item{best.sparse}{A 2-column \code{tibble} with the most-sparse attributions with
#'      the highest MAP, in the same format as element \code{MAP}.}
#'
#' \item{best.sparse.row}{{A 1-row \code{tibble} with various information on the
#'    most-sparse exposure with the best MAP.}}
#'
#' \item{all.tested}{A \code{tibble} of all the search results.}
#'
#' \item{messages}{Possibly empty character vector with messages.}
#'
#' \item{success}{\code{TRUE} is search was successful, \code{FALSE} otherwise.}
#'
#' \item{time.for.MAP.assign}{Value from \code{system.time} for
#'  \code{\link{MAPAssignActivityInternal}}}.
#'
#' \item{MAP.recon}{Reconstruction based on \code{MAP}}.
#'
#' \item{sparse.MAP.recon}{Reconstruction based on \code{best.sparse}}.
#'
#' \item{MAP.distances}{Various distances and similarities
#' between \code{spect} and \code{MAP.recon}}.
#'
#' \item{sparse.MAP.distances}{Various distances and similarities
#' between \code{spect} and \code{sparse.MAP.recon}}.
#'
#' }
#'
#' These elements will be \code{NULL} if \code{max.subsets} is exceeded.
#'
#' @export
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
           max.subsets             = 1000,
           max.presence.proportion = 0.99,
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
        max.subsets             = max.subsets,
        max.presence.proportion = max.presence.proportion,
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
#' @import ICMAS 
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
           max.subsets             = 1000,
           max.presence.proportion = 0.99,
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
      max.subsets             = max.subsets,
      max.presence.proportion = max.presence.proportion,
      progress.monitor        = progress.monitor,
      seed                    = seed)
    
    if (!isTRUE(retval$success)) {
      return(retval)
    }
    
    output.path <- file.path(output.dir, spect.name)
    dir.create(path = output.path, showWarnings = FALSE)
    save(retval, file = file.path(output.path, paste0(spect.name, ".MAP.Rdata")))
    
    distance.info <- retval$MAP.distances
    
    write.csv(distance.info, 
              file = file.path(output.path, paste0(spect.name, ".distances.csv")), 
              row.names = TRUE)
    
    ICAMS::WriteCatalog(catalog = ICAMS::as.catalog(spect), 
                        file = file.path(output.path, 
                                         paste0(spect.name , ".catalog.csv")))
    
    tmp <- dplyr::arrange(retval$MAP, dplyr::desc(.data$count))
    sig.names <- tmp$sig.id
    inferred.exposure <- matrix(tmp$count)
    rownames(inferred.exposure) <- sig.names
    colnames(inferred.exposure) <- colnames(spect)
    
    ICAMSxtra::WriteExposure(exposure = inferred.exposure, 
                             file = file.path(output.path, 
                                              paste0(spect.name , ".inferred.exposure.csv")))
    ICAMSxtra::PlotExposureToPdf(inferred.exposure, 
                                 file = file.path(output.path, 
                                                  paste0(spect.name , ".inferred.exposure.pdf")))
    
    sigs.names <- rownames(inferred.exposure)
    sigs1 <- sigs[, sigs.names, drop = FALSE]
    
    sig.type <- GetSigType(sigs)
    if (!is.null(sig.type)) {
      etiologies <- sigs.etiologies[[sig.type]]
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
    
    reconstructed.spectrum <- round(retval$MAP.recon)
    colnames(reconstructed.spectrum) <- 
      paste0("reconstructed (count = ", round(colSums(reconstructed.spectrum)),
             ", cosine similarity = ", round(distance.info$value[3], 5), ")")
    colnames(spect) <- paste0(colnames(spect), " (count = ",colSums(spect), ")")
    list.of.catalogs <- list(spect, reconstructed.spectrum, sigs1)
    PlotListOfCatalogsToPdf(list.of.catalogs,
                            file = file.path(output.path, 
                                             paste0(spect.name , 
                                                    ".reconstructed.spectrum.pdf")))
    return(retval)
  }

#' @keywords internal
SortSigId <- function(sig.id) {
  num <- ICAMSxtra::NumFromId(sig.id)
  sig.id2 <- sig.id[order(num)]
  return(sig.id2)
}

#' @keywords internal
GetSigType <- function(sigs) {
  if (nrow(sigs) == 96) {
    return("SBS96")
  } else if (nrow(sigs) == 192) {
    return("SBS192")
  } else if (nrow(sigs) == 78) {
    return("DBS78")
  } else if (nrow(sigs) == 83) {
    return("ID")
  } else {
    return(NULL)
  }
}

#' Plot List of catalogs to Pdf
#' 
#' @param list.of.catalogs List of catalogs in \code{\link{ICAMS}} format.
#'
#' @inheritParams ICAMS::PlotCatalogToPdf
#' 
#' @importFrom  ICAMS PlotCatalog
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
  catalog.type <- attr(list.of.catalogs[[1]], "catalog.type")
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
