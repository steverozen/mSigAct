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
#' The elements \code{proposed.assignment}, \code{proposed.reconstruction},
#' \code{reconstruction.distances} will be \code{NULL} if the algorithm could
#' not find the optimal reconstruction or there are errors coming out for
#' \strong{all} samples.
#' 
#' @md
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This is a long running example unless parallel computing is supported on your machine
#' indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
#' spectra <- PCAWG7::spectra$PCAWG$SBS96[, indices[1:2], drop = FALSE]
#' sigs <- PCAWG7::signature$genome$SBS96
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
           seed                    = NULL,
           max.subsets             = 1000) {
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
        seed                    = seed,
        max.subsets             = max.subsets)

      return(retval1)
    }

    num.parallel.samples <- Adj.mc.cores(num.parallel.samples)

    retval <- parallel::mclapply(1:ncol(spectra),
                                 f1,
                                 mc.cores = num.parallel.samples)
    
    names(retval) <- colnames(spectra)
    
    error.messages <- lapply(retval, FUN = function(x) {
      return(x$error.messages)
    })
    # Remove NULL elements from error.messages
    error.messages <- Filter(Negate(is.null), error.messages)
    
    # Check for samples which have NULL return for proposed.assignment
    null.assignment <- sapply(retval, FUN = function(x) {
      return(is.null(x$proposed.assignment))
    })
    retval.non.null <- retval[!null.assignment]
    
    if (length(retval.non.null) == 0) {
      # The case when all samples have NULL assignment
      return(list(proposed.assignment          = NULL,
                  proposed.reconstruction      = NULL,
                  reconstruction.distances     = NULL,
                  error.messages               = error.messages,
                  results.details              = retval))
    }
    
    proposed.assignment <- GetExposureInfo(list.of.MAP.out = retval.non.null)
    # Replace NA to 0 in proposed.assignment
    proposed.assignment[is.na(proposed.assignment)] <- 0
    
    proposed.reconstruction <- GetReconstructionInfo(list.of.MAP.out = retval.non.null)
    # Add attributes to proposed.reconstruction to be same as spectra
    proposed.reconstruction <- AddAttributes(proposed.reconstruction, spectra)

    reconstruction.distances <- GetDistanceInfo(list.of.MAP.out = retval.non.null)


    if (length(error.messages) == 0) {
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances,
                  results.details              = retval))
    } else {
      # The case when part of samples have NULL assignment
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances,
                  error.messages               = error.messages,
                  results.details              = retval))
    } 
  }

#' Add attributes to proposed.reconstruction to be same as spectra
#'
#' @keywords internal
AddAttributes <- function(proposed.reconstruction, spectra) {
  proposed.reconstruction <- ICAMS::as.catalog(proposed.reconstruction)
  attr(proposed.reconstruction, "ref.genome") <- attr(spectra, "ref.genome")
  attr(proposed.reconstruction, "region") <- attr(spectra, "region")
  attr(proposed.reconstruction, "abundance") <- attr(spectra, "abundance")
  return(proposed.reconstruction)
}

#' Retrieve exposure information from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#'
#' @keywords internal
GetExposureInfo <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      exposures <- x$proposed.assignment
      count <- matrix(exposures[, 1], ncol = nrow(exposures))
      colnames(count) <- rownames(exposures)
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

  retval2 <- retval1[SortSigId(rownames(retval1)), , drop = FALSE]
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
  
  if (length(tmp1) == 1) {
    return(tmp1[[1]])
  }
  
  retval <- do.call(cbind, tmp1)
  colnames(retval) <- names(list.of.MAP.out)[index.of.non.null]

  return(retval)
}

#' Retrieve distance information from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#'
#' @keywords internal
GetDistanceInfo <- function(list.of.MAP.out) {
  MAP.tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      distances <- x$reconstruction.distances
      values <- matrix(distances$proposed.assignment, ncol = nrow(distances))
      colnames(values) <- distances$method
      return(as.data.frame(values))
    } else {
      return(NULL)
    }
  })
  MAP.index.of.non.null <- sapply(MAP.tmp, FUN = Negate(is.null))

  MAP.tmp1 <- MAP.tmp[MAP.index.of.non.null]
  MAP.retval <- do.call(dplyr::bind_rows, MAP.tmp1)
  rownames(MAP.retval) <- names(list.of.MAP.out)[MAP.index.of.non.null]
  
  QP.tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      distances <- x$reconstruction.distances
      values <- matrix(distances$QP.assignment, ncol = nrow(distances))
      colnames(values) <- distances$method
      return(as.data.frame(values))
    } else {
      return(NULL)
    }
  })
  QP.index.of.non.null <- sapply(QP.tmp, FUN = Negate(is.null))
  
  QP.tmp1 <- QP.tmp[QP.index.of.non.null]
  QP.retval <- do.call(dplyr::bind_rows, QP.tmp1)
  rownames(QP.retval) <- names(list.of.MAP.out)[QP.index.of.non.null]
  
  return(list(MAP.distances = MAP.retval, QP.distances = QP.retval))
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
