#' Find Maximum A Posteriori (MAP) assignment of signature exposures that
#' explain multiple spectra
#' 
#' This function also can do sparse assignment by specifying \code{use.sparse.assign = TRUE}.
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
#' @param drop.low.mut.samples Whether to exclude low mutation samples from
#' the analysis. If \code{TRUE(default)}, samples with SBS total mutations less
#' than 100, DBS or ID total mutations less than 25 will be dropped.
#'
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: Proposed signature assignment for \code{spectra}
#' with the highest MAP found. If \code{use.sparse.assign = TRUE}, this will
#' be the most sparse set of signatures that can plausibly explain \code{spectra}.
#'
#' * \code{proposed.reconstruction}: Proposed reconstruction of \code{spectra} based on \code{MAP}.
#' If \code{use.sparse.assign = TRUE}, this will be the reconstruction based on
#' sparse assignment.
#'
#' * \code{reconstruction.distances}: Various distances and similarities
#' between \code{spectra} and \code{proposed.reconstruction}.
#' 
#' * \code{all.tested}: All tested possible ways to reconstruct each
#' sample in \code{spectra}.
#' 
#' * \code{alt.solutions}: A \code{tibble} showing all the alternative solutions
#' that are statistically as good as the \code{proposed.assignment} that can
#' plausibly reconstruct \code{spectra}.
#' 
#' * \code{time.for.assignment}: Value from \code{system.time} for running
#'  \code{MAPAssignActivity1} for each sample in \code{spectra}.
#' 
#' * \code{error.messages}: Only appearing if there are errors running
#' \code{MAPAssignActivity}.
#'
#' The elements \code{proposed.assignment}, \code{proposed.reconstruction},
#' \code{reconstruction.distances}, \code{all.tested},
#' \code{time.for.assignment} will be \code{NULL} if the algorithm could not
#' find the optimal reconstruction or there are errors coming out for
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
           max.subsets             = 1000,
           use.sparse.assign       = FALSE,
           drop.low.mut.samples    = TRUE) {
    if (drop.low.mut.samples) {
      spectra <- DropLowMutationSamples(spectra)
    } else {
      spectra <- spectra
    }
    
    if (ncol(spectra) == 0) {
      return(NullReturnForMAPAssignActivity1(msg = "No sample to analyse"))
    }
    
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
        max.subsets             = max.subsets,
        use.sparse.assign       = use.sparse.assign)

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
      return(NullReturnForMAPAssignActivity1(msg = error.messages))
    }
    
    proposed.assignment <- GetExposureInfo(list.of.MAP.out = retval.non.null)
    # Replace NA to 0 in proposed.assignment
    proposed.assignment[is.na(proposed.assignment)] <- 0
    
    proposed.reconstruction <- GetReconstructionInfo(list.of.MAP.out = retval.non.null)
    # Add attributes to proposed.reconstruction to be same as spectra
    proposed.reconstruction <- AddAttributes(proposed.reconstruction, spectra)
    
    if (use.sparse.assign == FALSE) {
      reconstruction.distances <- GetDistanceInfo(list.of.MAP.out = retval.non.null)
    } else if (use.sparse.assign == TRUE) {
      reconstruction.distances <- GetDistanceInfo(list.of.MAP.out = retval.non.null,
                                                  sparse.assign = TRUE)
    }
    
    all.tested <- GetAllTestedTables(list.of.MAP.out = retval.non.null)
    alt.solutions <- GetAltSolutionsTables(list.of.MAP.out = retval.non.null)
    time.for.MAP.assign <- GetTimeForMAPAssign(list.of.MAP.out = retval.non.null)

    if (length(error.messages) == 0) {
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances,
                  all.tested                   = all.tested,
                  alt.solutions                = alt.solutions,
                  time.for.assignment          = time.for.MAP.assign))
    } else {
      # The case when part of samples have NULL assignment
      return(list(proposed.assignment          = proposed.assignment,
                  proposed.reconstruction      = proposed.reconstruction,
                  reconstruction.distances     = reconstruction.distances,
                  all.tested                   = all.tested,
                  alt.solutions                = alt.solutions,
                  time.for.assignment          = time.for.MAP.assign,
                  error.messages               = error.messages))
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
#' \code{MAPAssignActivity1} or \code{SparseAssignActivity1} on multiple samples
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
#' \code{MAPAssignActivity1} or \code{SparseAssignActivity1} on multiple samples
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
#' \code{MAPAssignActivity1} or \code{SparseAssignActivity1} on multiple samples
#'
#' @keywords internal
GetDistanceInfo <- function(list.of.MAP.out, sparse.assign = FALSE) {
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
  
  if (sparse.assign == FALSE) {
    return(list(MAP.distances = MAP.retval, QP.distances = QP.retval))
  } else {
    return(list(sparse.assign.distances = MAP.retval, QP.distances = QP.retval))
  }
  
}

#' Retrieve all tested tables from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#'
#' @keywords internal
GetAllTestedTables <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      all.tested <- x$all.tested
      return(all.tested)
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  tmp1 <- tmp[index.of.non.null]
  
  return(tmp1)
}

#' Retrieve alternative solutions tables from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#'
#' @keywords internal
GetAltSolutionsTables <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      all.tested <- x$alt.solutions
      return(all.tested)
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  tmp1 <- tmp[index.of.non.null]
  
  return(tmp1)
}

#' Retrieve time.for.MAP.assign from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#'
#' @keywords internal
GetTimeForMAPAssign <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    if (is.null(x$error.messages)) {
      xx <- x$time.for.MAP.assign
      df <- data.frame(user = xx[1] + xx[4], # user.self + user.child
                       system = xx[2] + xx[5],
                       elapsed = xx[3])
      return(df)
    } else {
      return(NULL)
    }
  })
  index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  tmp1 <- tmp[index.of.non.null]
  
  tmp2 <- do.call("rbind", tmp1)
  return(tmp2)
}

#' @keywords internal
SortSigId <- function(sig.id) {
  num <- ICAMSxtra::NumFromId(sig.id)
  sig.id2 <- sig.id[order(num)]
  return(sig.id2)
}

#' @keywords internal
GetMutationType <- function(spect) {
  spect <- as.matrix(spect)
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

#' @keywords internal
DropLowMutationSamples <- function(spectra) {
  spectra <- as.matrix(spectra)
  mut.type <- GetMutationType(spectra)
  if (mut.type %in% c("SBS96", "SBS192")) {
    thresh.value <- 100
  } else if (mut.type %in% c("DBS78", "ID")) {
    thresh.value <- 25
  }
  
  indices <- which(colSums(spectra) < thresh.value)
  
  if (length(indices) == 0) {
    return(spectra)
  } else {
    message("Samples with total mutations less than ", thresh.value,
            " were excluded in the analysis for ", mut.type, "\n",
            paste(colnames(spectra)[indices], collapse = " "),
            "\nSet argument drop.low.mut.samples = FALSE to suppress")
    return(spectra[, -indices, drop = FALSE])
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
