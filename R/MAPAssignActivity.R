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
#' the analysis. If \code{TRUE}, samples with SBS total mutations less
#' than 100, DBS or ID total mutations less than 25 will be dropped.
#' 
#' @param sig.pres.test.q.thresh Test parameter.
#' 
#' @param save.files If \code{TRUE} save several files for each input sample
#'   in a directory named after the sample. 
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
#' SBS96.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96",
#'                                  cancer.type = "Lung-AdenoCA")
#' sigs.to.use <- SBS96.sigs[, names(sigs.prop), drop = FALSE]
#' MAP.out <- MAPAssignActivity(spectra = spectra,
#'                              sigs = sigs.to.use,
#'                              sigs.presence.prop = sigs.prop,
#'                              output.dir = file.path(tempdir(), "Lung-AdenoCA"),
#'                              max.level = ncol(sigs.to.use) - 1,
#'                              p.thresh = 0.05 / ncol(sigs.to.use),
#'                              num.parallel.samples = 2,
#'                              mc.cores.per.sample = 30,
#'                              seed = 2561)
#'}
MAPAssignActivity <-
  function(spectra,
           sigs,
           sigs.presence.prop,
           output.dir,
           max.level                  = 5,
           p.thresh                   = 0.05,
           m.opts                     = DefaultManyOpts(),
           num.parallel.samples       = 5,
           mc.cores.per.sample        = min(20, 2^max.level),
           progress.monitor           = NULL,
           seed                       = NULL,
           max.subsets                = 1000,
           use.sparse.assign          = FALSE,
           use.forward.search         = FALSE,
           drop.low.mut.samples       = TRUE,
           use.sig.presence.test      = FALSE,
           sig.pres.test.nbinom.size  = NULL,
           sig.pres.test.p.thresh     = 0.05,
           sig.pres.test.q.thresh     = NULL,
           save.files                 = TRUE) {
    
    null.assignment1 <- matrix(rep(0, ncol(sigs)))
    colnames(null.assignment1) <- "No samples"
    rownames(null.assignment1) <- colnames(sigs)
    null.spect1       <- matrix(rep(0, nrow(sigs)))
    colnames(null.spect1) <- "No samples"
    
    if (ncol(spectra) == 0) {
      return(NullReturnForMAPAssignActivity(signature.universe = sigs, 
                                            msg = "No samples to analyse"))
    }
    
    f1 <- function(i) {
      retval1 <- RunMAPOnOneSample(
        spect                       = spectra[ , i, drop = FALSE],
        sigs                        = sigs,
        sigs.presence.prop          = sigs.presence.prop,
        output.dir                  = output.dir,
        max.level                   = max.level,
        p.thresh                    = p.thresh,
        m.opts                      = m.opts,
        max.mc.cores                = mc.cores.per.sample,
        progress.monitor            = progress.monitor,
        seed                        = seed,
        max.subsets                 = max.subsets,
        use.sparse.assign           = use.sparse.assign,
        use.forward.search          = use.forward.search,
        drop.low.mut.samples        = drop.low.mut.samples, 
        use.sig.presence.test       = use.sig.presence.test,
        sig.pres.test.nbinom.size   = sig.pres.test.nbinom.size,
        sig.pres.test.p.thresh      = sig.pres.test.p.thresh,
        sig.pres.test.q.thresh      = sig.pres.test.q.thresh,
        save.files                  = save.files)

      return(retval1)
    }

    num.parallel.samples <- Adj.mc.cores(num.parallel.samples)

    retval <- parallel::mclapply(1:ncol(spectra),
                                 f1,
                                 mc.cores = num.parallel.samples)
    
    # retval is a list with each element a value returned from
    # RunMAPOnOneSample. Each element is a list with a proposed.assignment and a
    # proposed.solution, in addition to other slots.
    names(retval) <- colnames(spectra) 
    
    clean_retval <- function(nn) {

      null.retval <- NullReturnForMAPAssignActivity(
        signature.universe = sigs, 
        target.spectrum    = spectra[ , nn, drop = FALSE]) # Careful, it must be that names(retval) == colnames(spectra)
      
      rr <- retval[[nn]]
      
      if (is.null(rr)) {
        null.retval$error.messages <- "R in child process possibly crashed"
        return(null.retval)
      }
      if (is.null(rr$proposed.assignment)) {
        if (inherits(rr, "try-error")) {
          null.retval$error.messages <- rr[1]
        } else {
          null.retval$error.messages <- 
            paste("Unexpected return value from child, class is :",
                  class(rr))
        }
        return(null.retval)
      }
      return(rr)
    }

    retval1 <- lapply(names(retval), clean_retval)
    names(retval1) <- names(retval)
    rm(retval)
    
    error.messages <- unlist(lapply(retval1, '[[', "error.messages"))
    names(error.messages) <- names(retval1)

    proposed.assignment <- GetExposureInfo(retval1)
    proposed.assignment[is.na(proposed.assignment)] <- 0 # FIX TODO TO DO check w/ Nanhai if we need to remove 0-exposure signatures earlier
    
    # proposed.reconstruction <- GetReconstructionInfo(list.of.MAP.out = retval1)
    proposed.reconstruction <- do.call(cbind, lapply(retval1, '[[', "proposed.reconstruction", drop = FALSE))
    
    reconstruction.distances <- 
      GetDistanceInfo(list.of.MAP.out = retval1,
                      sparse.assign = use.sparse.assign, 
                      use.forward.search = use.forward.search)

    # browser()
    # all.tested <- GetAllTestedTables(list.of.MAP.out = retval1)
    all.tested <- lapply(retval1, `[[`, "all.tested")
    
    alt.solutions <- lapply(retval1, `[[`, "alt.solutions") # delete getaltsolutions
    time.for.assignment <- GetTimeForMAPAssign(retval1)

    return(list(proposed.assignment          = proposed.assignment,
                proposed.reconstruction      = proposed.reconstruction,
                reconstruction.distances     = reconstruction.distances,
                all.tested                   = all.tested,
                alt.solutions                = alt.solutions,
                time.for.assignment          = time.for.assignment,
                error.messages               = error.messages))
  }


#' Retrieve exposure information from the output generated by running
#' \code{MAPAssignActivity1} or \code{SparseAssignActivity1} on multiple samples
#'
#' @keywords internal
GetExposureInfo <- function(list.of.MAP.out) {
  # browser()
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
      exposures <- x$proposed.assignment
      count <- matrix(exposures[, 1], ncol = nrow(exposures))
      colnames(count) <- rownames(exposures)
      return(as.data.frame(count))
  })
  retval1 <- t(do.call(dplyr::bind_rows, tmp))
  colnames(retval1) <- names(list.of.MAP.out)
  retval2 <- retval1[SortSigId(rownames(retval1)), , drop = FALSE]
  return(retval2)
}


#' Retrieve distance information from the output generated by running
#' \code{MAPAssignActivity1} or \code{SparseAssignActivity1} on multiple samples
#'
#' @keywords internal
GetDistanceInfo <- function(list.of.MAP.out, sparse.assign = FALSE,
                            use.forward.search = FALSE) {
  MAP.tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    distances <- x$reconstruction.distances
    values <- matrix(distances$proposed.assignment, ncol = nrow(distances))
    colnames(values) <- distances$method
    return(as.data.frame(values))
  })

  MAP.retval <- do.call(dplyr::bind_rows, MAP.tmp)
  rownames(MAP.retval) <- names(list.of.MAP.out)
  
  QP.tmp <- lapply(list.of.MAP.out, FUN = function(x) {
      distances <- x$reconstruction.distances
      values <- matrix(distances$QP.assignment, ncol = nrow(distances))
      colnames(values) <- distances$method
      return(as.data.frame(values))
  })
  QP.retval <- do.call(dplyr::bind_rows, QP.tmp)
  rownames(QP.retval) <- names(list.of.MAP.out)
  
  if (use.forward.search) {
    return(list(forward.search.distances = MAP.retval, QP.distances = QP.retval))
  } else if (sparse.assign) {
    return(list(sparse.assign.distances = MAP.retval, QP.distances = QP.retval))
  } else {
    return(list(MAP.distances = MAP.retval, QP.distances = QP.retval))
  }
  
}


#' Retrieve time.for.assignment from the output generated by running
#' \code{MAPAssignActivity1} on multiple samples
#'
#' @keywords internal
GetTimeForMAPAssign <- function(list.of.MAP.out) {
  tmp <- lapply(list.of.MAP.out, FUN = function(x) {
    # if (is.null(x$error.messages)) {
      xx <- x$time.for.assignment
      df <- data.frame(user = xx[1] + xx[4], # user.self + user.child
                       system = xx[2] + xx[5],
                       elapsed = xx[3])
      return(df)
    # } else {
    #   return(NULL)
    # }
  })
  # index.of.non.null <- sapply(tmp, FUN = Negate(is.null))
  # tmp1 <- tmp[index.of.non.null]
  
  tmp2 <- do.call("rbind", tmp)
  return(tmp2)
}

#' @keywords internal
SortSigId <- function(sig.id) {
  return(gtools::mixedsort(x = sig.id))
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
LowMutationCountThresh <- function(mut.type) {
  if (is.null(mut.type)) {
    return(-1)
  } else if (mut.type %in% c("SBS96", "SBS192")) {
    return(100)
  } else if (mut.type %in% c("DBS78", "ID")) {
    return(25)
  } else {
    warning("Returning -1 because of unknown mut.type in LowMutationCountThresh: ", mut.type)
    return(-1)
  }
}

#' @keywords internal
DropLowMutationSamples <- function(spectra) {
  spectra <- as.matrix(spectra)
  mut.type <- GetMutationType(spectra)
  
  thresh.value <- LowMutationCountThresh(mut.type)

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

  if (exists("opar")) {
    on.exit(graphics::par(opar))
  }

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

#' Get the numerical parts of identifiers.
#'
#' @param s A character vector.
#'
#' @return A vector, each element of which is the integer
#' corresponding to the first string of digits of an element of s.
#'
#' @details Not very sophisticated.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x<- c("SBS22", "SBS2", "SBS7b", "SBS7a")
#' NumFromId(x)
#' x[order(NumFromId(x))]
#' }
NumFromId<- function(s) {
  return(
    # If there are no numerical parts from s, then there will
    # be warning message: NAs introduced by coercion
    # This warning message can cause the wrapper function 
    # calling NumFromId() very slow
    suppressWarnings(as.numeric(
      sub("[^0123456789]*(\\d+).*", "\\1", s, perl = TRUE)))
    )
}
