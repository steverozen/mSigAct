#' Find a Maximum A Posteriori (MAP) assignment of signature exposures that explain one spectrum.
#'
#' @export
#'
#' @inheritParams MAPAssignActivityInternal
#'
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: Proposed signature assignment for \code{spect}
#' with the highest MAP found.
#'    
#' * \code{proposed.reconstruction} :Reconstruction based on \code{MAP}.
#' 
#' * \code{reconstruction.distances}: Various distances and similarities
#' between \code{spect} and \code{proposed.reconstruction}.
#'
#' * \code{all.tested}: A \code{tibble} of all the search results.
#' 
#' * \code{time.for.MAP.assign}: Value from \code{system.time} for running
#'  \code{MAPAssignActivity1}.
#'
#' * \code{error.messages}: Only present if there were errors running
#' \code{MAPAssignActivity1}.
#'
#' The elements \code{proposed.assignment}, \code{proposed.reconstruction},
#' \code{reconstruction.distances}, \code{all.tested},
#' \code{time.for.MAP.assign} will be \code{NULL} if the algorithm could not
#' find the optimal reconstruction or there are errors coming out.
#' 
#' @md
#' 
#' @examples 
#' \dontrun{
#' # This is a long running example unless parallel computing is supported on your machine
#' indices <- grep("Lung-AdenoCA", colnames(PCAWG7::spectra$PCAWG$SBS96))
#' spect <- PCAWG7::spectra$PCAWG$SBS96[, indices[1], drop = FALSE]
#' sigs <- PCAWG7::signature$genome$SBS96
#' sigs.prop <- ExposureProportions(mutation.type = "SBS96", 
#'                                  cancer.type = "Lung-AdenoCA")
#' MAP.out <- MAPAssignActivity1(spect = spect, 
#'                               sigs = sigs, 
#'                               sigs.presence.prop = sigs.prop, 
#'                               max.level = length(sigs.prop) - 1)
#' }                               
MAPAssignActivity1 <-
  function(spect,
           sigs,
           sigs.presence.prop,
           max.level               = 5,
           p.thresh                = 0.05,
           m.opts                  = DefaultManyOpts(),
           max.mc.cores            = min(20, 2^max.level),
           progress.monitor        = NULL,
           seed                    = NULL,
           max.subsets             = 1000) {
    time.for.MAP.assign <- system.time(3)
    
    tryCatch({
      
      if (sum(spect) < 1) stop("< 1 mutation in spectrum ", colnames(spect))
      
      # If there are non integers in spect, round it first. Otherwise, there
      # will be a lot of warnings later calculating likelihood using
      # stats::dnbinom() and this will cause the program run very slowly
      if (any(spect != round(spect))) {
        warning("Round non integers in spectrum to integers")
        spect <- round(spect)
      }

      time.for.MAP.assign <- system.time(
        MAPout <- MAPAssignActivityInternal(
          spect                   = spect,
          sigs                    = sigs,
          sigs.presence.prop      = sigs.presence.prop,
          max.level               = max.level,
          p.thresh                = p.thresh,
          m.opts                  = m.opts,
          max.mc.cores            = max.mc.cores,
          max.subsets             = max.subsets,
          max.presence.proportion = 0.99,
          progress.monitor        = progress.monitor,
          seed                    = seed))
      
      xx <- ListOfList2Tibble(MAPout)
      
      best <- dplyr::arrange(xx, .data$MAP)[nrow(xx),  ]
      names.best <- names(best[["exp"]])
      best.exp <- best[["exp"]][[1]]
      if (is.null(names(best.exp))) {
        names(best.exp) <- names.best
      }
      MAP <- tibble::tibble(sig.id = names(best.exp), count = best.exp)
      
      if (FALSE) {
        sparse.best <-
          dplyr::arrange(xx, .data$df, .data$MAP)[nrow(xx), ]
        names.sparse.best <- names(sparse.best[["exp"]])
        most.sparse.exp <- sparse.best[["exp"]][[1]]
        if (is.null(names(most.sparse.exp))) {
          names(most.sparse.exp) <- names.sparse.best # Necessary if only 1 signature
        }
        
        most.sparse <-
          tibble::tibble(sig.id = names(most.sparse.exp), count = most.sparse.exp)
      }
      
      # Best MAP
      MAP.recon <-
        ReconstructSpectrum(sigs, exp = best.exp, use.sig.names = TRUE)
      
      # Internally set max.presence.proportion to be 0.99 in case there will be -Inf
      # for MAP distances
      max.presence.proportion <- 0.99
      sigs.presence.prop[sigs.presence.prop > max.presence.proportion] <- 
        max.presence.proportion
      MAP.distances <-
        DistanceMeasures(spect = spect, recon = MAP.recon, 
                         nbinom.size = m.opts$nbinom.size,
                         model = names(best.exp),
                         sigs.presence.prop = sigs.presence.prop)
      
      # Round the exposure and reconstruction
      exposure <- matrix(round(MAP$count), nrow = nrow(MAP))
      rownames(exposure) <- MAP$sig.id
      colnames(exposure) <- colnames(spect)
      
      # If there are signatures get assigned zero mutation counts after rounding,
      # remove these signatures from the exposure matrix
      non.zero.indices <- rowSums(exposure) > 0
      exposure <- exposure[non.zero.indices, , drop = FALSE]
      
      MAP.recon <- round(MAP.recon)
      colnames(MAP.recon) <- colnames(spect)
      
      # Add attributes to MAP.recon to be same as spect
      MAP.recon <- AddAttributes(MAP.recon, spect)
      
      if (FALSE) {
        # MAP most sparse
        sparse.MAP.recon <-
          ReconstructSpectrum(sigs, exp = most.sparse.exp, use.sig.names = TRUE)
        
        sparse.MAP.distances <-
          DistanceMeasures(spect, sparse.MAP.recon, m.opts$nbinom.size)
      }
      
    return(list(proposed.assignment          = exposure,
                proposed.reconstruction      = MAP.recon,
                reconstruction.distances     = MAP.distances,
                all.tested                   = xx,
                time.for.MAP.assign          = time.for.MAP.assign
                #MAP.row                     = best,
                #best.sparse                 = most.sparse,
                #best.sparse.row             = sparse.best,
                #error.messages              = c(),
                #success                     = TRUE,
                #sparse.MAP.recon            = sparse.MAP.recon,
                #sparse.MAP.distances        = sparse.MAP.distances
                ))
    },
    error = function(err.info) {
      if (!is.null(err.info$message)) err.info <- err.info$message
      message(err.info)
      return(NullReturnForMAPAssignActivity1(err.info, time.for.MAP.assign))
    })
  }

NullReturnForMAPAssignActivity1 <- function(msg, all.tested, 
                                            time.for.MAP.assign = NULL) {
  return(
    list(proposed.assignment           = NULL,
         proposed.reconstruction       = NULL,
         reconstruction.distances      = NULL,
         all.tested                    = NULL,
         time.for.MAP.assign           = time.for.MAP.assign,
         error.messages                = msg
         #success                      = FALSE,
         #MAP.row                      = NULL,
         #best.sparse                  = NULL,
         #best.sparse.row              = NULL,
         #sparse.MAP.recon             = NULL,
         #sparse.MAP.distances         = NULL
         ))
}

#' Find a Maximum A Posteriori assignment of signature exposures for one spectrum.
#'
#' @keywords internal
#'
#' @param spect A single spectrum.
#'
#' @param sigs A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being a subset of \code{colnames(sigs)}. See \code{\link{ExposureProportions}}
#'    for more details.
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param p.thresh If
#'  the p value for a better reconstruction with a
#'  set of signatures (as opposed to
#'  without that set of signatures)
#'  is > than this argument, then we can use exposures without this set.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.
#'
#' @param max.subsets This argument provides a way to 
#'   heuristically limit the
#'   amount of time spent by this function. Larger values of this
#'   argument will tend to allow longer running times.
#'   The algorithm
#'   successively tries to remove all subsets of 1 signature, 2
#'   signatures, 3 signatures, etc., down to \code{max.level}.
#'   (Not every subset is tested at each level; if a subset was
#'   already found to be necessary the algorithm does not test
#'   supersets of that subset.) If at any level the algorithm
#'   needs to test more than \code{max.subsets} this function will
#'   not proceed.
#'
#' @param max.presence.proportion The maximum value of the proportion
#'   of tumors that must have a given signature.
#'
#' @param progress.monitor Function called at the start of each
#'   new level (number of signatures to try excluding). Must
#'   take named arguments \code{value} and \code{detail}, and
#'   no others. Designed for a \code{\link[ipc]{AsyncProgress}}
#'   progress bar function.
#'   
#' @param seed Random seed; set this to get reproducible results. (The
#'   numerical optimization is in two phases; the first, global phase
#'   might rarely find different optima depending on the random
#'   seed.)

MAPAssignActivityInternal <-
  function(spect,
           sigs,
           sigs.presence.prop,
           max.level               = 5,
           p.thresh                = 0.05,
           m.opts                  = DefaultManyOpts(),
           max.mc.cores            = min(20, 2^max.level),
           max.subsets             = 1000,
           max.presence.proportion = 0.99,
           progress.monitor        = NULL,
           seed                    = NULL) {
  
  # Type checking
  if (is.null(sigs)) stop("MAPAssignActivityInternal: sigs is NULL")
  if (is.null(sigs.presence.prop)) 
    stop("MAPAssignActivityInternal: sigs.presence.prop is NULL")
  
  if (!is.null(seed)) set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  my.msg <- function(trace.level, ...)
    if (m.opts$trace >= trace.level) message("MAPAssignActivity1: ", ...)

  sigs.presence.prop[sigs.presence.prop > max.presence.proportion] <-
    max.presence.proportion

  sigs <- sigs[ , names(sigs.presence.prop), drop = FALSE]

  max.sig.index <- ncol(sigs)
  my.msg(10, "number of signatures = ", max.sig.index)
  mode(spect) <-  'numeric'
  if (!isTRUE(all.equal(colnames(sigs), names(sigs.presence.prop)))) {
    msg <- paste("class sigs =", 
    class(sigs), 
    "\nclass sigs.presence.prop =",
    class(sigs.presence.prop),
    "colnames(sigs) =",
    colnames(sigs),
    "names(sigs.presence.prop) =",
    names(sigs.presence.prop),
    "!isTRUE(all.equal(colnames(sigs), names(sigs.presence.prop)))")
  }
 
  cannot.generate <- setdiff(which(rowSums(sigs) == 0), which(as.vector(spect) == 0))
  if (length(cannot.generate) > 0)
    stop("Cannot generate spectrum from the specified signatures;\n",
         "No signatures has > 0 proportion for ",
         paste(rownames(sigs)[cannot.generate], collapse = ", ")) 
  
  message("Analyzing sample ", colnames(spect))
  start <- OptimizeExposure(spect, sigs, m.opts  = m.opts)

  lh.w.all <- start$loglh  # The likelihood with all signatures
  if (lh.w.all == -Inf) stop("Cannot generate spectrum from the specified signatures")

  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp >= 0.5) # indices of signatures with non-zero
                                             # exposures

  if (length(non.0.exp.index) < ncol(sigs)) {
    removed.sig.names <- colnames(sigs)[which(start.exp < 0.5)]
    my.msg(0, ncol(sigs) - length(non.0.exp.index),
           " signatures removed at beginning")
    my.msg(0, "removed: ",
           paste(removed.sig.names, collapse = ", "))
  }

  df0.sig.names <- colnames(sigs)[non.0.exp.index]
  my.msg(0, "Starting with ",
         paste(df0.sig.names, collapse = ","),
         "\nmax.level = ", max.level,
         "\nlog likelihood using all signatures = ", lh.w.all)

  df0.prob.of.model <- P.of.M(df0.sig.names, sigs.presence.prop)

  df0 <- list(sig.names            = paste(df0.sig.names, collapse = ","),
              p.for.removing.sigs  = NA,
              exp                  = start.exp[non.0.exp.index],
              sig.indices          = non.0.exp.index,
              removed.sig.names    = "",
              loglh.of.exp         = lh.w.all,
              prob.of.model        = df0.prob.of.model,
              MAP                  = lh.w.all + df0.prob.of.model,
              df                   = 0)

  out.list <- list(df0)

  if (length(non.0.exp.index) < 2) {
    my.msg(0, "returning, only ", length(non.0.exp.index), " non-0 exposures")
    return(out.list)
  }

  max.level <- min(max.level, length(non.0.exp.index) - 1)

  # c.r.s is "cannot remove subsets", i.e. subset of the signature indices that
  # cannot be removed
  c.r.s <- sets::set()

  # Internal function ====================================================
  info.of.removed.subset <- function(subset) {
    # subset is of class set (package sets)
    subset.to.remove.v <- as.numeric(subset)
    if (any(subset.to.remove.v > max.sig.index)) {
      stop("Got illegal signature index in subset.to.remove = ",
           paste(subset.to.remove.v, collapse = " "))
    }
    try.sigs.indices <- setdiff(non.0.exp.index, subset.to.remove.v)
    try.sigs <- sigs[ , try.sigs.indices, drop = FALSE]

    my.msg(10, "Testing removal of signature index(s) ",
           paste(subset.to.remove.v, collapse = ", "),
           "; name(s) = ",
           paste(colnames(sigs)[subset.to.remove.v], collapse = ", "))

    try.sigs.names <- colnames(sigs)[try.sigs.indices]

    # Get the maximum likelihood exposure for try.sigs
    try.exp <- OptimizeExposure(spect,
                                try.sigs,
                                m.opts = m.opts)

    if (is.infinite(try.exp$loglh)) {
      # It was not possible to generate the spectrum from the signatures, so
      # try.exp$loglh should be -Inf and the p value the test that the spectrum
      # can be better reconstructed with union(try.sigs, subset.to.remove.v)
      # than with try.sigs is 0.
       if (try.exp$loglh < 0) {
          chisq.p <- 0
       } else {
         stop("Probable coding error")
       }
    } else {
      statistic <- 2 * (lh.w.all - try.exp$loglh)
      chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
    }
    prob.of.model <- P.of.M(try.sigs.names, sigs.presence.prop)

    return(list(
      sig.names            = paste(colnames(sigs)[try.sigs.indices],
                                   collapse = ","),
      p.for.removing.sigs  = chisq.p,
      exp                  = try.exp[["exposure"]],
      sig.indices          = try.sigs.indices,
      removed.sig.names    = paste(colnames(sigs)[subset.to.remove.v],
                                   collapse = ","),
      loglh.of.exp         = try.exp[["loglh"]],
      prob.of.model        = prob.of.model,
      MAP                  = try.exp[["loglh"]] + prob.of.model,
      df                   = df))
  } # Internal function info.of.removed.subset ===========================


  for (df in 1:max.level) {
    on.exit(
      if (!is.null(progress.monitor)) {
        progress.monitor(
          value  = 1 - (df - 1) / max.level,
          detail = "Finished testing removal of subsets of signatures")
      }
    )

    my.msg(0, "\ndf = ", df)
    subsets <- as.list(sets::set_combn(non.0.exp.index, df))
    discard <- lapply(subsets, is.superset.of.any, background = c.r.s)
    subsets2 <- subsets[!(unlist(discard))]
    if (length(subsets2) == 0) break;

    my.msg(0, "Number of subsets to remove = ", length(subsets2))
    if (length(subsets2) > max.subsets) {
      my.msg(-1, "Number of subsets (", length(subsets2), ") > max.subsets (", max.subsets, ")")
      my.msg(-1, "Returning NULL")
      stop("There were too many ways to reconstruct the spectrum ", 
           colnames(spect),
          "; please try removing some of the less likely signatures. ",
          "Or you can try to increase the value of argument max.subsets, ",
          "which will take longer computing time")
    }

    if (!is.null(progress.monitor)) {
      progress.monitor(
        value  = 1 / max.level,
        detail = paste0("Testing removal of subsets of ", df, " signatures (",
                        length(subsets2), " subsets)"))
    }

    time.used <- system.time(
      check.to.remove <-
        parallel::mclapply(X = subsets2,
                           FUN = info.of.removed.subset,
                           # Set to mc.cores 1 if OS is MS Windows
                           mc.cores = Adj.mc.cores(max.mc.cores)
        )
    )
    if (m.opts$trace > 3) {
      message("Time used to compute p values for df ", df, ":")
      print(time.used)
    }

    check.mclapply.result(check.to.remove, "SparseAssignActivity1")

    p.to.remove <- unlist(lapply(check.to.remove, `[`, "p.for.removing.sigs"))
    names(p.to.remove) <-
      unlist(lapply(check.to.remove,
                    function(x) paste(x$removed.sig.names, collapse = ",")))
    if (m.opts$trace > 10) {
      message("MAPAssignActivity1: p.to.remove = ")
      for (ii in 1:length(p.to.remove)) {
        message(ii, " ", names(p.to.remove)[ii], " ", p.to.remove[ii])
      }
    }
    if (all(p.to.remove < p.thresh)) {
      my.msg(10, "Cannot remove any subsets at level ", df)
      break;
    }
    ok.after.removal <- check.to.remove[p.to.remove >= p.thresh]
    out.list <-c(out.list, ok.after.removal)

    cannot.remove <- subsets2[p.to.remove < p.thresh]
    c.r.s <- sets::set_union(c.r.s, sets::as.set(cannot.remove))

  }

  return(out.list)
}

#' Calculate \eqn{P(M)} -- the probability of a model of which signatures are present in a sample.
#'
#' @param model Names of sigs present in a trial exposure. Do not use indices.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being a superset of \code{model}.
#'
#' @keywords internal
#'
P.of.M <- function(model, sigs.presence.prop) {
  stopifnot(length(setdiff(model, names(sigs.presence.prop))) == 0)
  present <- sigs.presence.prop[model]
  not.present.names <- setdiff(names(sigs.presence.prop), model)
  not.present <- sigs.presence.prop[not.present.names]
  not.present <- 1 - not.present
  rr <- sum(log(c(present, not.present)))
  # This could be -Inf if one of the signatures is always present in the prior data but absent; probably want to fudge this slightly
  return(rr)
}

#' Calculate several distance measures between a spectrum and its reconstruction
#'
#' @keywords internal

DistanceMeasures <- 
  function(spect, recon, nbinom.size, model, sigs.presence.prop) {
    my.fn <- function(method) {
      df <- rbind(as.vector(spect),
                  as.vector(recon))
      return(philentropy::distance(x = df, method = method, test.na = FALSE))
    }
    
    vv <- unlist(lapply(c("euclidean", "manhattan","cosine"), my.fn))
    vv <- c(neg.log.likelihood =
              LLHSpectrumNegBinom(
                as.vector(spect),
                as.vector(recon),
                nbinom.size = nbinom.size),
            MAP = 
              LLHSpectrumMAP(spectrum = spect, expected.counts = recon,
                             nbinom.size = nbinom.size, model = model,
                             sigs.presence.prop = sigs.presence.prop),
            vv)
    return(tibble::tibble(method = names(vv), value = vv))
  }

#' @return A list with the elements \describe{
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

