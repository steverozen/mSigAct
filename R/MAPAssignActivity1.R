#' Find a Maximum A Posteriori (MAP) assignment of signature exposures that explain one spectrum.
#'
#' @export
#'
#' @inheritParams MAPAssignActivityInternal
#'
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
#'
#' }
#'
#' These elements will be \code{NULL} if \code{max.subsets} is exceeded.

MAPAssignActivity1 <-
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
    
    if (sum(spect) < 1) {
      return(NullReturnForMAPAssignActivity1("0 mutations in spectrum"))
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
      max.presence.proportion = max.presence.proportion,
      progress.monitor        = progress.monitor,
      seed                    = seed))

  if (is.null(MAPout)) {
    return(
      NullReturnForMAPAssignActivity1(
        paste("There were too many ways to reconstruct the spectrum; ",
              "please try removing some of the less likely signatures"), 
        time.for.MAP.assign))
  }
  
  xx <- ListOfList2Tibble(MAPout)

  best <- dplyr::arrange(xx, .data$MAP)[nrow(xx),  ]
  names.best <- names(best[["exp"]])
  best.exp <- best[["exp"]][[1]]
  if (is.null(names(best.exp))) {
    names(best.exp) <- names.best
  }
  MAP <- tibble::tibble(sig.id = names(best.exp), count = best.exp)

  sparse.best <-
    dplyr::arrange(xx, .data$df, .data$MAP)[nrow(xx), ]
  names.sparse.best <- names(sparse.best[["exp"]])
  most.sparse.exp <- sparse.best[["exp"]][[1]]
  if (is.null(names(most.sparse.exp))) {
    names(most.sparse.exp) <- names.sparse.best # Necessary if only 1 signature
  }

  most.sparse <-
    tibble::tibble(sig.id = names(most.sparse.exp), count = most.sparse.exp)

  # Best MAP
  MAP.recon <-
    ReconstructSpectrum(sigs, exp = best.exp, use.sig.names = TRUE)

  # MAP most sparse
  sparse.MAP.recon <-
    ReconstructSpectrum(sigs, exp = most.sparse.exp, use.sig.names = TRUE)

  MAP.distances <-
    DistanceMeasures(spect, MAP.recon, m.opts$nbinom.size)

  sparse.MAP.distances <-
    DistanceMeasures(spect, sparse.MAP.recon, m.opts$nbinom.size)


  return(list(MAP                  = MAP,
              MAP.row              = best,
              best.sparse          = most.sparse,
              best.sparse.row      = sparse.best,
              all.tested           = xx,
              messages             = c(),
              success              = TRUE,
              time.for.MAP.assign  = time.for.MAP.assign,
              MAP.recon            = MAP.recon,
              sparse.MAP.recon     = sparse.MAP.recon,
              MAP.distances        = MAP.distances,
              sparse.MAP.distances = sparse.MAP.distances))
  }

NullReturnForMAPAssignActivity1 <- function(msg, time.for.MAP.assign = NULL) {
  return(
    list(MAP                  = NULL,
         MAP.row              = NULL,
         best.sparse          = NULL,
         best.sparse.row      = NULL,
         all.tested           = NULL,
         messages             = msg,
         success              = FALSE,
         time.for.MAP.assign  = time.for.MAP.assign,
         MAP.recon            = NULL,
         sparse.MAP.recon     = NULL,
         MAP.distances        = NULL,
         sparse.MAP.distances = NULL))
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
#'    being a subset of \code{colnames(sigs)}.
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param p.thresh If
#'  the p value for a better reconstruction with as opposed to
#'  without a set of signatures
#'  is > than this argument, then we can use exposures without this set.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.
#'
#' @param max.subsets The maximum number of subsets that can be
#'   tested for removal from the set of signatures.
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

MAPAssignActivityInternal <- function(spect,
                                      sigs,
                                      sigs.presence.prop,
                                      max.level    = 5,
                                      p.thresh     = 0.05,
                                      m.opts       = DefaultManyOpts(),
                                      max.mc.cores = min(20, 2^max.level),
                                      max.subsets  = 1000,
                                      max.presence.proportion = 0.99,
                                      progress.monitor  = NULL) {
  
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
  start <- OptimizeExposure(spect,
                            sigs,
                            m.opts  = m.opts)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  if (lh.w.all == -Inf) {
    problem <- "Cannot generate spect from sigs"
    message(problem)
    rr <- rep(NaN, ncol(sigs))
    attr(rr, "log.likelihood") <- lh.w.all
    attr(rr, "problem") <- problem
    return(rr)
  }
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

  mc.cores <- Adj.mc.cores(max.mc.cores) # Set to 1 if OS is MS Windows

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
         stop("How did we here?")
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
      return(NULL)
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
                           mc.cores = max.mc.cores)
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

DistanceMeasures <- function(spect, recon, nbinom.size) {
  my.fn <- function(method) {
    df <- rbind(as.vector(spect),
                as.vector(recon))
    return(philentropy::distance(x = df, method = method, test.na = FALSE))
  }

  vv <- unlist(lapply(c("euclidean", "cosine"), my.fn))
  vv <- c(neg.log.likelihood =
            LLHSpectrumNegBinom(
              as.vector(spect),
              as.vector(recon),
              nbinom.size = nbinom.size),
          vv)
  return(tibble::tibble(method = names(vv), value = vv))
}

