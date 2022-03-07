#' Find a Maximum A Posteriori (MAP) assignment of signature exposures that explain one spectrum.
#' 
#' This function also can do sparse assignment by specifying \code{use.sparse.assign = TRUE}.
#' 
#' @param drop.low.mut.samples Whether to exclude low mutation samples from
#' the analysis. If \code{TRUE(default)}, samples with SBS total mutations less
#' than 100, DBS or ID total mutations less than 25 will be dropped.
#'
#' @keywords internal
#'
#' @inheritParams MAPAssignActivityInternal
#'
#' @return A list with the elements:
#'
#' * \code{proposed.assignment}: Proposed signature assignment for \code{spect}
#' with the highest MAP found. If \code{use.sparse.assign = TRUE}, this will
#' be the most sparse set of signatures that can plausibly explain \code{spect}.
#'    
#' * \code{proposed.reconstruction} :Reconstruction based on \code{MAP}. 
#' If \code{use.sparse.assign = TRUE}, this will be the reconstruction based on
#' sparse assignment.
#' 
#' * \code{reconstruction.distances}: Various distances and similarities
#' between \code{spect} and \code{proposed.reconstruction}.
#'
#' * \code{all.tested}: A \code{tibble} of all the search results.
#' 
#' * \code{alt.solutions}: A \code{tibble} showing all the alternative solutions
#' that are statistically as good as the \code{proposed.assignment} that can
#' plausibly reconstruct \code{spect}.
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
#' sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
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
           max.subsets             = 1000,
           use.sparse.assign       = FALSE,
           drop.low.mut.samples    = TRUE,
           use.sig.presence.test   = FALSE,
           q.thresh                = 0.05) {
    
    if (drop.low.mut.samples) {
      spect <- DropLowMutationSamples(spect)
    } else {
      spect <- spect
    }
    
    if (ncol(spect) == 0) {
      return(NullReturnForMAPAssignActivity1(msg = "No sample to analyse"))
    }
    
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
          seed                    = seed,
          use.sparse.assign       = use.sparse.assign,
          use.sig.presence.test   = use.sig.presence.test,
          q.thresh                = q.thresh))
      
      xx <- ListOfList2Tibble(MAPout)
      
      if (use.sparse.assign == FALSE) {
        best <- dplyr::arrange(xx, .data$MAP)[nrow(xx),  ]
      } else if (use.sparse.assign == TRUE) {
        best <- dplyr::arrange(xx, .data$df, .data$loglh.of.exp)[nrow(xx), ]
      }
      
      names.best <- names(best[["exp"]])
      best.exp <- best[["exp"]][[1]]
      if (is.null(names(best.exp))) {
        names(best.exp) <- names.best
      }
      MAP <- tibble::tibble(sig.id = names(best.exp), count = best.exp)
      
      MAP.recon <-
        ReconstructSpectrum(sigs, exp = best.exp, use.sig.names = TRUE)
      
      if (use.sparse.assign == FALSE) {
        # Internally set max.presence.proportion to be 0.99 in case there will be -Inf
        # for MAP distances
        max.presence.proportion <- 0.99
        sigs.presence.prop[sigs.presence.prop > max.presence.proportion] <- 
          max.presence.proportion
        MAP.distances <-
          DistanceMeasures(spect = spect, recon = MAP.recon, 
                           nbinom.size = m.opts$nbinom.size,
                           model = names(best.exp),
                           sigs.presence.prop = sigs.presence.prop,
                           likelihood.dist = m.opts$likelihood.dist,
                           signatures = sigs[, names(best.exp), drop = FALSE])
      } else if (use.sparse.assign == TRUE) {
        MAP.distances <-
          DistanceMeasuresSparse(spect = spect, recon = MAP.recon, 
                                 nbinom.size = m.opts$nbinom.size,
                                 likelihood.dist = m.opts$likelihood.dist,
                                 signatures = sigs[, names(best.exp), drop = FALSE])
      }
      
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
      
      all.tested <- TestAltSolutions(tibble = xx, 
                                     sparse.assign = use.sparse.assign)
      alt.solutions <- GetAltSolutions(tibble = all.tested,
                                       spectrum = spect,
                                       sigs = sigs, 
                                       mc.cores = max.mc.cores,
                                       sparse.assign = use.sparse.assign,
                                       wt.thresh = 0.95,
                                       q.thresh = 0.05)
      return(list(proposed.assignment          = exposure,
                  proposed.reconstruction      = MAP.recon,
                  reconstruction.distances     = MAP.distances,
                  all.tested                   = all.tested,
                  alt.solutions                = alt.solutions,
                  time.for.MAP.assign          = time.for.MAP.assign
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
         alt.solutions                 = NULL,
         time.for.MAP.assign           = time.for.MAP.assign,
         error.messages                = msg
         ))
}

#' Find a Maximum A Posteriori assignment of signature exposures for one spectrum.
#' 
#' This function also can do sparse assignment by specifying \code{use.sparse.assign = TRUE}.
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
#'
#' @param use.sparse.assign Whether to use sparse assignment. If \code{TRUE},
#'   arguments designed for Maximum A Posteriori assignment such as
#'   \code{sigs.presence.prop} will be ignored.
#'   
#' @param use.sig.presence.test Whether to use signature presence test first to
#'   filter out those signatures that are not needed in the reconstruction of
#'   the spectrum.

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
           seed                    = NULL,
           use.sparse.assign       = FALSE,
           use.sig.presence.test   = FALSE,
           q.thresh                = 0.05) {
    
    # Type checking
    if (missing(sigs)) stop("MAPAssignActivityInternal: sigs is NULL")
    
    if (use.sparse.assign == FALSE) {
      if (missing(sigs.presence.prop)) {
        stop("MAPAssignActivityInternal: sigs.presence.prop is NULL")
      }
    }
    
    if (!is.null(seed)) set.seed(seed, kind = "L'Ecuyer-CMRG")
    
    if (use.sparse.assign) {
      msg <- "SparseAssignActivity1: "
    } else {
      msg <- "MAPAssignActivity1: "
    }
    
    my.msg <- function(trace.level, ...)
      if (m.opts$trace >= trace.level) message(msg, ...)
    
    if (use.sparse.assign == FALSE) {
      sigs.presence.prop[sigs.presence.prop > max.presence.proportion] <-
        max.presence.proportion
      sigs <- sigs[ , names(sigs.presence.prop), drop = FALSE]
    }
    
    max.sig.index <- ncol(sigs)
    my.msg(10, "number of signatures = ", max.sig.index)
    mode(spect) <-  'numeric'
    
    if (use.sparse.assign == FALSE) {
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
    }
    
    cannot.generate <- setdiff(which(rowSums(sigs) == 0), which(as.vector(spect) == 0))
    if (length(cannot.generate) > 0)
      stop("Cannot generate spectrum from the specified signatures;\n",
           "No signatures has > 0 proportion for ",
           paste(rownames(sigs)[cannot.generate], collapse = ", ")) 
    
    message("Analyzing sample ", colnames(spect))
    
    if (use.sig.presence.test) {
      sigs.presence.tests <- parallel::mclapply(colnames(sigs), FUN = function(sig.name) {
        retval <- SignaturePresenceTest1(spectrum = spect, 
                                         sigs = sigs, 
                                         target.sig.index = sig.name, 
                                         m.opts = m.opts, 
                                         seed = seed)
        return(retval)
      }, mc.cores = max.mc.cores)
      
      names(sigs.presence.tests) <- colnames(sigs)
      
      p.values <- sapply(sigs.presence.tests, FUN = "[[", 4)
      # Those signatures with q.values < q.thresh are needed in the reconstruction
      needed.sigs <- names(p.values[p.values < q.thresh])
      
      print(p.values)
      
      two.sigs <- c("SBS18", "SBS24")
      if (all(two.sigs %in% needed.sigs)) {
        if (!"SBS29" %in% needed.sigs) {
          needed.sigs <- c(needed.sigs, "SBS29")
        }
      }
      my.msg(10, "Remained signatures after signature presence test ", paste(needed.sigs, collapse = " "))
      sigs <- sigs[, needed.sigs, drop = FALSE]
    }
    
    start <- OptimizeExposure(spect, sigs, m.opts  = m.opts)
    
    lh.w.all <- start$loglh  # The likelihood with all signatures
    if (lh.w.all == -Inf) stop("Cannot generate spectrum from the specified signatures")
    
    start.exp <- start$exposure
    non.0.exp.index <- which(start.exp >= 0.5) # indices of signatures with non-zero
    # exposures
    
    if (length(non.0.exp.index) < ncol(sigs)) {
      removed.sig.names <- colnames(sigs)[which(start.exp < 0.5)]
      my.msg(1, ncol(sigs) - length(non.0.exp.index),
             " signatures removed at beginning")
      my.msg(1, "removed: ",
             paste(removed.sig.names, collapse = ", "))
    }
    
    df0.sig.names <- colnames(sigs)[non.0.exp.index]
    my.msg(1, "Starting with ",
           paste(df0.sig.names, collapse = ","),
           "\nmax.level = ", max.level,
           "\nlog likelihood using all signatures = ", lh.w.all)
    
    ss <- list(sig.names            = paste(df0.sig.names, collapse = ","),
               p.for.sig.subset     = NA,
               exp                  = start.exp[non.0.exp.index],
               loglh.of.exp         = lh.w.all)
    
    if (use.sparse.assign == FALSE) {
      df0.prob.of.model <- P.of.M(df0.sig.names, sigs.presence.prop)
      
      df0 <- c(ss, list(prob.of.model        = df0.prob.of.model,
                        MAP                  = lh.w.all + df0.prob.of.model,
                        df                   = 0))
    } else if (use.sparse.assign == TRUE) {
      df0 <- c(ss, list(df = 0))
    }
    
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
      
      rr <- list(
        sig.names            = paste(colnames(sigs)[try.sigs.indices],
                                     collapse = ","),
        p.for.sig.subset     = chisq.p,
        exp                  = try.exp[["exposure"]],
        #sig.indices          = try.sigs.indices,
        #removed.sig.names    = paste(colnames(sigs)[subset.to.remove.v],
        #                             collapse = ","),
        loglh.of.exp         = try.exp[["loglh"]]
      ) 
      if (!use.sparse.assign) {
        prob.of.model <- P.of.M(try.sigs.names, sigs.presence.prop)
        return(
          c(rr,
            list(
              prob.of.model        = prob.of.model,
              MAP                  = try.exp[["loglh"]] + prob.of.model,
              df                   = df)))
      } else if (use.sparse.assign == TRUE) {
        return(c(rr, list(df = df)))
      }
      
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
      
      check.mclapply.result(check.to.remove, "MAPAssignActivityInternal")
      
      p.to.remove <- unlist(lapply(check.to.remove, `[`, "p.for.sig.subset"))
      names(p.to.remove) <-
        unlist(lapply(check.to.remove,
                      function(x) paste(x$removed.sig.names, collapse = ",")))
      if (m.opts$trace > 10) {
        message(msg, "p.to.remove = ")
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
#' This is the second factor in the product \eqn{P(M|D) = P(D|M)P(M)}.
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
  sigs.not.having.prop <- setdiff(model, names(sigs.presence.prop))
  if (length(sigs.not.having.prop) > 0) {
    stop("Signatures used in the model but don't have presence proportions: ",
         paste(sigs.not.having.prop, collapse = " "))
  }
  present <- sigs.presence.prop[model]
  not.present.names <- setdiff(names(sigs.presence.prop), model)
  not.present <- sigs.presence.prop[not.present.names]
  not.present <- 1 - not.present
  rr <- sum(log(c(present, not.present)))
  # This could be -Inf if one of the signatures is always present in the prior data but absent; probably want to fudge this slightly
  return(rr)
}

#' Calculate several distance measures between a spectrum and its reconstruction.
#' 
#' The idea is to provide several measures of how well the 
#' reconstruction matches the spectrum.
#' 
#' @param spect The spectrum we are trying to reconstruct.
#' 
#' @param recon The unrounded reconstruction.
#' 
#' @param nbinom.size \strong{Only} needed when \code{likelihood.dist =
#'   "neg.binom"}.The dispersion parameter for the negative binomial
#'   distribution; smaller is more dispersed. See
#'   \code{\link[stats]{NegBinomial}}.
#'        
#' @param model Names of sigs present in a trial exposure. Do not use indices.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'   signature. A numerical vector (values between 0 and 1), with names being a
#'   superset of \code{model}.
#'   
#' @param likelihood.dist The probability distribution used to calculate the
#'   likelihood, can be either "multinom" (multinomial distribution) or
#'   "neg.binom" (negative binomial distribution).
#'   
#' @param signatures \strong{Only} used to compute distances for quadratic
#'   programming assignment. Signature as a matrix or data frame, with each row
#'   one mutation type (g.e. CCT > CAT or CC > TT) and each column a signature.
#'   It should be the proposed signatures used by Maximum A Posteriori (MAP)
#'   assignment to reconstruct \code{spect}. It is needed to calculate distances
#'   of quadratic programming assignment.
#'   
#' @return A data frame whose first column indicates the distance method. The
#'   second column \code{proposed.assignment} shows the values of various
#'   distances using Maximum A Posteriori (MAP) assignment. 
#'   
#'   When \code{signatures} is \strong{not} NULL, there will be a third column
#'   \code{QP.assignment} shows the values of various distances using quadratic
#'   programming assignment.
#'   
#' @keywords internal

DistanceMeasures <- 
  function(spect, recon, nbinom.size, model, sigs.presence.prop, 
           likelihood.dist = "multinom", signatures = NULL) {
    my.fn <- function(method, spect, recon) {
      df <- rbind(as.vector(spect),
                  as.vector(recon))
      return(suppressMessages(philentropy::distance(x = df, 
                                                    method = method, 
                                                    test.na = FALSE)))
    }
    
    vv <- unlist(lapply(c("euclidean", "manhattan","cosine"), my.fn,
                        spect = spect, recon = recon))
    if (likelihood.dist == "multinom") {
      log.likelihood <- LLHSpectrumMultinom(as.vector(spect), as.vector(recon))
    } else if (likelihood.dist == "neg.binom") {
      log.likelihood <- 
        LLHSpectrumNegBinom(as.vector(spect), as.vector(recon), nbinom.size = nbinom.size)
    }
    
    MAP <- LLHSpectrumMAP(spectrum = spect, expected.counts = recon,
                          nbinom.size = nbinom.size, model = model,
                          likelihood.dist = likelihood.dist,
                          sigs.presence.prop = sigs.presence.prop)
    
    vv <- c(log.likelihood = log.likelihood, MAP = MAP, vv)
    
    if (!is.null(signatures)) {
      # Do signature assignment using QP
      QP.expo <- OptimizeExposureQP(spectrum = spect, signatures = signatures)
      QP.expo.non.zero <- QP.expo[QP.expo >= 0.5]
      QP.recon <- ReconstructSpectrum(sigs = signatures, exp = QP.expo.non.zero,
                                      use.sig.names = TRUE)
      QP.distances <- unlist(lapply(c("euclidean", "manhattan","cosine"), my.fn,
                                    spect = spect, recon = QP.recon))
      
      if (likelihood.dist == "multinom") {
        log.likelihood <- LLHSpectrumMultinom(as.vector(spect), as.vector(QP.recon))
      } else if (likelihood.dist == "neg.binom") {
        log.likelihood <- 
          LLHSpectrumNegBinom(as.vector(spect), as.vector(QP.recon), nbinom.size = nbinom.size)
      }
      
      MAP <- LLHSpectrumMAP(spectrum = spect, expected.counts = QP.recon,
                            nbinom.size = nbinom.size, model = model,
                            likelihood.dist = likelihood.dist,
                            sigs.presence.prop = sigs.presence.prop)
      
      QP.distances <- c(log.likelihood = log.likelihood, MAP = MAP, QP.distances)
      
      return(tibble::tibble(method = names(vv), proposed.assignment = vv,
                            QP.assignment = QP.distances))
    } else {
      return(tibble::tibble(method = names(vv), proposed.assignment = vv))
    }
  }

#' Test to try to find alternative solutions that are statistically as good as
#' the best solution from \code{MAPAssignActivitiy1}
#' 
#' @keywords internal
TestAltSolutions <- function(tibble, sparse.assign = FALSE) {
  all.tested <- tibble
  
  if (nrow(all.tested) == 1) {
    return(all.tested)
  }
  
  if (sparse.assign) {
    loglh.col.name <- "loglh.of.exp"
    all.tested.sorted <- 
      dplyr::arrange(all.tested, dplyr::desc(df), dplyr::desc(loglh.of.exp))
  } else {
    loglh.col.name <- "MAP"
    all.tested.sorted <- 
      dplyr::arrange(all.tested, dplyr::desc(MAP))
  }
  
  best.df <- all.tested.sorted$df[1]
  best.loglh <- all.tested.sorted[[loglh.col.name]][1]
  best.sig.names <- unlist(strsplit(all.tested.sorted$sig.names[1], ","))
  best.num.of.sigs <- length(best.sig.names)
  best.AIC <- -2 * best.loglh + 2 * best.num.of.sigs
  best.set <- sets::as.set(best.sig.names)
  
  # Get all the non nested solution indices
  non.nested.indices0 <- sapply(2:nrow(all.tested.sorted), FUN = function(x) {
    index <- x
    try.df <- all.tested.sorted$df[index]
    try.sig.names <- unlist(strsplit(all.tested.sorted$sig.names[index], ","))
    # If the model has the same df as the best solution, then it is non nested
    if (try.df == best.df) {
      return(index)
    } else {
      # If the model is not a superset of the best solution, then it is non nested
      try.set <- sets::as.set(try.sig.names)
      if (!sets::set_is_proper_subset(x = best.set, y = try.set)) {
        return(index)
      }
    }
  })
  non.nested.indices <- unlist(non.nested.indices0)
  nested.indices <- setdiff(2:nrow(all.tested.sorted), non.nested.indices)
  
  # Check whether the non nested models are statistically different from the 
  # best solution by Akaike weights
  # https://link.springer.com/content/pdf/10.3758/BF03206482.pdf
  non.nested.results <- sapply(non.nested.indices, FUN = function(x) {
    index <- x
    try.loglh <- all.tested.sorted[[loglh.col.name]][index]
    try.sig.names <- unlist(strsplit(all.tested.sorted$sig.names[index], ","))
    try.num.of.sigs <- length(try.sig.names)
    try.AIC <- -2 * try.loglh + 2 * try.num.of.sigs
    AIC.min <- min(best.AIC, try.AIC)
    
    # Calculate relative likelihood
    best.rv.lh <- exp((AIC.min - best.AIC) / 2)
    try.rv.lh <- exp((AIC.min - try.AIC) / 2)
    
    # Calculate Akaike weights
    try.akaike.weights <- try.rv.lh / sum(best.rv.lh + try.rv.lh)
    
    # If the Akaike weights of the alternative solution is greater than 0.95, then
    # it is significantly favored over the best solution
  })
  
  # Perform likelihood ratio test for nested models
  nested.results <- sapply(nested.indices, FUN = function(x) {
    index <- x
    try.df <- all.tested.sorted$df[index]
    try.loglh <- all.tested.sorted[[loglh.col.name]][index]
    
    LR.statistic <- 2 * (try.loglh - best.loglh)
    df <- best.df - try.df
    p.value <- stats::pchisq(q = LR.statistic, df = df, lower.tail = FALSE)
  })
  
  all.tested.sorted$akaike.weights <- NA
  all.tested.sorted$LRT.p.value <- NA
  all.tested.sorted$LRT.q.value <- NA
  
  df <- as.data.frame(all.tested.sorted)
  
  df[non.nested.indices, ]$akaike.weights <- non.nested.results
  df[nested.indices, ]$LRT.p.value <- nested.results
  # Adjust p-values for multiple comparisons to Benjamini-Hochberg false discovery rate
  df[nested.indices, ]$LRT.q.value <- 
    stats::p.adjust(p = nested.results, method = "BH") 
  return(df)
}

#' Get alternative solutions that are statistically as good as the best solution
#' from \code{MAPAssignActivitiy1} 
#'
#' @keywords internal
GetAltSolutions <- function(tibble, spectrum, sigs, mc.cores = 1, 
                            sparse.assign = FALSE, 
                            wt.thresh = 0.95, q.thresh = 0.05) {
  all.tested <- tibble
  
  if (nrow(all.tested) == 1) {
    return(all.tested[0, ])
  }
  
  all.tested.nested.models <- all.tested[!is.na(all.tested$LRT.q.value), ]
  all.tested.LR.OK <- 
    all.tested.nested.models[all.tested.nested.models$LRT.q.value < q.thresh, ]
  
  all.test.non.nested.models <- all.tested[!is.na(all.tested$akaike.weights), ]
  all.tested.akaike.weights.OK <-
    all.test.non.nested.models[all.test.non.nested.models$akaike.weights > wt.thresh, ]
  
  alt.solutions <- rbind(all.tested.LR.OK, all.tested.akaike.weights.OK)
  if (sparse.assign) {
    alt.solutions1 <- 
      dplyr::arrange(alt.solutions, dplyr::desc(df), dplyr::desc(loglh.of.exp))
  } else {
    alt.solutions1 <- dplyr::arrange(alt.solutions, dplyr::desc(MAP))
  }
  
  if (nrow(alt.solutions1) == 0) {
    return(alt.solutions1)
  }
  
  retval <- parallel::mclapply(1:nrow(alt.solutions1), FUN = function(x) {
    index <- x
    sig.names <- unlist(strsplit(alt.solutions1$sig.names[index], ","))
    sig.to.use <- sigs[, sig.names, drop = FALSE]
    QP.retval <- OptimizeExposureQP(spectrum = spectrum, signatures = sig.to.use)
    reconstuction <- ReconstructSpectrum(sigs = sig.to.use, exp = QP.retval)
    cosine <- cossim(v1 = spectrum, v2 = reconstuction)
    return(tibble::tibble(QP.exp = list(QP.retval), QP.cosine = cosine))
  }, mc.cores = Adj.mc.cores(mc.cores))
  
  retval2 <- do.call("rbind", retval)
  alt.solutions2 <- cbind(alt.solutions1, retval2)
  return(alt.solutions2)
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

