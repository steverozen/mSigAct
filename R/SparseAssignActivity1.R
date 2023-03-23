if (FALSE) {
#' Component of \code{\link{SparseAssignActivity}} for one spectrum.
#' @keywords internal
#'
#' @param spect A single spectrum.
#'
#' @param sigs A numerical matrix, possibly an \code{\link[ICAMS]{ICAMS}} catalog.
#'
#' @param max.level The maximum number of signatures to try removing.
#'
#' @param p.thresh The p value threshold for deciding if a set of signatures is necessary.
#'
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#'
#' @param max.mc.cores
#'   The maximum number of cores to use.
#'   On Microsoft Windows machines it is silently changed to 1.)
#'   
#' @param seed Random seed; set this to get reproducible results. (The
#'   numerical optimization is in two phases; the first, global phase
#'   might rarely find different optima depending on the random
#'   seed.)
#' @name tryingtoremove
SparseAssignActivity1 <- function(spect,
                                  sigs,
                                  max.level    = 5,
                                  p.thresh     = 0.05,
                                  m.opts       = DefaultManyOpts(),
                                  max.mc.cores = min(20, 2^max.level),
                                  seed         = NULL) {
  
  if (!is.null(seed)) set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  my.msg <- function(trace.level, ...)
    if (m.opts$trace >= trace.level) message("SparseAssignActivity1: ", ...)

  max.sig.index <- ncol(sigs)
  
  message("Analyzing sample ", colnames(spect))
  my.msg(10, "SparseAssignActivity1: number of signatures = ", max.sig.index)
  mode(spect) <-  'numeric'
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
  non.0.exp.index <- which(start.exp > 0.5) # indices of signatures with non-zero
                                            # exposures； TODO, possibly remove
  my.msg(1, "Starting with ",
         paste(names(start.exp)[non.0.exp.index], collapse = ","),
         "\nmax.level = ", max.level,
         "\nlog likelihood using all signatures = ", lh.w.all)

  if (length(non.0.exp.index) < 2) {
    my.msg(0, "returning, only ", length(non.0.exp.index), " non-0 exposures")
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)

  mc.cores <- Adj.mc.cores(max.mc.cores) # Set to 1 if OS is MS Windows

  # c.r.s is "cannot remove subsets", i.e. subset of the signature indices that
  # cannot be removed
  c.r.s <- sets::set()
  best.exp <- list() # Index by level, the expression of the assignments (exposures) at a given level
  best.sig.indices <- list() # Index by level, the signature indices for the exposures in best.exp
  # TODO, Should best actually be maximum likelihood or maximum a posteriori?

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

    # Get the maximum likelihood exposure for try.sigs
    try.exp <- OptimizeExposure(spect,
                                try.sigs,
                                m.opts = m.opts)

    # TODO -- deal with case when try$loglh is Inf

    statistic <- 2 * (lh.w.all - try.exp$loglh)
    chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
    return(list(p                 = chisq.p,
                exp               = try.exp[["exposure"]],
                sig.indices       = try.sigs.indices,
                removed.sig.names = colnames(sigs)[subset.to.remove.v],
                loglh.of.exp      = try.exp[["loglh"]]
                ))
  }

  for (df in 1:max.level) {
    my.msg(0, "\ndf = ", df)
    subsets <- as.list(sets::set_combn(non.0.exp.index, df))
    discard <- lapply(subsets, is.superset.of.any, background = c.r.s)
    subsets2 <- subsets[!(unlist(discard))]
    if (length(subsets2) == 0) break;

    my.msg(0, "Number of subsets to remove = ", length(subsets2))
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

    p.to.remove <- unlist(lapply(check.to.remove, `[`, "p"))
    names(p.to.remove) <-
      unlist(lapply(check.to.remove,
                    function(x) paste(x$removed.sig.names, collapse = ",")))
    if (m.opts$trace > 10) {
      message("SparseAssignActivity1: p.to.remove =")
      for (ii in 1:length(p.to.remove)) {
        message(ii, " ", names(p.to.remove)[ii], " ", p.to.remove[ii])
      }
    }
    if (all(p.to.remove < p.thresh)) {
      my.msg(10, "Cannot remove any subsets at level ", df)
      break;
    }
    cannot.remove <- subsets2[p.to.remove < p.thresh]
    c.r.s <- sets::set_union(c.r.s, sets::as.set(cannot.remove))
    xx <- which(p.to.remove == max(p.to.remove))
    my.msg(10, "xx = ", paste(xx, collapse = ","),
           "\nmax(p.to.remove = ", max(p.to.remove), ")")
    if (length(xx) > 1) {
      possible.sigs <-
        lapply(X = subsets[xx],
               function(subset) paste(colnames(sigs)[unlist(subset)], collapse = ","))
      xx <- min(xx)
      sig.name.to.remove <- colnames(sigs)[unlist(subsets[[xx]])]

      message("Temporary warning, > 1 signature can be removed; ",
              "selecting the first one arbitrarily (index = ", xx, ",
              name = ", paste(sig.name.to.remove, collapse = ","), ")\n",
              "Possibilities were\ ", paste(possible.sigs, collapse = "\n"))
    }
    best.exp[df] <- list(check.to.remove[[xx]]$exp)
    best.sig.indices[df] <- list(check.to.remove[[xx]]$sig.indices)
  }

  # Need to return the exposures in the context of the
  # original signatures matrix
  out.exp <- numeric(ncol(sigs)) #all zeros
  names(out.exp) <- colnames(sigs)
  max.df <- length(best.sig.indices)
  if (max.df == 0) {
    out.exp <- start.exp
  } else {
    out.exp[unlist(best.sig.indices[[max.df]])] <- unlist(best.exp[[max.df]])
  }
  stopifnot(abs(sum(out.exp) - sum(spect)) < 1)
  
  # Remove signature that have exposure less than 0.5
  sig.indices.non.zero.exposure <- which(out.exp >= 0.5)
  exposure <- as.matrix(out.exp[sig.indices.non.zero.exposure])
  colnames(exposure) <- colnames(spect)
  
  # Get the reconstructed spectrum
  recon <- ReconstructSpectrum(sigs = sigs, exp = exposure, use.sig.names = TRUE)
  
  distances <- 
    DistanceMeasuresSparse(spect = spect, recon = recon, 
                           nbinom.size = m.opts$nbinom.size,
                           likelihood.dist = m.opts$likelihood.dist,
                           signatures = sigs[, rownames(exposure), drop = FALSE])
  
  # Round the exposure and reconstruction
  exposure <- round(exposure)
  recon <- round(recon)
  
  # If there are signatures get assigned zero mutation counts after rounding,
  # remove these signatures from the exposure matrix
  non.zero.indices <- rowSums(exposure) > 0
  exposure <- exposure[non.zero.indices, , drop = FALSE]
  
  # Add attributes to recon to be same as spect
  recon <- CopyAttributes(to = recon, from = spect)
  
  return(list(proposed.assignment          = exposure,
              proposed.reconstruction      = recon,
              reconstruction.distances     = distances))
}
}

#' Calculate several distance measures between a spectrum and its reconstruction using
#' sparse assignment
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
#' @param likelihood.dist The probability distribution used to calculate the
#'   likelihood, can be either "multinom" (multinomial distribution) or
#'   "neg.binom" (negative binomial distribution).
#'   
#' @param signatures \strong{Only} used to compute distances for quadratic
#'   programming assignment. Signature as a matrix or data frame, with each row
#'   one mutation type (g.e. CCT > CAT or CC > TT) and each column a signature.
#'   It should be the proposed signatures used by sparse
#'   assignment to reconstruct \code{spect}. It is needed to calculate distances
#'   of quadratic programming assignment.
#'   
#' @return A data frame whose first column indicates the distance method. The
#'   second column \code{proposed.assignment} shows the values of various
#'   distances using sparse assignment. 
#'   
#'   When \code{signatures} is \strong{not} NULL, there will be a third column
#'   \code{QP.assignment} shows the values of various distances using quadratic
#'   programming assignment.
#'   
#' @keywords internal

DistanceMeasuresSparse <- 
  function(spect, recon, nbinom.size, likelihood.dist = "multinom", 
           signatures = NULL, m.opts = NULL) {
    my.fn <- function(method, spect, recon) {
      df <- rbind(as.vector(spect),
                  as.vector(recon))
      return(suppressMessages(philentropy::distance(x = df, 
                                                    method = method, 
                                                    test.na = FALSE)))
    }
    
    vv <- unlist(lapply(c("euclidean", "manhattan","cosine"), my.fn,
                        spect = spect, recon = recon))
    
    if (is.null(likelihood.dist)) {
      log.likelihood <- m.opts$loglh.fn(as.vector(spect), as.vector(recon))
    } else if (likelihood.dist == "multinom") {
      log.likelihood <- LLHSpectrumMultinom(as.vector(spect), as.vector(recon))
    } else if (likelihood.dist == "neg.binom") {
      log.likelihood <- 
        LLHSpectrumNegBinom(as.vector(spect), as.vector(recon), nbinom.size = nbinom.size)
    }
    
    vv <- c(log.likelihood = log.likelihood, vv)
    
    if (!is.null(signatures)) {
      # Do signature assignment using QP
      QP.expo <- OptimizeExposureQP(spectrum = spect, signatures = signatures)
      QP.expo.non.zero <- QP.expo[QP.expo >= 0.5]
      QP.recon <- ReconstructSpectrum(sigs = signatures, exp = QP.expo.non.zero,
                                      use.sig.names = TRUE)
      QP.distances <- unlist(lapply(c("euclidean", "manhattan","cosine"), my.fn,
                                    spect = spect, recon = QP.recon))
      
      if (is.null(likelihood.dist)) {
        log.likelihood <- m.opts$loglh.fn(as.vector(spect), as.vector(QP.recon))
      } else if (likelihood.dist == "multinom") {
        log.likelihood <- LLHSpectrumMultinom(as.vector(spect), as.vector(QP.recon))
      } else if (likelihood.dist == "neg.binom") {
        log.likelihood <- 
          LLHSpectrumNegBinom(as.vector(spect), as.vector(QP.recon), nbinom.size = nbinom.size)
      } 
      
      QP.distances <- c(log.likelihood = log.likelihood, QP.distances)
      
      return(tibble::tibble(method = names(vv), proposed.assignment = vv,
                            QP.assignment = QP.distances))
    } else {
      return(tibble::tibble(method = names(vv), proposed.assignment = vv))
    }
  }