DefaultGlobalOpts <- function() {
  return(
    list(algorithm    = "NLOPT_GN_DIRECT",
         xtol_rel      = 1e-9,
         # print_level = print_level,
         print_level   = 0,
         maxeval       = 10000))
}

DefaultLocalOpts <- function() {
  return(list(algorithm   = "NLOPT_LN_COBYLA",
              xtol_rel    = 1e-15,
              print_level = 0,
              maxeval     = 10000))
}

#' Set default options for many functions, especially \code{\link[nloptr]{nloptr}}.
#' 
#' @export
#' 
#' @return A list with the following elements
#' \describe{
#'   \item{global.opts}{Options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the global optimization phase.}
#'   \item{local.opts}{Options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the local optimization phase.}
#'   \item{nbinom.size}{The \code{size}
#'    parameter used by \code{\link[stats]{NegBinomial}}.}
#'   \item{trace}{If > 0 print progress messages.}
#' }
DefaultManyOpts <- function() {
  return(list(
    global.opts = DefaultGlobalOpts(),
    local.opts  = DefaultLocalOpts(),
    nbinom.size = 5,
    trace       = 0
  ))
}

#' Optimize assignment of a fixed set of signature activities for a \strong{single} tumor. 
#'
#' @param spectrum A single mutational spectrum a numeric vector or
#' single-column numeric matrix.
#' 
#' @param sigs A matrix of mutational signatures.
#' 
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#' 
#' @param m.opts See \code{\link{DefaultManyOpts}}.
#' 
#' @param local.opts See \code{\link[nloptr]{nloptr}}.
#' 
#' @param ... Additional arguments for \code{eval_f}.
#' 
#' @keywords internal
Nloptr1Tumor <- function(spectrum, 
                         sigs,
                         m.opts = NULL,
                         eval_f,
                         ... ) {
  if (!"matrix" %in% class(sigs)) {
    if (!"numeric" %in% class(sigs)) {
      stop("Unexpected class argument: ", class(sigs))
    }
    # Otherwise sigs is a numeric vector, not a matrix. We assume the caller
    # intended it as a single column matrix.
    sigs <- matrix(sigs, ncol = 1)
  }
  
  if (nrow(sigs) != length(spectrum)) {
    # If we get here there is an error. We save specturm and sigs, which seems
    # useful for debugging in a call to mclapply (multi-core lapply).
    save(spectrum, sigs, file = 'spectrum.sigs.debug.Rdata')
    stop("nrow(sigs) != length(spectrum), look in file spectrum.sigs.debug.Rdata")
  }
  
  stopifnot(mode(spectrum) == "numeric")
  foo <- spectrum
  storage.mode(spectrum) <- "double"
  stopifnot(foo == spectrum)
  
  if (is.null(m.opts)) m.opts <- DefaultManyOpts()
  
  # x0 is uniform distribution of mutations across signatures
  # (Each signature gets the same number of mutations)
  x0 <- rep(sum(spectrum) / ncol(sigs), ncol(sigs))
  
  if (!is.null(m.opts$global.opts)) { # WARNING, ADDITIONAL CODE WOULD NEED TO BE CHANGED DISABLE THIS BRANCH
    # set.seed(101010, kind = "L'Ecuyer-CMRG")
    global.res <- nloptr::nloptr(
      x0       = x0,
      eval_f   = eval_f,
      lb       = rep(0, ncol(sigs)),
      ub       = rep(sum(spectrum), ncol(sigs)), 
      opts     = m.opts$global.opts,
      spectrum = spectrum,
      sigs     = sigs,
      ...)
    if (m.opts$trace > 0) {
      message("globa.res$objective = ", global.res$objective)
    }
    if (global.res$iterations == m.opts$global.opts[["maxeval"]]) {
      # warning("reached maxeval on global optimization: ", global.res$iterations)
    }
  }

  local.res <- nloptr::nloptr(
    # x0=x0,
    x0       = global.res$solution,
    eval_f   = eval_f,
    lb       = rep(0, ncol(sigs)),
    ub       = rep(sum(spectrum) + 1e-2, ncol(sigs)), 
    opts     = m.opts$local.opts, 
    spectrum = spectrum,
    sigs     = sigs,
    ...)
  if (m.opts$trace > 0)
    message("local.res$objective = ", local.res$objective)
  
  if (local.res$iterations == m.opts$local.opts[["maxeval"]]) {
    warning("reached maxeval on local optimization: ", local.res$iterations)
  }
  

  names(local.res$solution) <- colnames(sigs)
  return(list(objective =  local.res$objective,
              solution   = local.res$solution,
              global.res = global.res,
              local.res  = local.res))
}

#' Returns a list with elements:
#'
#' \code{loglh}    - the log likelihood of the best solution (set of exposures) found
#' \code{exposure}  - the vector of exposures that generate \code{loglh}, in this case
#' this is the number of mutations ascribed to each signature.
#'
#' @keywords internal
one.lh.and.exp <- function(spect,
                           sigs, 
                           m.opts,
                           eval_f) {
  
  stopifnot(mode(spect) == "numeric")

  r <- Nloptr1Tumor(spect, 
                    sigs, 
                    m.opts = m.opts,
                    eval_f      = eval_f,
                    nbinom.size = m.opts$nbinom.size)
  loglh <- r$objective
  
  if (loglh == Inf && m.opts$trace > 0)
    message("Got -Inf in one.lh.and.exp\n")
  
  # sum(recon) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spect) / sum(r$solution)) # sum(recon))
  
  list(loglh = -loglh, exposure = exp, everything.else = r)
}


#' Is the set 'probe' a superset of any set in 'background'?
#' 
#' @keywords internal
is.superset.of.any <- function(probe, background) {
  for (b in background) {
    if (sets::set_is_proper_subset(b, probe)) return(TRUE)
  }
  return(FALSE)
}

#' Find known signatures that can most sparsely reconstruct each spectrum in a catalog.
#' 
#' @param spectra The spectra (multiple spectra) to be reconstructed.
#' 
#' @param sigs The known signatures to use in reconstruction.
#' 
#' @param max.level The largest number of signatures to consider discarding
#' in the reconstruction.
#' 
#' @param p.thresh The maximum p value based on which it is decided
#' to retain a signature in a reconstruction.
#' 
#' @param eval_f The objective function for
#'  \code{\link[nloptr]{nloptr}}. If \code{NULL} defaults
#'  to \code{\link{ObjFnBinomMaxLH}}.
#'  
#' @param m.opts If \code{NULL},
#'   defaults to \code{\link{DefaultManyOpts}()}.
#' 
#' @param mc.cores The number of cores to use; if \code{NULL} use
#'  \code{min(10, ncol(spectra))}, except on Windows, where \code{mc.cores}
#'  is always 1.
#' 
#' @return A list with the inferred exposure matrix as element \code{exposure}.
#' 
#' @export

SparseAssignActivity <- function(spectra,
                                 sigs,
                                 max.level = 5,
                                 p.thresh  = 0.5,
                                 eval_f    = NULL,
                                 m.opts    = NULL,
                                 mc.cores  = NULL) {
  f1 <- function(i) {
    retval1 <- SparseAssignActivity1(
      spect    = spectra[ , i, drop = FALSE],
      sigs     = sigs,
      p.thresh = p.thresh,
      m.opts   = m.opts,
      eval_f   = eval_f)
    
    return(retval1)
  }
  
  if (is.null(m.opts)) m.opts <- DefaultManyOpts()
  
  if (is.null(eval_f)) eval_f <- ObjFnBinomMaxLHNoRoundOK
  
  mc.cores <- Adj.mc.cores(
    ifelse(is.null(mc.cores), 
           min(10, ncol(spectra)),
           mc.cores))
  
  retval <- parallel::mclapply(1:ncol(spectra), f1, mc.cores = mc.cores)

  
  retval2 <- matrix(unlist(retval), ncol = length(retval))
  colnames(retval2) <- colnames(spectra)
  rownames(retval2) <- colnames(sigs)
  
  return(list(exposure = retval2))
  # return(retval)
  
}

#' Component of \code{\link{SparseAssignActivity}} for one spectrum.
#' @keywords internal

SparseAssignActivity1 <- function(
  spect, sigs, max.level = 5, p.thresh = 0.05, eval_f, m.opts) {
  mode(spect) <-  'numeric'
  start <- one.lh.and.exp(spect, 
                          sigs,
                          eval_f     = eval_f,
                          m.opts  = m.opts)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp > 0.5)
  if (m.opts$trace > 0) {
    message('Starting with',
        paste(names(start.exp)[non.0.exp.index], collapse = ","),
        '\n')
    message('max.level =', max.level, '\n')
  }
  if (length(non.0.exp.index) < 2) {
    if (m.opts$trace > 0)
      message("returning, only ", length(non.0.exp.index),
                       " non-0 exposures")
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)
  
  c.r.s <- sets::set() # subsets of the signature indices that cannot be removed
  
  max.p <- numeric(max.level)
  best.subset <- non.0.exp.index
  best.try <- start
  best.try.is.start <- T
  
  for (df in 1:max.level) {
    if (m.opts$trace > 1) message('df =', df, '\n')
    subsets <- sets::set_combn(non.0.exp.index, df)
    for (subset in subsets) {
      # subset is of class set (package sets)
      if (is.superset.of.any(subset, c.r.s)) next
      subset.to.remove.v <- as.numeric(subset)
      subset.name <- paste(names(start.exp)[subset.to.remove.v], collapse = ",")
      tmp.set <- setdiff(non.0.exp.index, subset.to.remove.v)
      try.sigs <- sigs[ , tmp.set]
      
      if (length(tmp.set) == 1) {
        try.sigs <- as.matrix(try.sigs)
        rownames(try.sigs) <- rownames(sigs)
        if (m.opts$trace > 0) message("Running new code\n")
      }
      
      # Get the max lh for try.sig / Get the maximum likelihood 
      # reconstruction using try.sig
      try <- one.lh.and.exp(spect, 
                            try.sigs,
                            eval_f    = eval_f,
                            m.opts = m.opts)
      # try contains maxlh and exposure
      statistic <- 2 * (lh.w.all - try$loglh)
      chisq.p <- stats::pchisq(statistic, df, lower.tail = FALSE)
      if (m.opts$trace > 0) {
        message("Trying" , subset.name, ", p = ", chisq.p)
        if (m.opts$trace > 1) 
          message("loglh = ", try$loglh, "; statistic  =", statistic)
      }
      if (chisq.p > p.thresh) {
        # This an acceptable set of exposures
        if (m.opts$trace > 0) message("acceptable")
        if (chisq.p > max.p[df]) {
          # This is the best so far
          max.p[df] <- chisq.p
          best.subset <- tmp.set
          best.try <- try
          best.try.is.start <- F
          if (m.opts$trace > 0) message('best\n')
        } else {
          if (m.opts$trace > 0) message('not best\n')
        }
      } else {
        c.r.s <- sets::set_union(c.r.s, sets::set(subset))
        if (m.opts$trace > 0)  {
          message('cannot remove\n')
        }
      }
    } # end for (subset in subsets)
    if (max.p[df] == 0) break
  } # end for df in 1:max.level
  
  # Need to return the exposures in the context of the 
  # orginal signatures matrix
  out.exp <- numeric(ncol(sigs)) #all zeros
  names(out.exp) <- colnames(sigs)
  if (best.try.is.start) {
    out.exp <- start.exp
  } else {
    out.exp[best.subset] <- best.try$exposure
  }
  if (m.opts$trace > 0) {
    message('max.p =', paste(max.p, collapse = ', '), '\n')
    message('Ending with',
        paste(names(start.exp)[best.subset], collapse = ","),
        '\n')
  }
  stopifnot(abs(sum(out.exp) - sum(spect)) < 1)
  return(out.exp)
}

# Helper function, given signatures (sigs) and exposures (exp), return a
# *proportional* reconstruction; in general, it is *not necessarily* scaled to
# the actual spectrum counts.
#' error checked (?) function to get reconstructed something?
#' 
#' @keywords internal
prop.reconstruct <- function(sigs, exp) {
  stopifnot(length(exp) == ncol(sigs))
  return(as.matrix(sigs) %*% exp)
}

#' The negative binomial maximum likelihood objective function.
#' 
#' For use by \code{\link[nloptr]{nloptr}}.
#'
#' @param exp The matrix of exposures ("activities").
#' @param spectrum The spectrum to assess.
#' @param sigs The matrix of signatures.
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'
#' The result is
#' 
#' -1 * log(likelihood(spectrum | reconstruction))
#'
#' (nloptr minimizes the objective function.)
#'
#' The lower the objective function, the better.
#'
#' Can be used as the
#' objective function for \code{\link{SparseAssignActivity}},
#' \code{\link{SparseAssignActivity1}}, 
#' and \code{\link{SignaturePresenceTest1}}.
#' 
#' @export
#' 
ObjFnBinomMaxLHMustRound <- function(exp, spectrum, sigs, nbinom.size) {
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE)
}

#' The negative binomial maximum likelihood objective function allowing no rounding.
#' 
#' For use by \code{\link[nloptr]{nloptr}}.
#' For details see \code{\link{ObjFnBinomMaxLHMustRound}}.
#'
#' @inheritParams ObjFnBinomMaxLHMustRound
#'
#' @export
#' 
ObjFnBinomMaxLHNoRoundOK <- function(exp, spectrum, sigs, nbinom.size) {
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = TRUE)
}

#' The negative binomial maximum likelihood objective function.
#' 
#' For use by \code{\link[nloptr]{nloptr}}
#'
#' @param exp The matrix of exposures ("activities").
#' @param spectrum The spectrum to assess.
#' @param sigs The matrix of signatures.
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#' @param no.round.ok If \code{TRUE}, allow use of unrounded
#'        reconstruction if some mutation types would have 0
#'        counts in the reconstructed spectrum.
#'
#' The result is
#' 
#' -1 * log(likelihood(spectrum | reconstruction))
#'
#' (nloptr minimizes the objective function.)
#'
#' The lower the objective function, the better.
#'
#' Can be used as the
#' objective function for \code{\link{SparseAssignActivity}},
#' \code{\link{SparseAssignActivity1}}, 
#' and \code{\link{SignaturePresenceTest1}}.
#' 
#' @export
#' 
ObjFnBinomMaxLH2 <- 
  function(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE,
           show.warning = FALSE) {
  
  if (any(is.na(exp))) return(Inf)
  
  reconstruction <-  prop.reconstruct(sigs = sigs, exp = exp)
  
  reconstruction2 <- round(reconstruction)
  # Will cause problems if round of the reconstruction is 0 for
  # any channel even if the reconstruction > 0, because then
  # log likelihood will be -inf. The situation is especial likely
  # to occur if mutation counts in the spectrum are low.
  if (any(reconstruction2 == 0) && no.round.ok) {
    if (show.warning) warning("unrounded reconstruction")
  } else {
    reconstruction <- reconstruction2
  }
  rm(reconstruction2)
  
  ## catch errors with NA in the input or in reconstruction.
  if (any(is.na(reconstruction))) {
    save(reconstruction, spectrum, sigs, exp, nbinom.size,
         file = "reconstruction.error.R")
  }
  stopifnot(!any(is.na(reconstruction)))
  
  loglh <- LLHSpectrumNegBinom(spectrum = spectrum, 
                                expected.counts = reconstruction,
                                nbinom.size = nbinom.size)

  return(-loglh)
}


#' Euclidean reconstruction error.
#' 
#' @keywords internal
EDist2SpectRounded <- function(exp, sig.names, spect) {
  reconstruction <- 
    prop.reconstruct(
      sigs = PCAWG7::signature$genome$SBS96[ , sig.names], 
      exp = round(exp))
  # TEST
  reconstruction <- round(reconstruction)
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


#' Euclidean reconstruction error.
#' 
#' @keywords internal
EDist2Spect <- function(exp, sig.names, spect) {
  reconstruction <- 
    prop.reconstruct(
      sigs = 
        PCAWG7::signature$genome$SBS96[ , sig.names], exp = exp)
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


#' "Polish" a solution by minimizing Euclidean distance to the input spectrum.
#' 
#' @keywords internal
Polish <- function(exp, sig.names, spect) {
  retval <- nloptr::nloptr(x0 = exp,
                            eval_f = EDist2Spect,
                            lb = rep(0, length(exp)),
                            ub = rep(sum(exp), length(exp)),
                            opts = list(algorithm   = "NLOPT_LN_COBYLA",
                                        maxeval     = 1000, 
                                        print_level = 0,
                                        xtol_rel    = 0.001,
                                        xtol_abs    = 0.0001),
                            sig.names = sig.names,
                            spect     = spect)
    
  names(retval$solution) <- names(exp)  
  return(retval$solution)
  
}


#' Framework for testing \code{\link{SparseAssignActivity1}}.
#' 
#' Tests only for single spectrum.
#' 
#' @keywords internal

SparseAssignTest1 <- function(sig.counts, trace = 0) {
  
  # set.seed(1449, kind = "L'Ecuyer-CMRG")
  
  sig.names <- names(sig.counts)
  
  some.sigs  <- 
    PCAWG7::signature$genome$SBS96[ , sig.names, drop = FALSE]
  ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
  region     <- attr(some.sigs, "region", exact = TRUE)
  if (is.null(region)) {
    message("Null region, why?")
    region <- "genome"
  }
  
  spect <- round(some.sigs %*% sig.counts)
  spect <-
    ICAMS::as.catalog(
      spect, 
      ref.genome   = ref.genome,
      region       = region,
      catalog.type = "counts")
  nbinom.size <- 5

  m.opts <- DefaultManyOpts()
  m.opts$trace <- trace
  
  SA.out <- SparseAssignActivity1(spect      = spect,
                                 sigs        = some.sigs,
                                 eval_f      = ObjFnBinomMaxLHMustRound,
                                 m.opts      = m.opts)
  
  zeros <- which(SA.out < 0.5)
  if (length(zeros) > 0) {
    SA.out2 <- SA.out[-zeros]
    sig.names2 <- sig.names[-zeros]
  } else {
    SA.out2 <- SA.out
    sig.names2 <- sig.names
  }
  
  polish.out <- Polish(exp    = SA.out2,
                       sig.names = sig.names2,
                       spect     = spect)
  
  recon1 <- round(prop.reconstruct(some.sigs, SA.out))
  recon2 <-
    round(
      prop.reconstruct(
        PCAWG7::signature$genome$SBS96[ , sig.names2, drop = FALSE],
        polish.out))
  
  
  return(list(soln1       = SA.out,
              soln2       = polish.out,
              truth       = sig.counts,
              edist1      = EDist2Spect(SA.out, sig.names, spect),
              edist1r     = EDist2SpectRounded(SA.out, sig.names, spect),
              LL1         = LLHSpectrumNegBinom(spect, recon1, nbinom.size),
              LL2         = LLHSpectrumNegBinom(spect, recon2, nbinom.size),
              edist2      = EDist2Spect(polish.out, sig.names2, spect),
              ed8st2r     = EDist2SpectRounded(polish.out, sig.names2, spect)
              #, input.spect = spect
              ))
  
}


#' Framework for testing \code{\link{SparseAssignActivity}}.
#' 
#' @param sig.counts A matrix of target exposures:
#' rows are signatures and columns are tumors.
#' 
#' @param trace If > 0 print progress messages.
#' 
#' @keywords internal

SparseAssignTest <- function(sig.counts, trace = 0, mc.cores = NULL) {
  
  if (!is.matrix(sig.counts)) {
    tmp.names <- names(sig.counts)
    sig.counts <- matrix(sig.counts, ncol = 1)
    rownames(sig.counts) <- tmp.names
  }
  sig.names <- rownames(sig.counts)
  colnames(sig.counts) <- paste0("tumor", 1:ncol(sig.counts))
  
  some.sigs  <- 
    PCAWG7::signature$genome$SBS96[ , sig.names, drop = FALSE]

  ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
  region     <- attr(some.sigs, "region", exact = TRUE)
  if (is.null(region)) {
    message("Null region, why?")
    region <- "genome"
  }
  
  spect <- round(some.sigs %*% sig.counts)
  spect <-
    ICAMS::as.catalog(
      spect, 
      ref.genome   = ref.genome,
      region       = region,
      catalog.type = "counts")
  
  m.opts <- DefaultManyOpts()
  m.opts$trace <- trace
  
  SA.out <- SparseAssignActivity(spectra    = spect,
                                  sigs      = some.sigs,
                                  eval_f    = ObjFnBinomMaxLHMustRound,
                                  m.opts    = m.opts,
                                  mc.cores  = mc.cores) 
  return(SA.out)
}


#' Test whether a given signature is plausibly present in a spectrum.
#' 
#' @param spectrum The spectrum to analyze
#' 
#' @param sigs A catalog of signatures from which to choose
#' 
#' @param target.sig.index The index of the signature the presence
#' of which we want to test.
#' 
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#' 
#' @param m.opts For documentation
#'    see \code{\link{DefaultManyOpts}}.
#' 
#' @export

SignaturePresenceTest1 <- function(
  spectrum, sigs, target.sig.index, m.opts, eval_f) {
  
  ret.with <- one.lh.and.exp(spect  = spectrum,
                             sigs   = sigs, 
                             m.opts = m.opts,
                             eval_f = eval_f)
  loglh.with <- ret.with$loglh

  ret.without <- one.lh.and.exp(spect  = spectrum, 
                                sigs   = sigs[ ,-target.sig.index],
                                eval_f = eval_f,
                                m.opts = m.opts)
  loglh.without <- ret.without$loglh
  
  statistic <- 2 * (loglh.with - loglh.without)
  chisq.p <- stats::pchisq(statistic, 1, lower.tail = FALSE)
  if (m.opts$trace > 0) 
    message("statistic  = ", statistic, "\nchisq p = ", chisq.p)
  
  list(with                    = loglh.with,
       without                 = loglh.without,
       statistic               = statistic,
       chisq.p                 = chisq.p,
       exp.with                = ret.with$exp,
       exp.without             = ret.without$exp,
       everything.else.with    = ret.with$everything.else,
       everything.else.without = ret.without$everything.else)
}


Adj.mc.cores <- function(mc.cores) {
  if (Sys.info()["sysname"] == "Windows" && mc.cores > 1) {
    message("On Windows, changing mc.cores from ", mc.cores, " to 1")
    return(1)
  }
  return(mc.cores)
}

#' Test whether a given signature is plausibly present in a catalog
#' 
#' @param spectra The catalog (matrix) to analyze. This could be
#'   an \code{\link[ICAMS]{ICAMS}} catalog or a numerical matrix.
#' 
#' @param sigs A catalog of signatures from which to choose.
#' This could be
#'   and \code{\link[ICAMS]{ICAMS}} catalog or a numerical matrix.
#' 
#' @param target.sig.index The index of the signature the presence
#' of which we want to test.
#' 
#' @param eval_f See \code{\link[nloptr]{nloptr}}.
#' 
#' @param m.opts If \code{NULL} use the return from 
#'    calling \code{\link{DefaultManyOpts}}. For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'    
#' @param mc.cores Number of cores to use. Always silently 
#'  changed to 1 on Microsoft Windows.
#' 
#' @export

SignaturePresenceTest <- function(
  spectra, sigs, target.sig.index, m.opts = NULL, eval_f, mc.cores = 10) {
  
  # check if signatures sum to 1
  all.col.sums <- colSums(sigs)
  
  if (!all.equal(all.col.sums,
                 rep(1, ncol(sigs)), 
                 tolerance = 1e-3,
                 check.names = FALSE)) {
    stop(
      "Argument sigs does not seem to be a signature ",
      "catalog because some columns do not sum to 1")
  }
  
  # Need to match exactly one signature name
  stopifnot(length(target.sig.index) == 1)
  
  spectra.as.list <- split(t(spectra), 1:ncol(spectra))
  
  names(spectra.as.list) <- colnames(spectra)
  
  if (is.null(m.opts)) m.opts <- DefaultManyOpts()

  mc.cores <- Adj.mc.cores(mc.cores)
  
  out.res <-
    parallel::mclapply(
      X                = spectra.as.list,
      FUN              = SignaturePresenceTest1,
      mc.cores         = mc.cores,
      sigs             = sigs,
      target.sig.index = target.sig.index,
      m.opts           = m.opts,
      eval_f           = eval_f)
  
  return(out.res)
}


#' Framework for testing \code{\link{SignaturePresenceTest1}}.
#' 
#' @keywords internal
#'
TestSignaturePresenceTest1 <- 
  function(sig.counts, 
           input.sigs = PCAWG7::signature$genome$SBS96, 
           trace      = 0,
           eval_f     = ObjFnBinomMaxLHNoRoundOK,
           m.opts     = NULL) {
  
    if (!requireNamespace("BSgenome.Hsapiens.1000genomes.hs37d5",
                          quietly = TRUE)) {
      stop("Please install Bioconductor library ",
           "BSgenome.Hsapiens.1000genomes.hs37d5")
    }
    
    sig.names <- names(sig.counts)
  
  if (sum(sig.names %in% colnames(input.sigs)) == 0) {
    stop("sig.names not all in input.sigs")
  }
  
  some.sigs  <- input.sigs[ , sig.names, drop = FALSE]
  
  ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
  region     <- attr(some.sigs, "region", exact = TRUE)
  if (is.null(region)) {
    message("Null region, why?")
    region <- "genome"
  }
  
  spectrum <- round(some.sigs %*% sig.counts)
  spectrum <-
    ICAMS::as.catalog(
      object       = spectrum, 
      ref.genome   = ref.genome,
      region       = region,
      catalog.type = "counts")
  
  if (is.null(m.opts)) {
    m.opts <- DefaultManyOpts()
  }
  m.opts$trace <- trace
  
  retval <- SignaturePresenceTest1(
    spectrum         = spectrum,
    sigs             = some.sigs, 
    target.sig.index = 1,
    eval_f           = eval_f,
    m.opts           = m.opts)
  
  return(retval)
}


# https://cran.r-project.org/web/packages/DirichletReg/DirichletReg.pdf


#' Deterimine whether any of several signatures in any combination are plausibly needed to reconsruct a given spectrum.
#' 
#' @param spect The spectrum to be reconstructed, as single column matrix or
#'  \code{\link[ICAMS]{ICAMS}} catalog.
#'  
#' @param all.sigs The matrix or catalog of singantures of possible interest;
#' \code{all.sigs} includes the signatures that might or might not be necessary
#' for plausibly reconstructing \code{spect}; \code{all.sigs} includes the
#' signatures in \eqn{H_0} plus the additional signatures in \eqn{H_a}.
#' 
#' @param target.sigs.index An integer vector of the indices of the signatures
#' of interest (those in \eqn{H_a - H_0}.
#' 
#' @param eval_f XXXXX
#' 
#' @param m.opts XXXX
#' 
#' Return XXXXXX 
#' 
#' 
# the P value of the null hypothesis that none of the signatures in Ha.sigs,
# either singly or in combination **** probabily of a likelihood as good or better
# generating the spectrum from the H0 signatures **** as opposed to compared to incluing the
# extra signatures.
#
#' WARNING: tests all non-empty subsets of \code{target.sigs}, so will get very slow for
#' large numbers of \code{target.sigs}.


# START HERE
# To Do: 1. remove empty set from more.sigs' subsets
#        3. see if it is reasonable to require all signatures in more.sigs
# @param spect A numerical matrix an \code{\link[ICAMS]{ICAMS}} catalog containg 
# mutational spectra.
AnySigSubsetPresent <- 
  function(spect, 
           all.sigs,          
           target.sigs.index, 
           # p.thresh = 0.05, 
           eval_f = mSigAct::ObjFnBinomMaxLHNoRoundOK,  
           m.opts) {
    
    H0.sigs <- all.sigs[ , -target.sigs.index, drop = FALSE]  # H0.sigs a matrix of sigs

    start <- one.lh.and.exp(spect  = spect, 
                            sigs   = H0.sigs,
                            eval_f = eval_f,
                            m.opts = m.opts)
    
    
    H0.lh <- start$loglh  # The likelihood with only H0.sigs
    start.exp <- start$exposure
    zero.exposures <- which(start.exp < 0.5)
    if (length(zero.exposures) > 0) {
      message("There were some signatures with no expsoures in H0\n")
    }
    
    # For information, in case some signatures are useless
    non.0.exp <- names(start.exp)[which(start.exp > 0.5)]

    if (m.opts$trace > 0) {
      message('H0 sigs are', paste(colnames(H0.sigs), collapse = ','),'\n')
    }
    
    # new.subsets contains all non-empty subsets of more.sigs (2^as.set(c()))
    # is the powerset of the empty set, i.e. {{}}).
    new.subsets <- 2^sets::as.set(target.sigs.index) - 2^sets::as.set(c())
    
    inner.fn <- function(sigs.to.add) {
      df <- sets::set_cardinality(sigs.to.add) # degrees of freedom for likelihood ratio test
      Ha.info <- one.lh.and.exp(spect  = spect,
                                sigs   = cbind(all.sigs[ , unlist(sigs.to.add), drop = FALSE ], H0.sigs),
                                eval_f = eval_f,
                                m.opts = m.opts)
      statistic <- 2 * (Ha.info$loglh - H0.lh)
      p <- pchisq(statistic, df, lower.tail = FALSE)
      return(list(sigs.added = paste(colnames(all.sigs)[unlist(sigs.to.add)], collapse = ","),
                  statistic  = statistic,
                  df         = df,
                  p          = p,
                  base.sigs  = paste(non.0.exp, collapse = ","),
                  Ha.info    = Ha.info))
    }
    
    out <- lapply(new.subsets, inner.fn)
    return(list(H0.info = start, all.Ha.info = out))
  }


TestEsoSpectra <- function(indices = NULL) {
  eso.index <- grep("Eso", colnames(PCAWG7::PCAWG.WGS.SBS.96), fixed = TRUE)
  spectra <- PCAWG7::PCAWG.WGS.SBS.96[ , eso.index, drop = FALSE]
  if (!is.null(indices)) {
    spectra <- spectra[ , indices, drop = FALSE]
  }
  return(spectra)
}


TestAny1 <- function(extra.sig, eso.index) {
  
  eso.spectra <- TestEsoSpectra(eso.index)
  
  m.opts <- DefaultManyOpts()

  sigs.plus <- TestEsoSigs(extra.sig) # The extra signatures are signature names, and will be the first columns of sigs.plus

  set.seed(101010, kind = "L'Ecuyer-CMRG")  
  out <- AnySigSubsetPresent(spect             = eso.spectra,
                             all.sigs          = sigs.plus,
                             target.sigs.index = 1:length(extra.sig),
                             eval_f            = mSigAct::ObjFnBinomMaxLHNoRoundOK,
                             m.opts            = m.opts)
  
  return(out)
}


if (FALSE) {
  any.retval <- TestAny1("SBS17a", 1)
}


TestEsoSigs <- function(extra.sigs = NULL) {
  sigs <- c(
    "SBS1",
    "SBS2",
    "SBS3",
    "SBS5",
    "SBS13",
    "SBS18",
    "SBS28",
    "SBS40")
  if (!is.null(extra.sigs)) {
    sigs <- c(extra.sigs, sigs)
  }
  return(PCAWG7::signature$genome$SBS96[ , sigs])
}

TestSignaturePresenceTest <- function(extra.sig, eso.indices) {

  eso.spectra <- TestEsoSpectra(eso.indices)
  
  m.opts <- DefaultManyOpts()
  
  sigs.plus <- TestEsoSigs(extra.sig)
  set.seed(101010, kind = "L'Ecuyer-CMRG") 
  retval1 <- mSigAct::SignaturePresenceTest(
    spectra          = eso.spectra, 
    sigs             = sigs.plus,
    target.sig.index = 1, 
    m.opts           = m.opts, 
    eval_f           = ObjFnBinomMaxLHNoRoundOK, 
    mc.cores         = 1)
  
  set.seed(101010, kind = "L'Ecuyer-CMRG")
  retval2 <- SignaturePresenceTest1(
    spectrum         = eso.spectra,
    sigs             = sigs.plus,
    target.sig.index = 1,
    m.opts           = m.opts,
    eval_f           = ObjFnBinomMaxLHNoRoundOK)
  
  stopifnot(all.equal(eso.17a[[1]][[1]], eso.17a[[2]]))
  
  return(list(retval1, retval2))
}


if (FALSE) {
  spt.retval <- TestSignaturePresenceTest("SBS17a", 1)
  stopifnot(all.equal(eso.17a$`Eso-AdenoCA::SP111062`$chisq.p, 0.1019716, tolerance = 1e-7))

  }

