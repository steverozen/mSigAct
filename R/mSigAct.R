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
#'   
#'   \item{local.opts}{Options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the local optimization phase.}
#'   
#'   \item{nbinom.size}{The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.}
#'    
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
#' @inheritParams OptimizeExposure 
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

#' Optimize the reconstruction of a spectrum from a set of signatures.
#' 
#' @param spectrum The spectrum to be reconsructed.
#' 
#' @param sigs The available signatures.
#' 
#' @param m.opts Options that govern the numerical optimizaiton. 
#'    For documentation see \code{\link{DefaultManyOpts}}.
#' 
#' @param eval_f The objective function for
#'  \code{\link[nloptr]{nloptr}}. We have only tested
#'  \code{\link{ObjFnBinomMaxLHNoRoundOK}} and
#'  \code{\link{ObjFnBinomMaxLHMustRound}}.
#'  
#' @param ... Additional arguments for \code{eval_f}.
#'
#' Returns a list with elements \describe{
#' \item{\code{loglh}}{The log likelihood of the best solution (set of exposures) found.
#'       For a more general objective function this might be \code{NULL}.}
#' \item{\code{exposure}}{The vector of exposures that generate \code{loglh}, i.e. 
#'    the number of mutations ascribed to each signature.}
#' \item{\code{obj.fn.value}}{The objective function value associated with \code{exposure}.}
#' \item{\code{everything.else}}{Everything returned by the function \code{\link{Nloptr1Tumor}.} }
#' }
#' 
#' @export
#' 
OptimizeExposure <- function(spectrum,
                             sigs, 
                             m.opts,
                             eval_f,
                             ...) {
  
  stopifnot(mode(spectrum) == "numeric")
  
  r <- Nloptr1Tumor(spectrum = spectrum, 
                    sigs, 
                    m.opts = m.opts,
                    eval_f      = eval_f,
                    nbinom.size = m.opts$nbinom.size,
                    ...)
  loglh <- r$objective
  
  if (loglh == Inf && m.opts$trace > 0)
    message("Got -Inf in one.lh.and.exp\n")
  
  # sum(r$solution) is likely to be close to, but not equal to the number
  # of mutations in the spectrum, so we scale exposures to get the
  # number of mutations in the spectrum
  exp <- r$solution * (sum(spectrum) / sum(r$solution)) # sum(recon))
  
  # TODO, there is redudant info in everything.else
  return(list(loglh = -loglh, exposure = exp, obj.fn.value = r$objective, everything.else = r))
}

#' @keywords  internal
one.lh.and.exp <- function(...) OptimizeExposure(...) 

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
#'  \code{\link[nloptr]{nloptr}}.
#'  
#' @param m.opts For documentation
#'    see \code{\link{DefaultManyOpts}}.
#'    
#' @param num.parallel.samples The (maximum) number of samples to run in parallel; each
#'    sample in turn can require multiple cores, as governed by
#'    \code{mc.cores.per.sample}.
#' 
#' @param mc.cores.per.sample 
#'    The maximum number of cores to use for each sample.
#'    If \code{NULL} defaults to \code{2^max.level} -- except on
#'    MS Windows machines, where it defaults to 1.
#' 
#' @return A list with the inferred exposure matrix as element \code{exposure}.
#' 
#' @export

SparseAssignActivity <- 
  function(spectra,
           sigs,
           max.level            = 5,
           p.thresh             = 0.05,
           eval_f               = ObjFnBinomMaxLHNoRoundOK,
           m.opts               = NULL,
           num.parallel.samples = 5,
           mc.cores.per.sample  = NULL) {
    f1 <- function(i) {
      retval1 <- SparseAssignActivity1(
      spect        = spectra[ , i, drop = FALSE],
      sigs         = sigs,
      p.thresh     = p.thresh,
      m.opts       = m.opts,
      eval_f       = eval_f,
      max.level    = max.level,
      max.mc.cores = mc.cores.per.sample)
    
    return(retval1)
  }
  
  if (is.null(m.opts)) m.opts <- DefaultManyOpts()
  
  num.parallel.samples <- Adj.mc.cores(num.parallel.samples)
  
  retval <- parallel::mclapply(1:ncol(spectra),
                               f1, 
                               mc.cores = num.parallel.samples)

  retval2 <- matrix(unlist(retval), ncol = length(retval))
  colnames(retval2) <- colnames(spectra)
  rownames(retval2) <- colnames(sigs)
  
  return(list(exposure = retval2))
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

#' A deprecated negative binomial maximum likelihood objective function.
#' 
#' Use \code{\link{ObjFnBinomMaxLHNoRoundOK}} instead.
#' 
#' This function will lead to errors in some situations
#' when the rounded reconstructed signature contains 0s for
#' mutations classes for which the target spectrum is > 0.
#'
#' @inheritParams ObjFnBinomMaxLHNoRoundOK
#'
#' @export
#' 
ObjFnBinomMaxLHMustRound <- function(exp, spectrum, sigs, nbinom.size) {
  ObjFnBinomMaxLH2(exp, spectrum, sigs, nbinom.size, no.round.ok = FALSE)
}

#' The preferred negative binomial maximum likelihood objective function.
#' 
#' Can be used as the
#' objective function for \code{\link{SparseAssignActivity}},
#' \code{\link{SparseAssignActivity1}}, 
#' and \code{\link{SignaturePresenceTest1}}.
#' (Internally used by by \code{\link[nloptr]{nloptr}}.)
#' 
#' @return 
#'
#' -1 * log(likelihood(spectrum | reconstruction))
#'
#' \code{\link[nloptr]{nloptr}} minimizes the objective function, so the 
#' lower the objective function, the better.
#'
#' @param exp The matrix of exposures ("activities").
#' @param spectrum The spectrum to assess.
#' @param sigs The matrix of signatures.
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.
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
# @param exp The matrix of exposures ("activities").
# @param spectrum The spectrum to assess.
# @param sigs The matrix of signatures.
# @param nbinom.size The dispersion parameter for the negative
#        binomial distribution; smaller is more dispersed.
#        See \code{\link[stats]{NegBinomial}}.
#'        
#' @inheritParams ObjFnBinomMaxLHNoRoundOK
#'
#' @param no.round.ok If \code{TRUE}, allow use of unrounded
#'        reconstruction if some mutation types would have 0
#'        counts in the reconstructed spectrum. Deprecated,
#'        in future will always be \code{TRUE}.
#'        
#' @param show.warning Deprecated; ignored.
#'
#' @keywords internal
#' 
NEW.VERSION.ObjFnBinomMaxLH2 <- 
  function(exp, spectrum, sigs, nbinom.size, no.round.ok = TRUE,
           show.warning = FALSE) {
    
    if (any(is.na(exp))) return(Inf)
    
    reconstruction.estimate <-  prop.reconstruct(sigs = sigs, exp = exp)
    
    # The reconstruction.estimate is an abstract model of the probabilities of 
    # the counts of each mutation class given the mixture of exposures in exp
    # (which are also probabilities) and the signatures. The probabilities of
    # the counts do not have to be integers.

    # catch errors with NA in the input or in reconstruction.
    if (any(is.na(reconstruction.estimate))) {
      save(reconstruction.estimate, spectrum, sigs, exp, nbinom.size,
           file = "reconstruction.NA.Rdata")
      stop("There arevNAs in reconstruction; see reconstruction.NA.Rdata")
    }
    
    if (!no.round.ok) {
      # There can be problems if the rounded reconstruction is 0 for
      # any channel (even if the unrounded reconstruction > 0), because then
      # log likelihood will be -inf. The situation is especial likely
      # to occur if mutation counts in the spectrum are low.
      if (any(round(reconstruction.estimate)[spectrum > 0] == 0)) {
        # warning("Possible problem in rounding reconstucted spectrum; ",
        #        "use eval_f = ObjFnBinomMaxLHNoRoundOK")
      } 
      # Use the rounded reconstruction estimate
      reconstruction.estimate <- round(reconstruction.estimate)      
    }
    
    loglh <- mSigBG::LLHSpectrumNegBinom(spectrum = spectrum, 
                                         expected.counts = reconstruction.estimate,
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
  reconstruction <- round(reconstruction)
  class(spect) <- "matrix"
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
  class(spect) <- "matrix"
  err <- stats::dist(t(cbind(reconstruction, spect)), method = "euclidean")
  return(err)
}


#' "Polish" a solution by minimizing Euclidean distance to the input spectrum.
#' 
#' This is experimental for testing.
#' 
#' @keywords internal
Polish <- function(exp, sig.names, spect) {
  class(spect) <- "matrix" # Otherwise cbind will use the ICAMS catalog
                           # method, which may complain that the reconstructed
                           # output is not a catalog.
  
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
#' @keywords internal

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


#' For each combination of several signatures, determine if the combination is plausibly needed to reconstruct a spectrum.
#' 
#' @description Please see \strong{Details}.
#' 
#' @param spect The spectrum to be reconstructed, as single column matrix or
#'  \code{\link[ICAMS]{ICAMS}} catalog.
#'  
#' @param all.sigs The matrix or catalog of all signatures of possible interest, 
#' which consist of the signatures for \eqn{H_0} and for the alternative
#' hypotheses.
#' 
#' @param Ha.sigs.indices An integer vector of the indices of the signatures
#' that are in the various \eqn{H_a}'s.
#' 
#' @param eval_f Usually one of \code{\link{ObjFnBinomMaxLHNoRoundOK}}
#'               or \code{\link{ObjFnBinomMaxLHMustRound}}. 
#'               For
#'               background see \code{\link[nloptr]{nloptr}}.
#' 
#' @param m.opts Controls the numerical search for maximum likelihood
#'    reconstructions of \code{spect} plus some additional
#'    flags; see \code{\link{DefaultManyOpts}}.
#'    
#' @param max.mc.cores 
#'   The maximum number of cores to use.
#'   If \code{NULL} defaults to \eqn{2^{n_a} - 1},
#'    where \eqn{n_a} is the length of \code{Ha.sigs.indices} -- except on
#'    MS Windows machines, where it defaults to 1.
#'    
#' @details 
#' Let \eqn{H_0} be the likelihood that
#' the signatures specified by \code{all.sigs[, -Ha.sigs.indicies, drop = FALSE]}
#' generated the observed spectrum, \code{spect}.  
#' For each non-empty subset, \eqn{S},
#' of \code{Ha.sigs.indices} let \eqn{H_a}
#' be the likelihood that all the signatures in \eqn{H_0}
#' plus the signatures specified by \eqn{S} generated \code{spect}.
#' Return a list with the results of likelihood ratio tests of 
#' all \eqn{H_a}'s against \eqn{H_0}.
#' 
#' @return A list with two elements: \describe{
#' 
#' \item{\code{H0.info}}{contains the sub-elements \describe{
#'    \item{\code{loglh}}{The log likelihood associated with \eqn{H_0}.}
#'    \item{\code{exposure}}{The signature attributions (exposures) corresponding
#'         to the \eqn{H_0} log likelihood.}
#'    \item{\code{everything.else}}{A sub-list with information on the output
#'         of the numerical optimization that provided \code{loglh}.}
#' }}
#' 
#' \item{\code{all.Ha.info}}{A list with one sub-element for each non-empty subset of
#'   \code{Ha.sigs.indices}. Each sub-element is a list with elements that include \describe{
#'   
#'   \item{\code{sigs.added}}{The identifiers of the (additional) signatures
#'        tested.}
#'        
#'   \item{\code{p}}{The \eqn{p} value for the likelihood-ratio test.} This
#'   this \eqn{p} value can be \code{NaN} when the likelihoods of (\eqn{H_0} and \eqn{H_a})
#'   are both \code{-Inf}. This can occur if there are are mutation classes in the spectra
#'   that are > 0 but that have 0 probability in all the available input signatures.
#'   This is unlikely to occur, since most spectra have non-0 (albeit very small)
#'   probabilities for most mutation classes. This is not an error is using 
#'   \code{eval_f = ObjFnBinomMaxLHNoRoundOK}. However, if \code{p == NaN}
#'   when using \code{eval_f = ObjFnBinomMaxLHMustRound}, switch to 
#'   \code{ObjFnBinomMaxLHNoRoundOK}.
#'   
#'   \item{\code{df}}{The degrees of freedom of the likelihood-ratio test 
#'      (equal to the number of signatures in \code{sigs.added}).}
#'      }
#'   
#'   
#' }}
#' 
#' 
#' WARNING: tests all non-empty subsets of \code{Ha.sigs.indices}, so will get very slow for
#' large numbers of \code{Ha.sigs.indices}.
#' 
#' @export

AnySigSubsetPresent <- 
  function(spect, 
           all.sigs,          
           Ha.sigs.indices, 
           eval_f = mSigAct::ObjFnBinomMaxLHNoRoundOK,  
           m.opts,
           max.mc.cores = NULL) {
    
    H0.sigs <- all.sigs[ , -Ha.sigs.indices, drop = FALSE]  # H0.sigs a matrix of sigs

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
    new.subsets <- 2^sets::as.set(Ha.sigs.indices) - 2^sets::as.set(c())
    
    inner.fn <- function(sigs.to.add) {
      df <- sets::set_cardinality(sigs.to.add) # degrees of freedom for likelihood ratio test
      Ha.info <- one.lh.and.exp(spect  = spect,
                                sigs   = cbind(all.sigs[ , unlist(sigs.to.add), drop = FALSE ], H0.sigs),
                                eval_f = eval_f,
                                m.opts = m.opts)
      statistic <- 2 * (Ha.info$loglh - H0.lh)
      p <- stats::pchisq(statistic, df, lower.tail = FALSE)
      return(list(sigs.added = paste(colnames(all.sigs)[unlist(sigs.to.add)], collapse = ","),
                  statistic  = statistic,
                  df         = df,
                  p          = p,
                  base.sigs  = paste(non.0.exp, collapse = ","),
                  Ha.info    = Ha.info))
    }
    
    if (is.null(max.mc.cores)) {
       max.mc.cores <- 2^length(Ha.sigs.indices) - 1      
    }
    
    mc.cores <- Adj.mc.cores(max.mc.cores) # Set to 1 if OS is MS Windows

    out <- parallel::mclapply(new.subsets, inner.fn, mc.cores = mc.cores)

    return(list(H0.info = start, all.Ha.info = out))
  }
