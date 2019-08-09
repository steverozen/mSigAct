#' Likelihood that observed spectrum was generated from a vector of expected counts.
#'  
#' Returns the sum of the likelihoods that each each element of the
#' spectrum was generated from the corresponding expected count.
#'  
#' @param spectrum An observed spectrum (a numeric vector)
#' 
#' @param expected.counts A vector of (integer) expected mutation counts, one
#' expected count for each mutation type. We want to know the 
#' likelihood that this model generated the observed 
#' spectrum, assuming each mutational types generates counts according to
#' a negative binomial distribution with 
#' the given \code{expected.counts} (argument \code{mu}
#' to \code{\link[stats]{dnbinom}}) and dispersion parameter
#' \code{nbinom.size}.
#' 
#' @param nbinom.size The \code{size} parameter that
#' governs dispersion. See \code{\link[stats]{dnbinom}}.
#' Smaller values correspond to larger dispersion.
#'
#' @return  \code{log(likelihood(spectrum | expected.counts))}
#' 
#' @importFrom stats dnbinom
#'
#' @export

LLHSpectrumNegBinom <-function(spectrum, expected.counts, nbinom.size) {
  
  stopifnot(length(spectrum) == length(expected.counts))
  loglh <- 0
  for (i in 1:length(spectrum)) { # Iterate over each channel in the
    # spectrum and sum the log
    # likelihoods.
    
    nbinom <- stats::dnbinom(x=spectrum[i],
                             mu = expected.counts[i],
                             size=nbinom.size,
                             log = TRUE)
    
    loglh <-loglh + nbinom
  }
  stopifnot(mode(loglh) == 'numeric' )
  return(loglh)
  # -loglh
}


#' Part of the objective function for one background-only spectrum.
#' 
#' Return the log likelihood of the observed spectrum
#' given the signature, and negative binomial dispersion
#' parameter \code{nbinom.size}.
#' 
#' @keywords internal
#' 
#' @param spectrum A single background spectrum
#' @param signature A signature as a numeric vector
#' @param nbinom.size The \code{size} argument for
#' \code{\link[stats]{dnbinom}}.
#' 
LLHOfSignatureGivenSpectrum <- function(spectrum, signature, nbinom.size) {
    expected.counts <- sum(spectrum) * signature
    LLH <- LLHSpectrumNegBinom(spectrum, expected.counts, nbinom.size)
    return(LLH)
  }


#' Objective function for nloptr for estimating a background signature.
#' 
#' Returns the negative log likelihood
#' that the signature with an
#' estimated negative binomial
#' dispersion parameter,
#' generated the observed spectra given the
#' likelihood of the observed total counts
#' associated with the signature.
#' 

# See line 49 in mSigAct.v0.11-start-here.R

#' @keywords internal
#' 
#' @param sig.and.nbinom.size Concatenation of signature as a vector and
#' the negative binomial dispersion parameter.
#' 
#' @param spectra The observed spectra as an \code{\link[ICAMS]{ICAMS}}
#' \code{catalog}.
#' 
#' @return -1 * the log likelihood that the the input signature,
#' dispersion parameter, and mutation count generated 
#' the observed spectrum

NegLLHOfSignature <- function(sig.and.nbinom.size, spectra) {
  len <- length(sig.and.nbinom.size)
  nbinom.size <- sig.and.nbinom.size[len]
  sig         <- sig.and.nbinom.size[-len]
  loglh <- 0
  for (i in 1:ncol(spectra)) {
    loglh.i <-
      LLHOfSignatureGivenSpectrum(spectra[ , i, drop = FALSE],
                                  sig,
                                  nbinom.size)
    loglh <- loglh + loglh.i
  }
  
  return(-1 * loglh)   
}

#' Build a signature for background extraction from a matrix of spectra.
#' 
#' This function not only produces a signature, but also an
#' estimate of the number of mutations usually generated
#' by the signature and an indication of variability around
#' that estimate.
#' 
#' Currently only workds on SBS 96 signatures.
#' 
#' @param spectra An \code{\link[ICAMS]{ICAMS}} \code{catalog} with 
#' \code{catalog.type = "counts"}.
#' 
#' @param algorithm See \code{link[nloptr]{nloptr}}.
#' @param maxeval See \code{link[nloptr]{nloptr}}.
#' @param print_level See \code{link[nloptr]{nloptr}}.
#' @param xtol_rel See \code{link[nloptr]{nloptr}}.
#' @param xtol_abs See \code{link[nloptr]{nloptr}}.
#' 
#' @return A list with the elements 
#' \enumerate{
#' \item \code{signature} An \code{\link[ICAMS]{ICAMS}}
#' \code{catalog} with
#' \code{catalog.type == "counts.signature"}.
#' \item \code{log10.counts} Mean log base 10 of the 
#' total counts in \code{spectra}
#' \item \code{sd.log10.counts.per.base} Standard deviation of 
#' \code{log10.counts.per.base}.
#' }
#' 
#' @export
#
# Internal notes: this function could estimate the
# dispersion for each channel separately, or we cold
# even find the maximum likelihood estimate of the
# generating signature, the number of 
# counts it produces, and its dispersion parameter.

EstimateSignatureFromSpectraLH <-
  function(spectra,
           algorithm="NLOPT_LN_COBYLA",
           maxeval=1000, 
           print_level=0,
           xtol_rel=0.001,  # 0.0001,)
           xtol_abs=0.0001)
  {
    # Maybe this is excessively "realistic"; maybe
    # just take mean of spectra.
    # It turns out that the result signature is exactly the mean
    # 
    # But this can estimate the negative binomial dispersion
    # parameter too, so perhaps useful
    
    spectra.as.sigs <-
      ICAMS::TransformCatalog(spectra, 
                              target.catalog.type = "counts.signature")
    x0.sig.vec <- rowSums(spectra.as.sigs) / ncol(spectra)
    # Start with a signature that is an average
     
    mean.sig <- matrix(x0.sig.vec, ncol = 1)
    rownames(mean.sig) <- rownames(spectra)
    mean.sig <- ICAMS::as.catalog(
      mean.sig,
      ref.genome   = attr(spectra, "ref.genome",   exact = TRUE),
      region       = attr(spectra, "region",       exact = TRUE),
      catalog.type = "counts.signature")

    x0.sig.and.size <- c(x0.sig.vec, 200)
    
    
    ret <- nloptr::nloptr(
      x0 = x0.sig.and.size,
      eval_f = NegLLHOfSignature,
      lb=rep(0, length(x0.sig.and.size)),
      opts=list(algorithm=algorithm,
                xtol_rel=xtol_rel,
                print_level=print_level,
                maxeval=maxeval),
      spectra=spectra)
    
    len <- nrow(spectra)
    sig <- matrix(ret$solution[1:len], ncol = 1)
    rownames(sig) <- ICAMS::catalog.row.order$SBS96
    sig <- sig / sum(sig)
    sig <- ICAMS::as.catalog(
      sig,
      region = attr(spectra.as.sigs, "region", exact = TRUE),
      ref.genome = attr(spectra.as.sigs,"ref.genome", exact = TRUE),
      catalog.type = "counts.signature")
    
    nbinom.size <- ret$solution[len + 1]
    
    return(list(background.sig = sig, mean.sig = mean.sig, nbinom.size = nbinom.size))
    
  }

#' Make a standard ICAMS SBS 96 catalog for the HepG2 background signature
#' 
#' @keywords internal

MakeHepG2BackgroundPart1 <- function() {
  spectra <-
    ICAMS:::ReadDukeNUSCat192(
      file = system.file(
        "data-raw/spectra.for.background.signatures/HepG2_SC2_background.txt",
        package = "mSigAct"),
      catalog.type = "counts",
      region = "unknown")$cat96
  ICAMS::WriteCatalog(
    spectra,
    file = system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct")
  )
}


#' Read the specified spectra file and estimate the HepG2 background signature
#' 
#' @keywords internal
#' 
MakeHepG2BackgroundPart2 <- function(maxeval) {
  set.seed(3214)
  s1 <- ICAMS::ReadCatalog(
    system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct"),
    region = "genome",
    catalog.type = "counts")

  # Calculate parameters for a negative bionomial 
  # distribution modeling the
  # number of mutations generated by the signature.
  count.nbinom.mu <- mean(colSums(s1))
  count.nbinom.size <- 20 # determined manually TODO(Steve): analytical formula?

  
  ret2 <-
    EstimateSignatureFromSpectraLH(s1, maxeval = maxeval, print_level = 1)

  return(list(background.sig    = ret2[["background.sig"]],
              mean.sig          = ret2[["mean.sig"]],
              sig.nbinom.size   = ret2[["nbinom.size"]],
              count.nbinom.mu   = count.nbinom.mu, 
              count.nbinom.size = count.nbinom.size))
}

#' Estimate a signature minus a background signature
#' 
#' Let the input spectra be s1, s2, ...
#' 
#' Fig target.sig that maximize the likelihoods
#' max_(b1, target.sig){s1 | 
#'             b1, target.sig, background.sig, 
#'             background.sig.nbinom.size, background.sig.count.mu,
#'      background.sig.count.nbinom.size}
#' b1 * background.sig + (total-mut1 - b1) * target.sig * prob(b1), \cr
#' b2 * background.sig + (total-mut2 - b2) * target.sig * prob(b2), \cr
#' ... \cr
#' \cr
#' 
#' Note: I guess you could estimate "background" and "target.sig"
#' at the same time, but then background might be slightly (?)
#' different for each set of input spectra. To do that would
#' include the likelhoods
#' background.s1 | background.sig, background.sig.nbinom.size
#' background.s2 | background.sig, background.sig.nbinom.size
#' and include a maximization over background.sig and
#' background.sig.nbinom.size
#' ....
#' 
#' @keywords internal

ObjFn1 <- function(
  est.target.sig.and.b, # Parameters to optimize
  obs.spectra,     
  bg.sig.info   # E.g. HepG2.background.info
) {
  bg.sig.profile <- bg.sig.info$background.sig
  len.sig <- nrow(bg.sig.profile)
  est.target.sig <- est.target.sig.and.b[1:len.sig]
  b <- est.target.sig.and.b[(1 + len.sig):length(est.target.sig.and.b)]

  loglh <- 0
  for (i in 1:ncol(obs.spectra)) {
    obs.spectrum <- obs.spectra[ , i, drop = FALSE]
    total.obs.count <- sum(obs.spectrum)
    expected.counts <- 
      (bg.sig.profile * b[i]) + (est.target.sig * (total.obs.count - b[i]))
    loglh1.i <- LLHSpectrumNegBinom(
      spectrum        = obs.spectrum,
      expected.counts = expected.counts,
      nbinom.size     = 10)
        # bg.sig.info$sig.nbinom.size) # TODO(Steve) need to test different values as hyperparameter
    
    loglh2.i <- dnbinom(x    = round(b[i]), 
                        mu   = bg.sig.info$count.nbinom.mu,
                        size = bg.sig.info$count.nbinom.size,
                        log  = TRUE)
    
    loglh <- loglh + loglh1.i + loglh2.i
  }
  
  return(-1 * loglh)
}

MeanOfSpectraAsSig <- function(spectra) {
  ctype <- attr(spectra, "catalog.type", exact = TRUE)
  if (ctype == "counts") {
    tctype <- "counts.signature"
  } else if (ctype == "density") {
    tctype <- "density.signature"
  } else {
    stop("Cannot run MeanOfSpectraAsSig when catalog.type is ", ctype)
  }
  
  sigs <- ICAMS::TransformCatalog(spectra, target.catalog.type = tctype)
  
  mean.sig <- apply(X = sigs, MARGIN = 1, mean)

  mean.sig <- matrix(mean.sig, ncol=1)
  rownames(mean.sig) <- ICAMS::catalog.row.order[["SBS96"]]
  
  mean.sig <- 
    ICAMS::as.catalog(mean.sig, catalog.type = tctype, region = "genome")
  
  colnames(mean.sig) <- "mean.spectra.based.sig"
  
  return(mean.sig)
}


FindSignatureMinusBackground <-
  function(spectra,
           bg.sig.info,
           algorithm='NLOPT_LN_COBYLA',
           maxeval=1000, 
           print_level=0,
           xtol_rel=0.001,  # 0.0001,)
           xtol_abs=0.0001,
           start.b.fraction = 0.1) {
  
  sig0 <- rep(1, nrow(spectra)) / nrow(spectra)
  
  # Test
  sig0 <- MeanOfSpectraAsSig(spectra)
  
  # Test 2 was not good ?
  if (FALSE) {
  what.to.subtract <-
    bg.sig.info$count.nbinom.mu * bg.sig.info$background.sig
  sub2 <- matrix(rep(what.to.subtract, ncol(spectra)), ncol = ncol(spectra))
  spectra.remain <- spectra - sub2 
  sig0 <- MeanOfSpectraAsSig(spectra.remain)[ , 1]
  sig0[sig0 < 0] <- 0
  }
  
  b.x0 <- start.b.fraction * colSums(spectra)
  est.target.sig.and.b.x0 <- c(sig0, b.x0)
  
  ret <- nloptr::nloptr(
    x0          = est.target.sig.and.b.x0,
    eval_f      = ObjFn1,
    lb          = rep(0, length(est.target.sig.and.b.x0)),
    ub          = c(rep(1, nrow(spectra)), # Each element of the singature <= 1
                    colSums(spectra)),     # The contribution of the background 
                                           # should not exceed the total count (not sure if this exactly correct)
    opts        = list(algorithm   = algorithm,
                       xtol_rel    = xtol_rel,
                       print_level = print_level,
                       maxeval     = maxeval),
    obs.spectra = spectra,     
    bg.sig.info = bg.sig.info)
  
  return(ret)    
  }


#' Test \code{FindSignatureMinusBackground} on background-only spectra.
#' @keywords internal
Test0 <- function(start.b.fraction = 0.9, maxeval = 10000) {
  
  spectra <- ICAMS::ReadCatalog(
    system.file(
      "data-raw/spectra.for.background.signatures/HepG2.background.96.csv",
      package = "mSigAct"))
    
  res <- FindSignatureMinusBackground(
    spectra = spectra,
    bg.sig.info = mSigAct::HepG2.background.info,
    maxeval=maxeval, 
    print_level=1,
    start.b.fraction = start.b.fraction)
  
  return(res)
  
}

#' Make spectrum catalog from VCFs from cisplatin exposed HepG2
#' @keywords internal
#' @return The catalog
MakeCisplatinCatalogs <- function() {
  files <- dir("tests/testthat/test.data/HepG2_Cis/", full.names = TRUE)
  cats <- ICAMS::StrelkaSBSVCFFilesToCatalog(
    files = files,
    ref.genome = "hg19",
    trans.ranges = 
      ICAMS::trans.ranges.GRCh37,
    region = "genome")
  ICAMS::WriteCatalog(cats$catSBS96, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS96.csv")
  ICAMS::WriteCatalog(cats$catSBS192, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS192.csv")
  ICAMS::WriteCatalog(cats$catDBS78, 
                      file = "tests/testthat/test.data/HepG2_Cis/DBS78.csv")
  ICAMS::PlotCatalogToPdf(cats$catSBS96, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS96.pdf")
  ICAMS::PlotCatalogToPdf(cats$catSBS192, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS192.pdf")
  ICAMS::PlotCatalogToPdf(cats$catDBS78, 
                      file = "tests/testthat/test.data/HepG2_Cis/DBS78.pdf")
  
  return(cats)
}

Solution2SigVec <- function(solution, sig.number = 96) {
  sig.vec <- solution[1:sig.number]
  return(sig.vec / sum(sig.vec))
}


Solution2Signature <- function(solution, 
                                sig.number = 96,
                                ref.genome = NULL,
                                region = "genome") {
  sig <- matrix(Solution2SigVec(solution), ncol = 1)
  rownames(sig) <- ICAMS::catalog.row.order[["SBS96"]]
  colnames(sig) <- "Inferred.sig"
  sig <- ICAMS::as.catalog(
    sig, ref.genome = ref.genome,
    region = region,
    catalog.type = "counts.signature"
  )
  return(sig)
}


Nloptr2BGMutationCounts <- function(nloptr.retval, sig.number = 96) {
  solution <- nloptr.retval$solution
  return(solution[(sig.number + 1):length(solution)])
}


Nloptr2Signature <- function(nloptr.retval, sig.number = 96) {
  return(Solution2Signature(nloptr.retval$solution, sig.number))
}

Nloptr2ObjFnValue <- function(nloptr.retval) {
  return(nloptr.retval$objective)
}

PlotFactorizations <- function(out.dir,
                               spectra,
                               bg.sig.info,
                               solution,
                               sig.number = 96,
                               ref.genome = NULL,
                               region = "genome")
{
  if (!dir.exists(out.dir)) {
    if (!dir.create(out.dir, recursive = TRUE)) {
      stop("Cannot create ", out.dir)
    }
  }
  sig <- Solution2Signature(solution,
                             sig.number,
                             ref.genome,
                             region)
  
  b <- solution[(sig.number + 1):length(solution)]
  if (length(b) != ncol(spectra)) {
    stop("The number of estimates of the contribution of the target sequence (",
         length(b), ")\n",
         "does not match the number of input spectra (", ncol(spectra), ")")
  }
  total.counts <- colSums(spectra)
  for (i in 1:ncol(spectra)) {
    bg.counts <- round(b[i] * mSigAct::HepG2.background.info$background.sig)
    attr(bg.counts, "catalog.type") <- "counts"
    sig.counts <- round((total.counts[i] - b[i]) * sig)
    attr(sig, "catalog.type") <- "counts"
    tmp <- cbind(spectra[ , i, drop = FALSE],
                 bg.counts,
                 sig.counts,
                 spectra[ , i, drop = FALSE] - bg.counts)
    
    # TODO(Steve) average the spectra minus the counts and see
    # what they look like
    # TODO(get the pcawg signatures and add them in at different
    # concentrations, with and without noise)
    name <- colnames(spectra)[i]
    colnames(tmp) <- c("Orig", "BG", "Exp*Sig", "Orig-BG")
    ICAMS::PlotCatalogToPdf(tmp, paste0(out.dir, "/", name, ".pdf"))
  }
  return(data.frame(sample = colnames(spectra),
                    spectrum.count = total.counts,
                    bg.count  = b,
                    target.sig.count = total.counts - b))
}

PlotTest0 <- function(out.dir, retval)  {
  PlotFactorizations(out.dir,
                     spectra = mSigAct::HepG2.background.spectra,
                     bg.sig.info = mSigAct::HepG2.background.info,
                     solution = retval$solution)
}

Test1 <- function(start.b.fraction = 0.1, maxeval = 10000) {
  res <- FindSignatureMinusBackground(
    spectra = mSigAct::cisplatin.exposed.HepG2.96,
    bg.sig.info = mSigAct::HepG2.background.info,
    maxeval=maxeval, 
    print_level=1,
    start.b.fraction = start.b.fraction)
  
  return(res)
}

PlotTest1 <- function(out.dir, test1.retval)  {
  PlotFactorizations(out.dir,
                     spectra = mSigAct::cisplatin.exposed.HepG2.96,
                     bg.sig.info = mSigAct::HepG2.background.info,
                     solution = test1.retval$solution)
}

#' Helper function: is the set 'probe' a superset of any set in 'background'?
#' 
#' @keywords internal
is.superset.of.any <- function(probe, background) {
  for (b in background) {
    if (sets::set_is_proper_subset(b, probe)) return(TRUE)
  }
  return(FALSE)
}

#' Find known signatures that can most sparsely reconstruct a spectrum.
#' 
#' @param spect The spectrum to be reconstructed.
#' 
#' @param sigs The signatures to use in reconstruction.
#' 
#' @param max.level The larges number of signatures to consider discarding
#' (this is not quite correct)
#' 
#' @param p.thresh The maximum p value based on which it is decided
#' to retain a signature in a reconstruction.
#' 
#' @param trace If > 0 print out information on progress.
#' 
#' @param obj.fun The objective function for nloptr ????
#' 
#' @param nbinom.size The \code{size} parameter to \code{XXXXX}.
#' 
#' @return An assignment of signature as XXXXXXXX

SparseAssignActivity <- function(spect, sigs,
                                 max.level=5,
                                 p.thresh=0.05, trace=0,
                                 obj.fun,
                                 nbinom.size) {
  mode(spect) <-  'numeric'
  start <- one.lh.and.exp(spect, sigs, trace=0,
                          obj.fun=obj.fun,
                          nbinom.size=nbinom.size)
  lh.w.all <- start$loglh  # The likelihood with all signatures
  start.exp <- start$exposure
  non.0.exp.index <- which(start.exp > 0.5)
  if (trace > 0) {
    cat('Starting with',
        paste(names(start.exp)[non.0.exp.index], collapse=','),
        '\n')
    cat('max.level =', max.level, '\n')
  }
  if (length(non.0.exp.index) < 2) {
    if (trace > 0) cat('returning, only', length(non.0.exp.index),
                       'non-0 exposures\n')
    return(start.exp)
  }
  max.level <- min(max.level, length(non.0.exp.index) - 1)
  
  c.r.s <- sets::set() # subsets of the signature indices that cannot be removed
  
  max.p <- numeric(max.level)
  best.subset <- non.0.exp.index
  best.try <- start
  best.try.is.start <- T
  
  for (df in 1:max.level) {
    if (trace > 0) cat('df =', df, '\n')
    subsets <- sets::set_combn(non.0.exp.index, df)
    for (subset in subsets) {
      # subset is of class set (package sets)
      if (is.superset.of.any(subset, c.r.s)) next
      subset.to.remove.v <- as.numeric(subset)
      subset.name <- paste(names(start.exp)[subset.to.remove.v], collapse=',')
      tmp.set <- setdiff(non.0.exp.index, subset.to.remove.v)
      try.sigs <- sigs[ , tmp.set]
      
      if (length(tmp.set) == 1) {
        try.sigs <- as.matrix(try.sigs)
        rownames(try.sigs) <- rownames(sigs)
        if (trace > 0) cat("New code\n")
      }
      
      # Get the max lh for try.sig / Get the maximum likelihood reconstruction using try.sig
      try <- one.lh.and.exp(spect, try.sigs, trace=0,
                            obj.fun=obj.fun,
                            nbinom.size=nbinom.size)
      # try contains maxlh and exposure
      statistic <- 2 * (lh.w.all - try$loglh)
      chisq.p <- pchisq(statistic, df, lower.tail=F)
      if (trace > 0) {
        cat('Trying', subset.name, 'p =', chisq.p, '; ')
        if (trace > 1) cat('loglh =', try$loglh, '; statistic  =', statistic, ';')
      }
      if (chisq.p > p.thresh) {
        # This an acceptable set of exposures
        if (trace > 0) cat (' acceptable;')
        if (chisq.p > max.p[df]) {
          # This is the best so far
          max.p[df] <- chisq.p
          best.subset <- tmp.set
          best.try <- try
          best.try.is.start <- F
          if (trace > 0) cat('best\n')
        } else {
          if (trace > 0) cat('not best\n')
        }
      } else {
        c.r.s <- sets::set_union(c.r.s, sets::set(subset))
        if (trace > 0)  {
          cat('cannot remove\n')
        }
      }
    } # end for (subset in subsets)
    if (max.p[df] == 0) break
  } # end for df in 1:max.level
  
  # Need to return the exposures in the context of the orginal signatures matrix
  out.exp <- numeric(ncol(sigs)) #all zeros
  names(out.exp) <- colnames(sigs)
  if (best.try.is.start) {
    out.exp <- start.exp
  } else {
    out.exp[best.subset] <- best.try$exposure
  }
  if (trace > 0) {
    cat('max.p =', paste(max.p, collapse = ', '), '\n')
    cat('Ending with',
        paste(names(start.exp)[best.subset], collapse=','),
        '\n')
  }
  stopifnot(abs(sum(out.exp) - sum(spect)) < 1)
  out.exp
}

