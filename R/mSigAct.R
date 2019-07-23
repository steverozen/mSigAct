#' Likelihood that observed spectrum was generated from a vector of expected counts.
#'  
#' Returns the sum of the likelihoods that each each element of the
#' spectrum was generated from the corresponding expected count.
#'  
#' @param spectrum An observed spectrum (a numeric vector)
#' 
#' @param expected.counts A vector of expected mutation counts, one
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
           algorithm='NLOPT_LN_COBYLA',
           maxeval=1000, 
           print_level=0,
           xtol_rel=0.001,  # 0.0001,)
           xtol_abs=0.0001)
  {
    # Maybe this is excessively "realistic"; maybe
    # just take mean of spectra.
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
      catalog.type = attr(spectra, "catalog.type", exact = TRUE))

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
    sig <- ICAMS::as.catalog(
      sig,
      region = attr(spectra.as.sigs, "region", exact = TRUE),
      ref.genome = attr(spectra.as.sigs,"ref.genome", exact = TRUE),
      catalog.type = attr(spectra.as.sigs,"catalog.type", exact = TRUE))
    
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
      ref.genome = "hg19",
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
    ref.genome = "hg19",
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
#' max_(b1, target.sig){s1 | b1, target.sig, background.sig, background.sig.nbinom.size, background.sig.count.mu,
#'      background.sig.count.nbinom.size)
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

ThisObjFn <- function(
  est.target.sig.and.total.mut.and.b, # Parameters to optimize
  target.spectra, # Matrix of s1, s2, ....
  background.sig.info # E.g. HepG2.background.info
) {
  return(NULL)
  
}
