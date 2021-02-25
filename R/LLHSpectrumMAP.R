#' Likelihood that 1 observed spectrum was generated from a vector of expected
#' mutation counts using prior information of the signature presence proportions
#' 
#' @param spectrum An observed spectrum (a numeric vector).
#'
#' @param expected.counts A vector of (integer) expected mutation counts, one
#' expected count for each mutation type. We want to know the
#' likelihood that this model generated the observed
#' spectrum, assuming each mutational types generates counts according to
#' a negative binomial distribution with
#' the given \code{expected.counts} (argument \code{mu}
#' to \code{\link[stats]{NegBinomial}}) and dispersion parameter
#' \code{nbinom.size}.
#'
#' @param nbinom.size The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.
#'        
#' @param model Names of sigs present in the MAP exposure. Do not use indices.
#'
#' @param sigs.presence.prop The proportions of samples that contain each
#'    signature. A numerical vector (values between 0 and 1), with names
#'    being a superset of \code{model}.
#'
#' @param verbose If \code{TRUE} print messages under some circumstances.
#'
#' @return \code{log(likelihood(spectrum | expected.counts)) +
#'   log(probability(model | sigs.presence.prop))}, or, in more detail, the sum
#'   of the negative binomial likelihoods that each element of the spectrum
#'   (i.e., the count for each mutation type e.g. ACT > AAT) was generated from
#'   the expected count for that mutation type plus the probability of the
#'   signature model used in the reconstruction given the prior
#'   \code{sigs.presence.prop} .
#'
#' @importFrom stats dnbinom
#'
#' @export
LLHSpectrumMAP <-
  function(spectrum, expected.counts, nbinom.size, model, sigs.presence.prop, 
           verbose = FALSE) {
    
    loglh.of.exp <- LLHSpectrumNegBinom(spectrum = spectrum,
                                        expected.counts = expected.counts,
                                        nbinom.size = nbinom.size,
                                        verbose = verbose)
    
    prob.of.model <- P.of.M(model = model,
                            sigs.presence.prop = sigs.presence.prop)
    
    loglh.MAP <- loglh.of.exp + prob.of.model
    
    return(loglh.MAP)
  }

# https://cran.r-project.org/web/packages/DirichletReg/DirichletReg.pdf
