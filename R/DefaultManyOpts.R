DefaultGlobalOpts <- function() {
  return(
    list(algorithm   = "NLOPT_GN_DIRECT",
         xtol_rel    = 1e-9,
         print_level = 0,
         maxeval     = 10000))
}

DefaultLocalOpts <- function() {
  return(list(algorithm   = "NLOPT_LN_COBYLA",
              xtol_rel    = 1e-15,
              print_level = 0,
              maxeval     = 10000))
}

#' Set default options for many functions, especially \code{\link[nloptr]{nloptr}}.
#' 
#' @param likelihood.dist The probability distribution used to calculate the
#'   likelihood, can be either "multinom" (multinomial distribution) or
#'   "neg.binom" (negative binomial distribution).
#'   
#' @export
#'
#' @return A list with the following elements
#' \describe{
#'
#'   \item{global.opts}{A sub-list with several options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the global optimization phase.}
#'
#'   \item{local.opts}{A sub-list with several options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the local optimization phase.}
#'
#'   \item{nbinom.size}{The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.} 
#'
#'   \item{trace}{If > 0 print progress messages.}
#'   
#'   \item{global_eval_f}{The objective function for the global optimization phase.}
#'   
#'   \item{local_eval_f}{The objective function for the local optimization phase.}
#'   
#'   \item{local_eval_g_ineq}{The inequality constraint function for 
#'   the local optimization phase.}
#' }

DefaultManyOpts <- function(likelihood.dist = "multinom") {
  if (!likelihood.dist %in% c("multinom", "neg.binom")) {
    stop("The value for argument likelihood.dist should be either multinom or neg.binom")
  }
  
  if (likelihood.dist == "multinom") {
    return(list(
      global.opts       = DefaultGlobalOpts(),
      local.opts        = DefaultLocalOpts(),
      trace             = 0,
      global_eval_f     = ObjFnMultinomMaxLH,
      local_eval_f      = ObjFnMultinomMaxLH,
      local_eval_g_ineq = g_ineq_for_ObjFnMultinomMaxLH))
  } else if (likelihood.dist == "neg.binom") {
    return(list(
      global.opts       = DefaultGlobalOpts(),
      local.opts        = DefaultLocalOpts(),
      nbinom.size       = 5,
      trace             = 0,
      global_eval_f     = ObjFnBinomMaxLHRound,
      local_eval_f      = ObjFnBinomMaxLHRound,
      local_eval_g_ineq = g_ineq_for_ObjFnBinomMaxLH2))
  }
}
