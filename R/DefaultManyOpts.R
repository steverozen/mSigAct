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
#' @export
#'
#' @return A list with the following elements
#' \describe{
#'
#'   \item{global.opts}{A sub-list with several options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the global optimization phase,
#'   including \code{eval_f}, the objective function.
#'   }
#'
#'   \item{local.opts}{A sub-list with several options for \code{\link[nloptr]{nloptr}},
#'   q.v., for the local optimization phase,
#'   including \code{eval_f}, the objective function and the
#'   inequality constraint function \code{eval_g_ineq}}
#'
#'   \item{nbinom.size}{The dispersion parameter for the negative
#'        binomial distribution; smaller is more dispersed.
#'        See \code{\link[stats]{NegBinomial}}.}
#'
#'   \item{trace}{If > 0 print progress messages.}
#' }

DefaultManyOpts <- function() {
  return(list(
    global.opts       = DefaultGlobalOpts(),
    local.opts        = DefaultLocalOpts(),
    nbinom.size       = 5,
    trace             = 0,
    global_eval_f     = ObjFnBinomMaxLHRound,
    local_eval_f      = ObjFnBinomMaxLHRound,
    local_eval_g_ineq = g_ineq_for_ObjFnBinomMaxLH2))
}
