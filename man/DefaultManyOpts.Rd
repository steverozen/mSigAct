% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DefaultManyOpts.R
\name{DefaultManyOpts}
\alias{DefaultManyOpts}
\title{Set default options for many functions, especially \code{\link[nloptr]{nloptr}}}
\usage{
DefaultManyOpts(likelihood.dist = "multinom")
}
\arguments{
\item{likelihood.dist}{The probability distribution used to calculate the
likelihood, can be either "multinom" (multinomial distribution) or
"neg.binom" (negative binomial distribution).}
}
\value{
A list with the following elements
\describe{

  \item{global.opts}{A sub-list with several options for \code{\link[nloptr]{nloptr}},
  q.v., for the global optimization phase.}

  \item{local.opts}{A sub-list with several options for \code{\link[nloptr]{nloptr}},
  q.v., for the local optimization phase.}

  \item{nbinom.size}{Only appearing if \code{likelihood.dist = "neg.binom"}.
  The dispersion parameter for the negative binomial distribution; smaller is
  more dispersed. See \code{\link[stats]{NegBinomial}}.}

  \item{trace}{If > 0 print progress messages.}
  
  \item{global_eval_f}{The objective function for the global optimization phase.}
  
  \item{local_eval_f}{The objective function for the local optimization phase.}
  
  \item{local_eval_g_ineq}{The inequality constraint function for 
  the local optimization phase.}
  
  \item{likelihood.dist}{The probability distribution used to calculate the likelihood.}
}
}
\description{
Set default options for many functions, especially \code{\link[nloptr]{nloptr}}
}
\examples{
my.opts <- DefaultManyOpts()
my.opts$trace <- 10
}
