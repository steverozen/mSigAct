#' Specifications of background signatures
#'
#
#'
#' @format A list with the elements \enumerate{
#' \item \code{background.sig} The background signature profile.
#' 
#' \item \code{mean.sig} Deprecated -- the background signature
#' profile computed as the mean of the background spectra; this
#' is essentially the same as \code{background.sig}.
#' 
#' \item \code{sig.nbinom.size} The \code{size} argument for
#' \code{\link[stats]{dnbinom}} for sampling error around
#' the components of \code{background.sig}.
#' 
#' \item \code{count.nbinom.mu} The \code{mu} argument
#'  for
#' \code{\link[stats]{dnbinom}} for the distribution
#' of total counts due to \code{background.sig} across
#' replicate exposed clones.
#' 
#' \item \code{count.nbinom.size} The \code{size} argument
#'  for
#' \code{\link[stats]{dnbinom}} for the distribution
#' of total counts due to \code{background.sig} across
#' replicate exposed clones.
#' 
#' }
#' 
#' @name background.info
#' 
#' @examples 
#' # MORE MORE MORE
#' 
"HepG2.background.info"

#' Test experimental spectra
#' 
#' @format XXXX
#' 
#' @examples 
#' # XXXX
"cisplatin.exposed.HepG2.96"


#' Test background spectra
#' 
#' @format XXXX
#' 
#' 
#' @examples 
#' # XXXX
"HepG2.background.spectra"

#' SigProfiler 96 SBS signatures
#' @format XXXXXXX
#' 
#' @examples 
#' # XXXXX
"sp.sigs"

#' Tests of HepG2 background subtraction
#' 
#' @format A data table
#' 
#' @examples 
#' 
#' HepG2.bg.tests.no.noise[1:2, 1:10]
"HepG2.bg.tests.no.noise"

#' Temporary test data
#' 
#' @format Complicated
#' 
#' @examples 
#' # XXXXXX
#' 
"simple.40000.HepG2.tests"


#' Temporary test data
#' 
#' @format Complicated
#' 
#' @examples 
#' # XXXXXX
#' 
"simple.40000.new.sig0"

#' Temporary test data
#' 
#' @format Complicated
#' 
#' @examples 
#' # XXXXXX
"simple.40000.remainder"

#' Control spectra from Kucab et al, 2019
#' 
#' @format An ICAMS catalog with 35 samples.
"kucab.controls"

#' Resampled mutations counts for combinations of 2:4 \code{kucab.controls}.
#' 
#' @format A list with elements 2:4, each of which is a vector of 
#' 10,000 elements, each of which is the mean of total 
#' mutation counts of 2 to 4 control spectra.
"kucab.control.dist"

#' Another temp test ouput
"simple.100000.NLOPT_LN_COBYLA"

#' Another temp test ouput
"simple.200000.NLOPT_GN_DIRECT_L"

#' Another temp test output
"simple.40000.NLOPT_GN_DIRECT_L"