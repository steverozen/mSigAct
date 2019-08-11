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

#' Mutational spectra of cisplatin exposed HepG2 cells.
#' 
#' @format An \code{\link[ICAMS]{ICAMS}} \code{counts} catalog.
#' 
"cisplatin.exposed.HepG2.96"


#' Background spectra for HepG2.
#' 
#' @format An \code{\link[ICAMS]{ICAMS}} \code{counts} catalog.
#' 
"HepG2.background.spectra"

#' SigProfiler 96 SBS \code{counts.signatures} under genome abundance.
#' @format An \code{[ICAMS]} \code{counts.spectrum} catalog;
#'  see Alexandrov et al. 
#'  \url{https://www.biorxiv.org/content/10.1101/322859v2}.
#' 
"sp.sigs"

#' SigProfiler 96 SBS \code{counts.signatures} under exome abundance.
#' @format An \code{[ICAMS]} \code{counts.spectrum} catalog;
#'  see Alexandrov et al. 
#'  \url{https://www.biorxiv.org/content/10.1101/322859v2}.
#' 
"sp.sigs.exome"

#' Tests of HepG2 background subtraction
#' 
#' @format A data table
#' 
#' @examples 
#' 
#' HepG2.bg.tests.no.noise[1:2, 1:10]
"HepG2.bg.tests.no.noise"


#' Control spectra from Kucab et al., 2019
#' 
#' @format An \code{\link[ICAMS]{ICAMS}} \code{counts} catalog with 35 samples.
"kucab.controls"

#' Resampled mutations counts for combinations of 2:4 \code{kucab.controls}.
#' 
#' @format A list with elements 2:4, each of which is a vector of 
#' 10,000 elements, each of which is the mean of total 
#' mutation counts of 2 to 4 control spectra.
"kucab.control.dist"


# ==== BELOW HERE EVERYTHING IS TEMPORARY TEST OUTPUT =====

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
"simple.40000.new.sig0"

#' Temporary test data
#' 
#' @format Complicated
#' 
"simple.40000.remainder"

#' Another temp test output
"simple.100000.NLOPT_LN_COBYLA"

#' Another temp test output
"simple.200000.NLOPT_GN_DIRECT_L"

#' Another temp test output
"simple.40000.NLOPT_GN_DIRECT_L"