#' Specifications of background signatures
#'
#
#'
#' @format A list with the elements \enumerate{
#' \item \code{background.sig} The background signature profile.
#' 
#' \item \code{sig.nbinom.size} The \code{size} argument for
#' \code{\link[stats]{NegBinomial}} for sampling error around
#' the components of \code{background.sig}.
#' 
#' \item \code{count.nbinom.mu} The \code{mu} argument
#'  for
#' \code{\link[stats]{NegBinomial}} for the distribution
#' of total counts due to \code{background.sig} across
#' replicate exposed clones.
#' 
#' \item \code{count.nbinom.size} The \code{size} argument
#'  for
#' \code{\link[stats]{NegBinomial}} for the distribution
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

#' Demonstration data: an \code{\link[ICAMS]{ICAMS}} catalog.
#' 
#' @source Cite Alexandrov et al, 2019
#' 
#' @format An \code{\link[ICAMS]{ICAMS}} catalog with 
#' mutations from 35 biliary tract tumors.
"BTSG.WGS.PCAWG"

#' @section Data sets from Kucab et al.

#' Kucab spectra -- 2 mutation label columns and 324 data columns.
#' 
#' The original format used in Xueqing ZOU's code.
#'
#' @source https://github.com/xqzou/Cell_MutagenSig
#' 
"K_sub_catalogue"

#' Kucab et al., 2019 spectra as a standard \code{ICAMS} count catalog.
#'
#' @source https://github.com/xqzou/Cell_MutagenSig
#' 
"kucab.spectra"

#' Kucab et al., 2019 signatures as a standard \code{ICAMS} counts.signature catalog.
#'
#' @source https://github.com/xqzou/Cell_MutagenSig
#' 
"kucab.sigs"


#' Background information from the controls in Kucab et al., 2019.
#'
#' Estimated from the control spectra in \code{\link{kucab.spectra}}.
#' 
"kucab.control.bg"

#' A deprecated version of \code{\link{kucab.control.bg}}
#' 
"kucab.bg"

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

