#' Specifications of background signatures
#'
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
#' @source \code{HepG2.background.info} was estimated from \code{\link{HepG2.background.spectra}}.
#'
#' @name background.info
#' 
#' @examples 
#' HepG2.background.info$count.nbinom.mu
#' HepG2.background.info$count.nbinom.size
#' HepG2.background.info$sig.nbinom.size
#' HepG2.background.info$background.sig[1:3, ]
#' 
#' 
#' 
"HepG2.background.info"

#' Background information for MCF-10A cells.
#' 
#' @name background.info
#' 
#' @source \code{MCF10A.background.info} was estimated from \code{\link{MCF10A.background.spectra}}
#' 
"MCF10A.background.info"

#' Background spectra for HepG2 and MCF-10A
#' 
#' @format An \code{\link[ICAMS]{ICAMS}} \code{counts} catalog.
#' 
#' @name background.spectra
#' 
"HepG2.background.spectra"

#' Background spectra for MCF-10A cells
#' 
#' @name background.spectra
#' 
"MCF10A.background.spectra"

#' Example nitrosamine data for background subtraction
#'
#' @format A list of spectra catalogs for nitrosamines.
#'  The names of the catalogs are self-explanatory.
#'  Each catalog has 2 spectra from each of four different
#'  nitrosamines, "NDEA", "NDMA", "NPIP", "NPYR". The
#'  samples names also should be self-explanatory.
#'
"nitrosamine.examples"