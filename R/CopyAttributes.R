#' If from is an ICAMS catalog copy its attributes to to to, otherwise just make to a counts catalog

#' @keywords internal
CopyAttributes <- function(from, to) {
  if (inherits(from, "catalog")) {
  new.catalog <- 
    ICAMS::as.catalog(to,
                      ref.genome     = attr(from, "ref.genome"),
                      region         = attr(from, "region"),
                      abundance      = attr(from, "abundance"),
                      catalog.type   = attr(from, "catalog.type"),
                      infer.rownames = TRUE)
  } else {
    new.catalog <- ICAMS::as.catalog(to)
  }
  return(new.catalog)
}
