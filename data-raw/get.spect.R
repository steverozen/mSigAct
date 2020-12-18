#!Rscript

# Script to create a 1-sample spectrum catalog from the PCAWG7 
# data package, from the PCAWG platinum data.

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2) {
  stop("Usage: get.spect.R <mutation type> <sample id 1> [<sample id 2> ...]")
}

mutation.type <- args[1]
id <- args[2:length(args)]
message(mutation.type)

known.sig.types <-c("DBS78","ID", "SBS96", "SBS192", "SBS1536")
if (!mutation.type %in% known.sig.types) {
  stop("First argument must be one of ",
          paste(known.sig.types, collapse = ", "))
}
ss <- PCAWG7::spectra$PCAWG[[mutation.type]]

for (dd in id) {
  if (!dd %in% colnames(ss)) stop(dd, " not found")
}

cc <- ss[ , id, drop = FALSE]

message("dim(cc) = ", dim(cc))

message(paste(id, collapse = ", "))

file.id <- gsub("::", "_", id[1])

message("file.id = ", file.id)
ff <- paste0(mutation.type, "_", file.id, ".csv")

cc <- ICAMS::as.catalog(cc)

ICAMS::WriteCatalog(cc, file = ff)