if (FALSE) {
# Spectra from HepG2 exposed to various nitrosamines

options(error = browser)
options(warn = 3)
showMethods("getSeq")
# Function: getSeq (package Biostrings)
# x="BSgenome"
# x="FaFile"
# x="FaFileList"
# x="TwoBitFile"
# x="XStringSet"


assign("last.warning", NULL, envir = baseenv())


# Important
SBS.dir <- "data-raw/HepG2-nitrosamines-2019-12-26/SBS"
tmp.files <- grep("\\.vcf$",
                 list.files(SBS.dir,
                            full.names = TRUE), 
                 ignore.case = TRUE, value = TRUE) 
# Important
tmp.cats3 <- 
  ICAMS::StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
  files = tmp.files,
  ref.genome = 
    BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
  trans.ranges = ICAMS::trans.ranges.GRCh37,
  region = "genome",
  names.of.VCFs = sub(".*(N..._cl.).*", "\\1", tmp.files),
  output.file = file.path(SBS.dir, "nitrosamines")
)
}

# Spectra from cisplatin exposed HepG2

#' Make spectrum catalog from VCFs from cisplatin exposed HepG2
#' @keywords internal
#' @return The catalog
MakeCisplatinCatalogs <- function() {
  files <- dir("tests/testthat/test.data/HepG2_Cis/", full.names = TRUE)
  cats <- ICAMS::StrelkaSBSVCFFilesToCatalog(
    files = files,
    ref.genome = "hg19",
    trans.ranges = 
      ICAMS::trans.ranges.GRCh37,
    region = "genome")
  ICAMS::WriteCatalog(cats$catSBS96, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS96.csv")
  ICAMS::WriteCatalog(cats$catSBS192, 
                      file = "tests/testthat/test.data/HepG2_Cis/SBS192.csv")
  ICAMS::WriteCatalog(cats$catDBS78, 
                      file = "tests/testthat/test.data/HepG2_Cis/DBS78.csv")
  ICAMS::PlotCatalogToPdf(cats$catSBS96, 
                          file = "tests/testthat/test.data/HepG2_Cis/SBS96.pdf")
  ICAMS::PlotCatalogToPdf(cats$catSBS192, 
                          file = "tests/testthat/test.data/HepG2_Cis/SBS192.pdf")
  ICAMS::PlotCatalogToPdf(cats$catDBS78, 
                          file = "tests/testthat/test.data/HepG2_Cis/DBS78.pdf")
  
  return(cats)
}

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}
