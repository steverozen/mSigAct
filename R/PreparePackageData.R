

Make.sp.sigs.exome <- function() {
  tmp <- mSigAct::sp.sigs
  attr(tmp, "abundance") <- ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$genome$'96'
  sp.sigs.exome <-
    ICAMS::TransformCatalog(
      tmp, 
      target.region = "exome",
      target.abundance = ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$exome$'96')
  usethis::use_data(sp.sigs.exome)
  
}



# Part 2, Spectra from cisplatin exposed HepG2 cells.

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
