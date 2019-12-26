if (FALSE) {
# Spectra from HepG2 exposed to various nitrosamines

# Create the spectra catalogs and plots
# Cannot find the output -- suspect it is in tempdir(), but cannot find it
# Actually, creates the files in dir, and then unlinks them. If run in debugger
# can quit before the unlink and grab the files from dir
tmp.cats <- 
  ICAMS::StrelkaSBSVCFFilesToZipFile(
  dir = "data-raw/nitrosamine-example-data.2019.12.14/SBS",
  file = tempdir(),
  zipfile.name = "~/nitrosamine.SBS",
  ref.genome = "hg19",
  region = "genome",
  trans.ranges = ICAMS::trans.ranges.GRCh37,
  output.file = "nitrosamines",
  names.of.VCFs = c("NDEA.cl1", "NDEA.cl2", 
                    "NMDA.cl1", "NMDA.cl2", 
                    "NPIP.cl1", "NPIP.cl2", 
                    "NPYR.cl1", "NPYR.cl2")
  
  )

debug(ICAMS::StrelkaIDVCFFilesToZipFile)
# debug(BSgenome::getSeq)
options(error = browser)
options(warn = 3)
showMethods("getSeq")
# Function: getSeq (package Biostrings)
# x="BSgenome"
# x="FaFile"
# x="FaFileList"
# x="TwoBitFile"
# x="XStringSet"
ICAMS::StrelkaIDVCFFilesToZipFile(
  dir = "data-raw/nitrosamine-example-data.2019.12.14/ID",
  file = tempdir(),
  zipfile.name = "~/nitrosamine.ID",
  ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5, # "hg19",
  region = "genome",
  output.file = "nitrosamines",
  names.of.VCFs = c("NDEA.cl1", "NDEA.cl2", 
                    "NMDA.cl1", "NMDA.cl2", 
                    "NPIP.cl1", "NPIP.cl2", 
                    "NPYR.cl1", "NPYR.cl2")
  
)

debug(ICAMS:::NormalizeGenomeArg)
### testing binding to "two bit file" name on Arnoud's computer
ICAMS::StrelkaIDVCFFilesToZipFile(
  dir = "data-raw/nitrosamine-example-data.2019.12.14/ID",
  file = tempdir(),
  zipfile.name = "~/nitrosamine.ID",
  ref.genome = "hg19",
  region = "genome",
  output.file = "nitrosamines",
  names.of.VCFs = c("NDEA.cl1", "NDEA.cl2", 
                    "NMDA.cl1", "NMDA.cl2", 
                    "NPIP.cl1", "NPIP.cl2", 
                    "NPYR.cl1", "NPYR.cl2")
  
)


?# Not used
tmp.files <- grep("\\.vcf$",
                 list.files("data-raw/nitrosamine-example-data.2019.12.14/SBS",
                            full.names = TRUE), 
                 ignore.case = TRUE, value = TRUE) 
# Not used
tmp.cats3 <- 
  ICAMS::StrelkaSBSVCFFilesToCatalog( # AndPlotToPdf(
  files = tmp.files,
  ref.genome = "hg19",
  trans.ranges = ICAMS::trans.ranges.GRCh37,
  region = "genome",
  names.of.VCFs = sub(".*(N..._cl.).*", "\\1", tmp.files) # ,
  # output.file = "data-raw/tmp-n/nitrosamines"
)

tmp.cats4 <- 
  ICAMS::StrelkaSBSVCFFilesToCatalog( # AndPlotToPdf(
    files = c("data-raw/spectra.for.background.signatures/MCF-10A-HepG2-background/background_vcfs/HepG2_SC2_cl1_SNVresult.vcf"),
    ref.genome = "hg19",
    trans.ranges = ICAMS::trans.ranges.GRCh37,
    region = "genome",
    names.of.VCFs = c("HepG2_SC2_cl1")) # ,
    # output.file = "data-raw/tmp-n/nitrosamines"

tmp.root <- "data-raw/nitrosamine-example-data.2019.12.14"
tmp.cats5 <- 
  ICAMS::StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
    files = c(ile.path(tmp.root, "/NDEA_cl1.results/", "passed.somatic.snvs.vcf")),
    ref.genome = "hg19",
    trans.ranges = ICAMS::trans.ranges.GRCh37,
    region = "genome",
    names.of.VCFs = c("NDEA_cl1"),
    output.file = tmp.root
  )

tmp.cat6 <- 
  ICAMS:::StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
    files = c(file.path(tmp.root, 
                        "NDEA_cl1.results", "HepG2_NDEA_cl1_SNVresult.vcf"), # passed.somatic.snvs.vcf"),
              file.path(tmp.root, 
                        "NDEA_cl2.results", "HepG2_NDEA_cl2_SNVresult.vcf") # passed.somatic.snvs.vcf")
              
              ),
    ref.genome = "hg19",
    trans.ranges = ICAMS::trans.ranges.GRCh37,
    region = "genome", 
    names.of.VCFs = c("NDEA_cl1", "NDEA_cl1"),
    output.file = file.path(tmp.root, "NDEA2"))

rm(tmp.files)
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
