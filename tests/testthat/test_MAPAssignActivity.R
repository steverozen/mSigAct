context("MAPAssignActivity")

test_that("MAPAssignActivity for ID Catalog", {
  skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  
  catalog <- ICAMS::ReadCatalog(file = "testdata/PCAWG7-Prost-AdenoCA-ten-samples.csv")
  sample.index <- 1:2
  catID <- catalog[, sample.index, drop = FALSE]
  ID.sigs <- ICAMS::ReadCatalog(file = "testdata/COSMIC-v3-genome-ID-sigs.csv")
  mutation.type <- "ID"
  cancer.type <- "Prost-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- ID.sigs[, names(sigs.prop), drop = FALSE]
  output.dir <- file.path(tempdir(), "ID")
  
  retval <- MAPAssignActivity(
    spectra                 = catID,
    sigs                    = sigs,
    sigs.presence.prop      = sigs.prop,
    output.dir              = output.dir,
    max.level               = ncol(sigs) - 1,
    p.thresh                = 0.01,
    num.parallel.samples    = 2,
    mc.cores.per.sample     = 30,
    seed                    = 8787)
  
  expect_equal(length(retval), 2)
  
  if (FALSE) {
    tb1 <- tibble::tibble(sig.id = c("ID1", "ID2", "ID3", "ID5", "ID8",  
                                     "ID9", "ID10"), 
                          count = c(ID1 = 60.0957928396838, 
                                    ID2 = 12.5753318347774,  
                                    ID3 = 20.5025441148954, 
                                    ID5 = 74.6080111814531, 
                                    ID8 = 19.7204896759284,  
                                    ID9 = 15.131816682705, 
                                    ID10 = 13.366013670557))
    expect_equal(retval[[1]]$MAP, tb1)
    
    tb2 <- tibble::tibble(sig.id = c("ID1", "ID2", "ID3", "ID4", "ID5",  
                                     "ID8"), 
                          count = c(ID1 = 86.4049290774173, 
                                    ID2 = 5.66462263579074,  
                                    ID3 = 20.901676243643, 
                                    ID4 = 96.6166407779654, 
                                    ID5 = 155.562140415472,  
                                    ID8 = 47.8499908497115))
    expect_equal(retval[[2]]$MAP, tb2)
  }
 
  unlink(output.dir, recursive = TRUE)
  
})