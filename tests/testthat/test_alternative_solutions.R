context("Testing alternative solutions")

test_that("Identifying ID samples with low reconstruction accuracy", {
  ID.spectra <- PCAWG7::spectra$PCAWG$ID
  ID.sigs <- PCAWG7::signature$genome$ID
  
  SP.ids <- c("SP102561", "SP79391", "SP135278")
  tumor.ids <- PCAWG7::map_SP_ID_to_tumor_type(SP.IDs = SP.ids)
  ID.three.spectra <- ID.spectra[, tumor.ids]
  
  output.dir <- file.path(tempdir(), "ID")
  
  panc.sparse.out <- 
    mSigAct::MAPAssignActivity(spectra = ID.three.spectra[, 1, drop = FALSE],
                               sigs = ID.sigs, 
                               output.dir = output.dir, 
                               max.level = ncol(ID.sigs) - 1,
                               p.thresh = 0.05 / ncol(ID.sigs), 
                               num.parallel.samples = 1, 
                               mc.cores.per.sample = 50, 
                               seed = 5196, 
                               max.subsets = 1e15,
                               use.sparse.assign = TRUE
    )
  
  alt.solutions <- panc.sparse.out$alt.solutions$`Panc-Endocrine::SP102561`
  
  # View all the alternative solutions and check the QP cosine
  # View(alt.solutions)
  # hist(alt.solutions$QP.cosine)
})