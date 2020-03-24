context("test-SparseAssignActivity.R")

# Framework for testing \code{\link{SparseAssignActivity1}}.
# 
#' Tests only for single spectrum.
# 


SparseAssignTest1 <- function(sig.counts, 
                              trace = 0,
                              max.mc.cores = 1) {
  
  set.seed(101010, kind = "L'Ecuyer-CMRG")  
  
  sig.names <- names(sig.counts)
  
  some.sigs  <- 
    PCAWG7::signature$genome$SBS96[ , sig.names, drop = FALSE]
  ref.genome <- attr(some.sigs, "ref.genome", exact = TRUE)
  region     <- attr(some.sigs, "region", exact = TRUE)
  if (is.null(region)) {
    message("Null region, why?")
    region <- "genome"
  }
  
  spect <- round(some.sigs %*% sig.counts)
  spect <-
    ICAMS::as.catalog(
      spect, 
      ref.genome   = ref.genome,
      region       = region,
      catalog.type = "counts")
  nbinom.size <- 5
  
  m.opts <- DefaultManyOpts()
  m.opts$trace <- trace
  
  SA.out <- SparseAssignActivity1(spect       = spect,
                                  sigs         = some.sigs,
                                  eval_f       = ObjFnBinomMaxLHMustRound,
                                  m.opts       = m.opts,
                                  max.mc.cores = max.mc.cores
  )
  
  zeros <- which(SA.out < 0.5)
  if (length(zeros) > 0) {
    SA.out2 <- SA.out[-zeros]
    sig.names2 <- sig.names[-zeros]
  } else {
    SA.out2 <- SA.out
    sig.names2 <- sig.names
  }
  
  if (FALSE) {
    
    polish.out <- Polish(exp       = SA.out2,
                         sig.names = sig.names2,
                         spect     = spect)
  }
  
  recon1 <- round(prop.reconstruct(some.sigs, SA.out))
  
  if (FALSE) {
    recon2 <-
      round(
        prop.reconstruct(
          PCAWG7::signature$genome$SBS96[ , sig.names2, drop = FALSE],
          polish.out))
  }
  
  return(list(soln1       = SA.out,
              # soln2       = polish.out,
              truth       = sig.counts,
              edist1      = EDist2Spect(SA.out, sig.names, spect),
              edist1r     = EDist2SpectRounded(SA.out, sig.names, spect),
              LL1         = LLHSpectrumNegBinom(spect, recon1, nbinom.size) #,
              # LL2         = LLHSpectrumNegBinom(spect, recon2, nbinom.size),
              # edist2      = EDist2Spect(polish.out, sig.names2, spect),
              # ed8st2r     = EDist2SpectRounded(polish.out, sig.names2, spect)
              #, input.spect = spect
  ))
  
}



test_that("SparseAssignActivity (multiple spectra) test 1", {
  input <- c(SBS1 = 1000, SBS22 = 2000)
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <- SparseAssignTest(sig.counts = input, 
                             num.parallel.samples = 1,
                             mc.cores.per.sample  = 1) 
  expected <- matrix(c(999.2467, 1998.7533), ncol = 1)
  rownames(expected) <- names(input)
  colnames(expected) <- "tumor1"
  testthat::expect_equal(retval$exposure, expected, tolerance = 1e-2)
})


test_that("SparseAssignActivity (multiple spectra) test 2", {
  testthat::skip_if_not(Sys.getenv("MSIGACT_TEST_LENGTH") == "long")
  input.exp <- matrix(c(1000, 2000, 0, 2000, 10, 1000), ncol = 3)
  rownames(input.exp) <- c("SBS1", "SBS22")
  # Automated testing on Travic-CI does not allow spawning processes.
  retval <-  SparseAssignTest(sig.counts = input.exp, 
                              num.parallel.samples = 1,
                              mc.cores.per.sample  = 1) 
  expected <- matrix(c(999.2467, 1998.7533, 0, 1996, 10.40123, 1000.59877),
                     ncol = 3)
  colnames(expected) <- paste0("tumor", 1:3)
  rownames(expected) <- rownames(input.exp)
  testthat::expect_equal(retval$exposure, expected, tolerance = 1e-2)
})

