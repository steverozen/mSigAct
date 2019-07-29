# tests with synthetic data
# 
# Find signatures at different distances from 
# HepG2.background.info$background.sig

DistancesToSPSigs <- function() {
  sim <-
    apply(mSigAct::sp.sigs, MARGIN = 2, 
          function(ref.sig) {
            lsa::cosine(
              ref.sig, 
              mSigAct::HepG2.background.info$background.sig)})
  return(sort(sim, decreasing = TRUE))
  
}

TestEstimateBackground <-
  function(rep.num, bg.sig.info, bg.contribution) {
  # To figure out the non-background, do we want
  # to have different background levels based
  # on distribution of intensities, and then
  # one level of non-background? Or a range of
<<<<<<< HEAD
  # non background signatures? Probably the second. 
=======
  # non background signatures? Probably the second
>>>>>>> cf69dddfed0f856ab8ffac78c6d8a718a8b0cca4
    
  # Maybe don't run, just put in 3D array of spectra?
  
}

TestOneSigAndContribution <- function(
  bg.sig.info, bg.contribution, num.replicates) {
  
  mc.cores.to.use <-
    ifelse(Sys.info()["sysname"] == "Windows", 1, num.replicates)

  mc.output <-
    parallel::mclapply(1:num.replicates,
             FUN = TestEstimateBackground,
             bg.sig.info     = bg.sig.info,
             bg.contribution = bg.contribution,
             mc.cores = mc.cores.to.use)
  
  # Not done, need to get the important part of 
  # mc.output and check for errors
  return(mc.output)
}

# Set up a grid of HepG2 mixture of 0.1, 0.5, 0.9
# Do 10 replicates at each mixture, draw
# the number HepG2 mutations from the values in HepG2.background.info,
# then back into the total number of mutations.

AllSynTests <- function() {
  sigs.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
  contribution <- c(0.1, 0.5, 0.9)
  dist <- round(DistancesToSPSigs(), digits = 1)
  num.replicates <- 10
  sig.names <- paste0(sigs.to.test, dist[sigs.to.test])
    
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(441441)
  # check that this works -- otherwise clusterSetRNGStream
  
  output <- 
    array(data = NA, 
          dim = c(num.replicates,
                  length(contribution),
                  length(sigs.to.test)),
          dimnames = list(
            replicate = paste0("replicate",
                               formatC(1:num.replicates,
                                       format = "d",
                                       width = 2,
                                       flag = "0")),
            contribution = as.character(contribution),
            signature    = names))
  
  for (ref.sig in sigs.to.test) {
    for (cont in contribution) {
      output[ref.sig, cont, ] <-
        TestOneSigAndContribution(
          bg.sig.info     = mSigAct::HepG2.background.info, 
          bg.contribution = cont,
          num.replicates  = num.replicates)
    }
  }
  return(output)
}
