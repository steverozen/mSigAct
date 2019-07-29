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


MakeTest <- function(rep.num, 
                     bg.sig.info,
                     bg.contribution,
                     ref.sig) {
  # To figure out the non-background, do we want
  # to have different background levels based
  # on distribution of intensities, and then
  # one level of non-background? Or a range of
  # non background signatures? Probably the second. 
  
  # Maybe don't run, just put in 3D array of spectra?
  
  # Draw the background signature contribution from
  # the negative binomial distribution 
  # with mean bg.sig.info$count.nbinom.mu and
  # size parameters bg.sig.info$count.nbinom.size
  
  mu <- bg.sig.info$count.nbinom.mu
  size <- bg.sig.info$count.nbinom.size
  
  bg.count <- stats::rnbinom(1, mu = mu, size = size)
  
  # Draw the "target" signature contribution from
  # the negative binomial distribution with mean
  # 1/(1 - bg.contribution) * bg.sig.info$count.nbinom.mu
  # and size bg.sig.info$count.nbinom.size
  
  target.count <-
    stats::rnbinom(1, mu = 1/(1 - bg.contribution) * mu, size = size)
  
  spectrum <- bg.count * bg.sig.info$background.sig
              + target.count * ref.sig
  
  list(bg.count = bg.count, target.count = target.count, spectrum = spectrum)
  
}


TestOneSigAndContribution <- function(
  bg.sig.info, bg.contribution, ref.sig, num.replicates) {
  
  mc.cores.to.use <-
    ifelse(Sys.info()["sysname"] == "Windows", 1, num.replicates)

  mc.output <-
    parallel::mclapply(1:num.replicates,
             FUN = MakeTest,
             bg.sig.info     = bg.sig.info,
             bg.contribution = bg.contribution,
             ref.sig         = ref.sig,
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
  sig.names.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
  contribution <- c(0.1, 0.5, 0.9)
  dist <- round(DistancesToSPSigs(), digits = 1)
  num.replicates <- 10
  sig.names <- paste0(sig.names.to.test, dist[sig.names.to.test])
  sigs.to.test <- mSigAct::sp.sigs[ , sig.names.to.test, drop = FALSE]
    
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(441441)
  # check that this works -- otherwise clusterSetRNGStream
  
  output <- 
    array(data = NA, 
          dim = c(num.replicates,
                  length(contribution),
                  length(sig.names)),
          dimnames = list(
            replicate = paste0("replicate",
                               formatC(1:num.replicates,
                                       format = "d",
                                       width = 2,
                                       flag = "0")),
            contribution = as.character(contribution),
            signature    = sig.names))
  
  for (i in 1:length(sigs.to.test)) {
    for (cont in contribution) {
      output[sig.names[i], cont, ] <-
        TestOneSigAndContribution(
          bg.sig.info     = mSigAct::HepG2.background.info, 
          bg.contribution = cont,
          target.sig      = sigs.to.test[i],
          num.replicates  = num.replicates)
    }
  }
  return(output)
}
