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
                     target.sig.name) {
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
  
  bg.mu   <- bg.sig.info$count.nbinom.mu
  size <- bg.sig.info$count.nbinom.size
  
  bg.count <- stats::rnbinom(1, mu = bg.mu, size = size)
  
  # Draw the "target" signature contribution from
  # the negative binomial distribution with mean
  # 1/(1 - bg.contribution) * bg.sig.info$count.nbinom.mu
  # and size bg.sig.info$count.nbinom.size
  
  target.mu <- bg.mu * (1 - bg.contribution)/bg.contribution
  
  target.count <- stats::rnbinom(1, mu = target.mu, size = size)
  
  ref.sig <- mSigAct::sp.sigs[ , target.sig.name, drop = FALSE]
  
  spectrum <-
    bg.count * bg.sig.info$background.sig + target.count * ref.sig
  
  c(target.sig.name = target.sig.name,
       bg.contribution = bg.contribution,
       bg.mu           = bg.mu,
       bg.count        = bg.count, 
       target.mu       = target.mu,
       target.count    = target.count,
       spectrum        = spectrum)
  
}

TestMakeTest <- function() {
  r1 <- MakeTest(1, HepG2.background.info, 0.1, "SBS6")
  r2 <- MakeTest(2, HepG2.background.info, 0.1, "SBS6")
  rx <- rbind(r1, r2)
  return(rx)
}

# Note: https://www.r-bloggers.com/populating-data-frame-cells-with-more-than-one-value/

# Set up a grid of HepG2 mixture of 0.1, 0.5, 0.9
# Do 10 replicates at each mixture, draw
# the number HepG2 mutations from the values in HepG2.background.info,
# then back into the total number of mutations.

MakeAllHepGSynTests <- function() {
  sig.names.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
  contribution <- c(0.1, 0.5, 0.9)
  dist <- round(DistancesToSPSigs(), digits = 1)
  num.replicates <- 10
  sig.names <- paste0(sig.names.to.test, dist[sig.names.to.test])
  sigs.to.test <- mSigAct::sp.sigs[ , sig.names.to.test, drop = FALSE]
    
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(441441)
  # check that this works -- otherwise clusterSetRNGStream
  
  output <- numeric(102)

  for (i in 1:length(sig.names.to.test)) {
    for (cont in contribution) {
      for (replicate in 1:num.replicates) {
        next.row <-
          MakeTest(rep.num         = replicate, 
                   bg.sig.info     = HepG2.background.info,
                   bg.contribution = cont,
                   target.sig.name = sig.names.to.test[i])
        output <- rbind(output, next.row)
      }
    }
  }
  return(output)
}
