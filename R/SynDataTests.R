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


# Set up a grid of HepG2 mixture of 0.1, 0.5, 0.9
# Do 10 replicates at each mixture, draw
# the number HepG2 mutations from the values in HepG2.background.info

AllSynTests <- function() {
  sigs.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
  contribution <- c(0.1, 0.5, 0.9)
  dist <- DistancesToSPSigs()
  num.replicates <- 10
  output <- 
    array(data = NA, 
          dim = c(num.replicates,
                  length(contribution),
                  length(sigs.to.test)),
          dimnames = list(
            replicate     = paste0("replicate",
                                   formatC(1:num.replicates,
                                           format = "d",
                                           width = 2,
                                           flag = "0")),
            contribution = as.character(contribution),
            signature    = sigs.to.test))
  for (ref.sig in sigs.to.test) {
    for (HepG2.contribution in contribution) {
      # output[ref.sig, contribution, ] <-
        # TestOneSigAndContribution(
          # ref.sig, HepG2.contribution, num.replicates)
    }
  }
  return(output)
}
