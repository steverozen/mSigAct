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


MakeTest <- function(replicate, 
                     bg.sig.info,
                     bg.contribution,
                     target.sig.name,
                     dist.to.bg) {
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
  
  return(cbind(data.frame(target.sig.name = target.sig.name),
               matrix(c( dist.to.bg      = dist.to.bg,
                         bg.contribution = bg.contribution,
                         replicate       = replicate,
                         bg.mu           = bg.mu,
                         bg.count        = bg.count, 
                         target.mu       = target.mu,
                         target.count    = target.count,
                         spectrum        = spectrum),
                      nrow = 1)))
}

TestMakeTest <- function() {
  r1 <- MakeTest(1, HepG2.background.info, 0.1, "SBS6")
  r2 <- MakeTest(2, HepG2.background.info, 0.1, "SBS6")
  rx <- rbind(r1, r2)
  return(rx)
}

# Note: https://www.r-bloggers.com/populating-data-frame-cells-with-more-than-one-value/

# TODO Update this comment Set up a grid of HepG2 mixture of 0.1, 0.5, 0.9
# Do 10 replicates at each mixture, draw
# the number HepG2 mutations from the values in HepG2.background.info,
# then back into the total number of mutations.

#' Make a "grid" of tests, with various "target signatures" and contribution from background
#'
#' @keywords internal
#' 
#' @return A \code{data.frame}, of which each row is a synthetic
#' spectrum, and the last 96 columns are the signatures.
MakeSynTestGrid <-
  function(sig.names.to.test, contribution, bg.sig.info,
           num.replicates = 10) {
  dist <- round(DistancesToSPSigs(), digits = 3)
  sigs.to.test <- mSigAct::sp.sigs[ , sig.names.to.test, drop = FALSE]
    
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(441441)
  # check that this works -- otherwise clusterSetRNGStream
  
  ncol <- 104
  output <- cbind(data.frame("dummy"), matrix(rep(0, ncol - 1), nrow = 1))
  # Each row is a synthetic spectrum.
  # We will have 104 columns
  colnames(output) <- c("target.sig",
                        "dist.to.bg",
                        "bg.contribution",
                        "replicate",
                        "bg.mu",
                        "bg.count", 
                        "target.mu",
                        "target.count",
                        ICAMS::catalog.row.order[["SBS96"]])

  for (i in 1:length(sig.names.to.test)) {
    for (cont in contribution) {
      for (replicate in 1:num.replicates) {
        sig.name.i <- sig.names.to.test[i]
        next.row <-
          MakeTest(replicate       = replicate, 
                   bg.sig.info     = bg.sig.info,
                   bg.contribution = cont,
                   target.sig.name = sig.name.i,
                   dist.to.bg      = dist[sig.name.i])
        colnames(next.row) <- colnames(output)
        output <- rbind(output, next.row)
      }
    }
  }
  output <- output[-1, ]
  rownames(output) <- 
    paste0("row", formatC(1:nrow(output), width = 3, flag = "0"))
  return(output)
}

MakeAndSaveHepG2Tests <- function(num.replicates = 10) {

    sig.names.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
    contribution <- c(0.1, 0.5, 0.9)
    
    HepG2.bg.tests.no.noise <-
      MakeSynTestGrid(sig.names.to.test = sig.names.to.test,
                      contribution      = contribution,
                      bg.sig.info       = HepG2.background.info,
                      num.replicates    = num.replicates)
    usethis::use_data(HepG2.bg.tests.no.noise, overwrite = TRUE)
}