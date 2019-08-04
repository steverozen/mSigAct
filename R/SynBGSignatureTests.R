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
  
  spectrum <- round(spectrum)
  
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
  output <- cbind(data.frame("dummy", stringsAsFactors = FALSE),
                  matrix(rep(0, ncol - 1), nrow = 1))
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

MakeAndSaveHepG2Tests <- function(num.replicates = 20) {

    sig.names.to.test <- c("SBS40", "SBS4", "SBS58", "SBS6")
    contribution <- c(0.1, 0.5, 0.9)
    
    HepG2.bg.tests.no.noise <-
      MakeSynTestGrid(sig.names.to.test = sig.names.to.test,
                      contribution      = contribution,
                      bg.sig.info       = mSigAct::HepG2.background.info,
                      num.replicates    = num.replicates)
    usethis::use_data(HepG2.bg.tests.no.noise, overwrite = TRUE)
}

TestTable2TestInput <- function(test.table, 
                                num.replicates,
                                num.spectra.per.replicate) {
  sig.names <- unique(test.table$target.sig)
  contributions <- unique(test.table$bg.contribution)
  output <- list()
  for (sig in sig.names) {
    for (cont in contributions) {
      sub.table <-
        test.table[test.table$target.sig == sig &
                     test.table$bg.contribution == cont, ]
      # Temporary, only one replicate
      s.range <- 1:num.spectra.per.replicate
      spectra <-
        t(sub.table[s.range, 9:ncol(sub.table)])
      spectra <-
        ICAMS::as.catalog(spectra, region = "genome",
                          catalog.type = "counts")
      output[[paste0(sig, "_", cont, "_", paste(s.range, collapse = ","))]] <- 
        list(spectra = spectra, test.rows = sub.table[s.range, ])
    }
  }
  return(output)
}

RunTests <- function(test.table, 
                     num.replicates,
                     num.spectra.per.replicate,
                     bg.info,
                     out.dir,
                     mc.cores = 10,
                     maxeval = 1000,
                     xtol_rel=0.001/10,  # 0.0001,)
                     xtol_abs=0.0001/10,
                     print_level = 0,
                     # algorithm = "NLOPT_LN_COBYLA", #
                     algorithm = "NLOPT_GN_DIRECT_L", #
                     verbose = FALSE) {
  
  if (!dir.exists(out.dir)) {
    if (!dir.create(out.dir, recursive = TRUE))
      stop("Cannot create ", out.dir)
  }
  

  if (Sys.info()["sysname"] == "Windows" &&
      mc.cores > 1) {
    message("On Windows, changing mc.cores from ", mc.cores, " to 1")
    mc.cores <- 1
  }
  
  message("Using ", mc.cores, " cores\n")
  
  test.input.list <-
    TestTable2TestInput(test.table,
                        num.replicates,
                        num.spectra.per.replicate)
  
  Run1Test <- function(test.name) {

    if (verbose) message(test.name)
    
    spectra <- test.input.list[[test.name]][["spectra"]]
    
    retval <-
      FindSignatureMinusBackground(spectra,
                                   bg.sig.info = bg.info,
                                   maxeval = maxeval, 
                                   print_level = print_level,
                                   xtol_rel=xtol_rel,  # 0.0001,)
                                   xtol_abs=xtol_abs,
                                   algorithm = algorithm,
                                   start.b.fraction = 0.1)
    
    return(list(nloptr.retval = retval, 
                test.name = test.name, 
                algorithm = algorithm,
                input = test.input.list[[test.name]]))
  } # function Run1Test
  
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(411411)
  output <- parallel::mclapply(
    X = names(test.input.list), FUN = Run1Test, mc.cores = mc.cores)
  names(output) <- names(test.input.list)
  
  return(output)
}

RunRunTests <- function(maxeval = 10000) {
  output <- RunTests(
    test.table = mSigAct::HepG2.bg.tests.no.noise, # [1:120, ],
    num.replicates = 1,
    num.spectra.per.replicate = 2,
    bg.info = mSigAct::HepG2.background.info,
    out.dir = "data-raw/big-test-output",
    mc.cores = 10,
    maxeval = maxeval,
    print_level = 0,
    verbose = TRUE,
    algorithm = 
  )
  invisible(output)
}



TestOutput2TestRows <- function(test.output.item) {
  return(test.output.item[["input"]][["test.rows"]])
}

TestOutput2GroundTruthSignature <- function(test.output.item) {
  test.rows <- TestOutput2TestRows(test.output.item)
  sig.names <- test.rows[ , "target.sig"]
  sig.name <- unique(sig.names)
  stopifnot(length(sig.name) == 1)
  return(mSigAct::sp.sigs[ , sig.name, drop = FALSE])
}


EvalOneTest <- function(test.output, bg.info) {
    nloptr.retval <- test.output[["nloptr.retval"]]
    algorithm     <- test.output$algorithm
    if (is.null(algorithm)) algorithm <- "NULL"
    iterations <- nloptr.retval$iterations
    inferred.sig <- Nloptr2Signature(nloptr.retval)
    
    ground.truth.sig <- TestOutput2GroundTruthSignature(test.output)
    
    cos.sim <- lsa::cosine(inferred.sig[ , 1], ground.truth.sig[ , 1])
    colnames(inferred.sig) <- 
      paste(colnames(inferred.sig), round(cos.sim, digits = 3), sep = "_")  
    
    # Compare background and inferred target signature counts 
    # to ground truth
    test.rows <- TestOutput2TestRows(test.output)
    inferred.bg.count <- Nloptr2BGMutationCounts(nloptr.retval)
    precis <- test.rows[ , c(1,3,4,6,8)]
    precis <- cbind(row.name = row.names(precis), precis)
    precis$inferred.bg.count <- inferred.bg.count
    precis$cos.sim           <- rep(cos.sim, nrow(precis))
    precis$nloptr.iterations <- rep(iterations, nrow(precis))
    precis$nloptr.obj.fn.valaue <-
      rep(Nloptr2ObjFnValue(nloptr.retval), nrow(precis))
    precis$algoritm <- rep(algorithm, nrow(precis))
    input.spectra <- t(as.matrix(test.rows[ , 9:ncol(test.rows)]))
    input.spectra <- ICAMS::as.catalog(input.spectra,
                                       region = "genome",
                                       catalog.type = "counts")
    bg.sig <- bg.info$background.sig
    to.subtract <- bg.sig %*% inferred.bg.count
    input.spectra.minus.inferred.bg <- input.spectra - to.subtract
    input.spectra.based.sigs <-
      ICAMS::TransformCatalog(input.spectra.minus.inferred.bg,
                              target.catalog.type = "counts.signature")
    colnames(input.spectra.based.sigs) <- 
      paste0(colnames(input.spectra.based.sigs), "-sub-spec-sig")
    
    mean.spectra.based.sig <- 
      MeanOfSpectraAsSig(input.spectra.minus.inferred.bg)
    
    mean.cos.sim <- 
      lsa::cosine(ground.truth.sig[ , 1], mean.spectra.based.sig[ , 1])
    
    precis$cos.sim.of.mean.sig <- rep(mean.cos.sim, nrow(precis))
    
    colnames(mean.spectra.based.sig) <-
      paste0("mean.spectra.based.sig_",
             round(mean.cos.sim, digits = 3))
    
    message("TODO, try just taking the average of the orignal spectra")

    return(
      list(
        cos.sim = cos.sim, 
        inferred.sig = inferred.sig,
        ground.truth.sig = ground.truth.sig,
        bg.sig  = bg.sig,
        precis = precis,
        input.spectra                   = input.spectra,
        input.spectra.minus.inferred.bg = input.spectra.minus.inferred.bg,
        input.spectra.based.sigs        = input.spectra.based.sigs,
        mean.spectra.based.sig          = mean.spectra.based.sig,
        mean.cos.sim                    = mean.cos.sim))
}


# We need to get the ground truth sig from the test output
EvalMultiTest <- function(test.output, bg.info) {
  stopifnot(!is.null(bg.info))
  retval <- lapply(X = test.output, 
                   FUN = EvalOneTest, 
                   bg.info = bg.info)
  return(retval)
  
}

SaveEvaluatedOuput <- function(out.dir, ev.output) {
  
  if (!dir.exists(out.dir)) {
    if (!dir.create(out.dir, recursive = TRUE)) {
      stop("Cannot create ", out.dir)
    }
  }

  csv.append <- FALSE
  csv.colnames <- TRUE
  for (i in 1:length(ev.output)) {
    test.name <- names(ev.output)[i]
    # mydir <- file.path(out.dir, test.name)
    

    ev <- ev.output[[i]]
    sigs <- cbind(ev$bg.sig,
                  ev$ground.truth.sig, 
                  ev$inferred.sig,
                  ev$mean.spectra.based.sig,
                  ev$input.spectra.based.sigs)
    attr(sigs, "catalog.type") <- "counts.signature"
    
    sigs.file <- paste0(out.dir, "/", test.name, ".sigs.pdf" )
    ICAMS::PlotCatalogToPdf(sigs, sigs.file)
    
    spectra <-
      cbind(ev$input.spectra, ev$input.spectra.minus.inferred.bg)
    
    spectra.file <- paste0(out.dir, "/", test.name, ".spectra.pdf")
    ICAMS::PlotCatalogToPdf(spectra, spectra.file)
    
    utils::write.table(ev$precis, 
                      file.path(out.dir, "precis.csv"),
                      col.names = csv.colnames,
                      append    = csv.append,
                      sep       = ",",
                      row.names = FALSE)
    csv.append <- TRUE
    csv.colnames <- FALSE
  }
}

# mfoo <- EvalMultiTest(simple.40000.HepG2.tests, HepG2.background.info)
# SaveEvaluatedOuput("data-raw/ev3", mfoo)

# foo2X <- EvalMultiTest(simple.40000.new.sig0, HepG2.background.info)
# SaveEvaluatedOuput("data-raw/ev.new.sig0x", foo2X)

# simple.40000.remainder <- RunRunTests(maxeval = 40000)
# foo.remainder <- EvalMultiTest(simple.40000.remainder, HepG2.background.info)
# SaveEvaluatedOuput("data-raw/ev.remainder", foo.remainder)

# simple.100000.NLOPT_GN_DIRECT_L <- RunRunTests(maxeval = 100000)
# usethis::use_data(simple.40000.NLOPT_GN_DIRECT_L)
# foo.NLOPT_GN_DIRECT_L <- EvalMultiTest(simple.40000.NLOPT_GN_DIRECT_L, HepG2.background.info)
# SaveEvaluatedOuput("data-raw/ev.NLOPT_GN_DIRECT_L", foo.NLOPT_GN_DIRECT_L)
