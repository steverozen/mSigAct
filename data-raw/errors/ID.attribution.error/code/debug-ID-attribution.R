#' Plot List of catalogs to Pdf
#' 
#' @param list.of.catalogs List of catalogs in \code{\link{ICAMS}} format.
#'
#' @inheritParams ICAMS::PlotCatalogToPdf
#'
#' @keywords internal
PlotListOfCatalogsToPdf <- function(list.of.catalogs, 
                                    file, 
                                    plot.SBS12 = FALSE, 
                                    cex     = 0.8,
                                    grid    = TRUE, 
                                    upper   = TRUE, 
                                    xlabels = TRUE,
                                    ylim    = NULL) {
  old.par.tck.value <- graphics::par("tck")
  # Setting the width and length for A4 size plotting
  grDevices::pdf(file, width = 8.2677, height = 11.6929, onefile = TRUE)
  graphics::par(tck = old.par.tck.value)
  
  num.of.catalogs <- length(list.of.catalogs)
  catalog.type <- attr(list.of.catalogs[[1]], "catalog.type")
  if (nrow(list.of.catalogs[[1]]) == 96) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(4, 5.5, 2, 1), oma = c(1, 1, 2, 1))
  } else if (nrow(list.of.catalogs[[1]]) == 192) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 2, 1, 1))
  } else if (nrow(list.of.catalogs[[1]]) == 78) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(2, 4, 2, 2), oma = c(3, 3, 2, 2))
  } else if (nrow(list.of.catalogs[[1]]) == 83) {
    opar <- graphics::par(mfrow = c(8, 1), mar = c(3, 4, 2.5, 2), oma = c(3, 3, 2, 2))
  } 
  
  on.exit(graphics::par(opar))
  
  for (i in 1:num.of.catalogs) {
    catalog <- list.of.catalogs[[i]]
    num.of.samples <- ncol(catalog)
    
    for (j in 1:num.of.samples) {
      cat <- catalog[, j, drop = FALSE]
      ICAMS::PlotCatalog(cat, plot.SBS12 = plot.SBS12, cex = cex, grid = grid, 
                         upper = upper, xlabels = xlabels, ylim = ylim)
    }
    
  }
  
  grDevices::dev.off()
  invisible(list(plot.success = TRUE))
}

#' @keywords internal
PlotMAPResultToPdf <- function(spect, sigs, MAP.out, file) {
  recon.spect <- round(MAP.out$proposed.reconstruction)
  cossim <- 
    round(MAP.out$reconstruction.distances$value[MAP.out$reconstruction.distances$method == "cosine"], 5)
  colnames(recon.spect) <- paste0("reconstructed (cosine similarity = ", cossim, ")")
  
  MAP.exp <- 
    MAP.out$proposed.assignment[order(MAP.out$proposed.assignment,  decreasing = TRUE), ]
  selected.sigs.names <- names(MAP.exp)
  selected.sigs <- sigs[, selected.sigs.names, drop = FALSE]
  colnames(selected.sigs) <- paste0(colnames(selected.sigs), " (exposure = ", 
                                    round(MAP.exp), ")")
  list.of.catalogs <- list(spect, recon.spect, selected.sigs)
  PlotListOfCatalogsToPdf(list.of.catalogs, file = file)
}

DebugIDAttribution <- function(num.fake = 0) {
  file.ID <- "data-raw/errors/ID.attribution.error/data/catID.counts.csv"
  catID <- ICAMS::ReadCatalog(file = file.ID, ref.genome = "hg19", region = "genome")
  
  ID.expo.prop <- mSigAct::ExposureProportions(mutation.type = "ID", 
                                               cancer.type = "Lung-AdenoCA")

  su <- 
    PCAWG7::signature$genome$ID[, names(ID.expo.prop), drop = FALSE]
  
  if (num.fake > 0) {
    fake.sig <- matrix(rep(0, nrow(su)), ncol = 1)
    rownames(fake.sig) <- rownames(su)
    colnames(fake.sig) <- "fake.ins.t.1.0"
    fake.sig["INS:T:1:0", ] <- 1
    ID.expo.prop <- c(ID.expo.prop, fake.ins.t.1.0 = 0.5)
    su <- cbind(su, 
                ICAMS::as.catalog(fake.sig, 
                                  catalog.type =
                                    attr(su, "catalog.type")))
    if (num.fake > 1) {
      fake.sig <- matrix(rep(1/nrow(su), nrow(su)), ncol = 1)
      rownames(fake.sig) <- rownames(su)
      colnames(fake.sig) <- "fake.ins.rep.2.0"
      # fake.sig["INS:repeats:2:0", ] <- 1
      ID.expo.prop <- c(ID.expo.prop, fake.ins.rep.2.0 = 0.5)
      su <- cbind(su, 
                  ICAMS::as.catalog(fake.sig, 
                                    catalog.type =
                                      attr(su, "catalog.type")))
      
    }
  }
  
  px <- "data-raw/errors/ID.attribution.error.mult/"
  ID.MAP.out <- mSigAct::MAPAssignActivity(spect = catID,
                                            sigs = su, 
                                            sigs.presence.prop = ID.expo.prop, 
                                            max.level = length(ID.expo.prop) - 1,
                                            p.thresh = 0.01,
                                            m.opts = mSigAct::DefaultManyOpts(),
                                            output.dir =   px
  )

  foo <- OptimizeExposureQP(spectrum = catID, signatures = su)
  foo.recon <- ReconstructSpectrum(sigs = su, exp = foo)
  ICAMS::PlotCatalogToPdf(foo.recon, file.path(px, "QPrecon.pdf"))
  cat("foo recon cossim:", cossim(foo.recon, catID))
    
  inferred.exp.max.lh <- ID.MAP.out$proposed.assignment
  # names(inferred.exp.max.lh) <- ID.MAP.out$MAP$sig.id
  cat("Max LH exposures:")
  print(round(inferred.exp.max.lh))
  rr <- ID.MAP.out$proposed.reconstruction
  ICAMS::PlotCatalogToPdf(round(rr), paste0(px, num.fake, file = "lh.recon.pdf"))
  
  spect.v <- catID[ , 1, drop = T] # spectrum as a vector
  trr.fake <- tibble::tibble(rownames(catID), 
                             rr, 
                             spect.v, 
                             rr  -  spect.v)
  data.table::fwrite(trr.fake, file = paste0(px, num.fake, "comp.table.lh.csv"))
  message("Max lh cosine = ", 
          philentropy::distance(
            rbind(spect.v, rr[ , 1, drop = TRUE]), "cosine"))
  negll <- LLHSpectrumNegBinom(spect.v, expected.counts = rr[, 1, drop = TRUE],
                               nbinom.size = DefaultManyOpts()$nbinom.size)
  message("max LL log likelihood = ", negll)
  max.LH.lh.debug <- LLHSpectrumNegBinomDebug(spect.v, expected.counts = rr[, 1, drop = TRUE],
                                      nbinom.size = DefaultManyOpts()$nbinom.size)
  View(max.LH.lh.debug)
  
  
  x2 <- OptimizeExposureQP(spectrum = catID, signatures = su)
  cat("QP exposures:\n")
  print(round(x2))
  ss <- ReconstructSpectrum(su, x2)
  ICAMS::PlotCatalogToPdf(round(ss), paste0(px, num.fake, file = "qp.recon.pdf"))
  tss.fake <- tibble::tibble(rownames(catID), ss, spect.v, ss  -  spect.v)
  data.table::fwrite(tss.fake, file = paste0(px, num.fake, "comp.table.qp.csv"))
  message("QP cossim = ", 
          philentropy::distance(rbind(spect.v, ss[ , 1, drop = TRUE]), "cosine"))
  
  negll <- LLHSpectrumNegBinom(spect.v, expected.counts = ss[, 1, drop = TRUE],
                                nbinom.size = DefaultManyOpts()$nbinom.size)
  message("QP log likelihood = ", negll)
  QP.lh.debug <- LLHSpectrumNegBinomDebug(spect.v, expected.counts = ss[, 1, drop = TRUE],
                                              nbinom.size = DefaultManyOpts()$nbinom.size)
  View(QP.lh.debug)
  
  
  save(ID.MAP.out, file = paste0(px, "ID.MAP.out.Rdata"))
  # PlotMAPResultToPdf(spect = catID, sigs = su, MAP.out = ID.MAP.out,
  #                   file = paste0(px, num.fake, "catID.pdf"))
  
  return(list(ID.MAP.out, max.LH.lh.debug, QP.lh.debug))
}
