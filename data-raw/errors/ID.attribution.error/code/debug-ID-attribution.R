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
  recon.spect <- round(MAP.out$MAP.recon)
  cossim <- 
    round(MAP.out$MAP.distances$value[MAP.out$MAP.distances$method == "cosine"], 5)
  colnames(recon.spect) <- paste0("reconstructed (cosine similarity = ", cossim, ")")
  
  MAP.exp <- dplyr::arrange(MAP.out$MAP, dplyr::desc(count))
  selected.sigs.names <- MAP.exp$sig.id
  selected.sigs <- sigs[, selected.sigs.names, drop = FALSE]
  colnames(selected.sigs) <- paste0(colnames(selected.sigs), " (exposure = ", 
                                    round(MAP.exp$count), ")")
  list.of.catalogs <- list(spect, recon.spect, selected.sigs)
  PlotListOfCatalogsToPdf(list.of.catalogs, file = file)
}

DebugIDAttribution <- function() {
  file.ID <- "data-raw/errors/ID.attribution.error/data/catID.counts.csv"
  catID <- ICAMS::ReadCatalog(file = file.ID, ref.genome = "hg19", region = "genome")
  
  ID.expo.prop <- mSigAct::ExposureProportions(mutation.type = "ID", 
                                               cancer.type = "Lung-AdenoCA")
  ID.sig.universe <- PCAWG7::signature$genome$ID[, names(ID.expo.prop), drop = FALSE]
  
  if (.Platform$OS.type == "windows") {
    num.of.cores <- 1
  } else {
    num.of.cores <- 50
  }
  
  ID.MAP.out <- mSigAct::MAPAssignActivity1(spect = catID,
                                            sigs = ID.sig.universe, 
                                            sigs.presence.prop = ID.expo.prop, 
                                            max.level = length(ID.expo.prop) - 1,
                                            p.thresh = 0.01,
                                            m.opts = mSigAct::DefaultManyOpts(),
                                            max.mc.cores = num.of.cores)
  
  save(ID.MAP.out, file = "data-raw/errors/ID.attribution.error/output/ID.MAP.out.Rdata")
  PlotMAPResultToPdf(spect = catID, sigs = ID.sig.universe, MAP.out = ID.MAP.out,
                     file = "data-raw/errors/ID.attribution.error/output/catID.pdf")
  
  return(ID.MAP.out)
}
