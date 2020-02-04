my.PlotCatalog.SBS96Catalog<-function (catalog, plot.SBS12, cex = 0.8, grid = TRUE, upper = TRUE, extraY=1,
                                       xlabels = TRUE) 
{
  stopifnot(dim(catalog) == c(96, 1))
  stopifnot(rownames(catalog) == ICAMS::catalog.row.order$SBS96)
  class.col <- c("#0000ff", "#000000", "#ff4040", "#838383", 
                 "#40ff40", "#ff667f")
  cols <- rep(class.col, each = 16)
  maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  num.classes <- length(catalog)
  if (attributes(catalog)$catalog.type == "density") {
    bp <- barplot(catalog[, 1] * 1e+06, xaxt = "n", yaxt = "n", 
                  xaxs = "i", xlim = c(-1, 230), lwd = 3, space = 1.35, 
                  border = NA, col = cols, ylab = "mut/million", cex.lab = 0.8)
    ymax <- max(catalog[, 1] * 1e+06)
  }
  else if (attributes(catalog)$catalog.type == "counts") {
    ymax <- 4 * ceiling(max(max(catalog[, 1]), 10)/4)
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", xlim = c(-1, 
                                                                 230), ylim = c(0, ymax), xaxs = "i", lwd = 3, space = 1.35, 
                  border = NA, col = cols, ylab = "counts", cex.lab = 0.8)
    for (i in 1:6) {
      j <- 16 + 16 * (i - 1)
      k <- 1 + 16 * (i - 1)
      text(bp[j], ymax * 1.2, labels = sum(catalog[k:(16 * 
                                                        i), ]), adj = c(1, 1), xpd = NA, cex = cex)
    }
  }
  else if (attributes(catalog)$catalog.type %in% c("counts.signature", 
                                                   "density.signature")) {
    yaxislabel <- ifelse(attributes(catalog)$catalog.type == 
                           "counts.signature", "counts proportion", "density proportion")
    ymax <- max(catalog[, 1])*extraY
    bp <- barplot(catalog[, 1], xaxt = "n", yaxt = "n", xaxs = "i", 
                  xlim = c(-1, 230), lwd = 3, space = 1.35, border = NA,ylim=c(0,ymax), 
                  col = cols, ylab = yaxislabel, cex.lab = 0.8)
  }
  if (grid) {
    segments(bp[1] - 0.5, seq(ymax/4, ymax, ymax/4), bp[num.classes] + 
               0.5, seq(ymax/4, ymax, ymax/4), col = "grey35", lwd = 0.25)
  }
  segments(bp[1] - 0.5, 0, bp[num.classes] + 0.5, 0, col = "grey35", 
           lwd = 0.25)
  y.axis.values <- seq(0, ymax, ymax/4)
  if (attributes(catalog)$catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  }
  else {
    y.axis.labels <- y.axis.values
  }
  if (grid) {
    text(-0.5, y.axis.values, labels = y.axis.labels, las = 1, 
         adj = 1, xpd = NA, cex = cex)
  }
  else {
    Axis(side = 2, at = y.axis.values, las = 1, cex.axis = cex, 
         labels = FALSE)
    text(-3.5, y.axis.values, labels = y.axis.labels, cex = cex, 
         las = 1, adj = 1, xpd = NA)
  }
  text(bp[1], ymax * 1.08, labels = colnames(catalog), xpd = NA, 
       cex = cex, font = 2, adj = c(0, 0))
  if (xlabels) {
    xlabel.idx <- seq(1, 96, by = 4)
    label <- c("A", "C", "G", "T")
    text(bp[xlabel.idx], -ymax/7, labels = label, cex = cex, 
         adj = 0.5, xpd = NA)
    x <- list(bp[xlabel.idx], bp[xlabel.idx + 1], bp[xlabel.idx + 
                                                       2], bp[xlabel.idx + 3])
    y <- c(-ymax/3.5, -ymax/2.8, -ymax/2.39, -ymax/2.1)
    for (i in 1:4) {
      text(x[[i]], y[i], labels = label[i], cex = cex, 
           adj = 0.5, xpd = NA)
    }
    text(1.5, -ymax/7, labels = "preceded by 5'", pos = 2, 
         xpd = NA, cex = cex)
    text(1.5, -ymax/3.5, labels = "followed by 3'", pos = 2, 
         xpd = NA, cex = cex)
  }
  if (upper) {
    x.left <- bp[seq(1, 81, 16)]
    x.right <- bp[seq(16, 96, 16)]
    rect(xleft = x.left, ymax * 1.28, xright = x.right, ymax * 
           1.3, col = class.col, border = NA, xpd = NA, adj = 0.5)
    text((x.left + x.right)/2, ymax * 1.38, labels = maj.class.names, 
         xpd = NA)
  }
  return(list(plot.success = TRUE, plot.object = bp))
}

My.PlotCatalog.IndelCatalog<-function (catalog, plot.SBS12, cex, grid, upper, xlabels,extraY=1) {
  stopifnot(dim(catalog) == c(83, 1))
  indel.class.col <- c("#fdbe6f", "#ff8001", "#b0dd8b", "#36a12e", 
                       "#fdcab5", "#fc8a6a", "#f14432", "#bc141a", "#d0e1f2", 
                       "#94c4df", "#4a98c9", "#1764ab", "#e2e2ef", "#b6b6d8", 
                       "#8683bd", "#61409b")
  num.classes <- length(catalog)
  cols <- rep(indel.class.col, c(6, 6, 6, 6, 6, 6, 6, 6, 6, 
                                 6, 6, 6, 1, 2, 3, 5))
  if (attributes(catalog)$catalog.type == "counts") {
    ymax <- 4 * ceiling(max(max(catalog[, 1]) * 1.3, 10)/4)*extraY
    bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, 
                  xaxt = "n", lwd = 3, space = 1.35, border = NA, col = cols, 
                  xpd = NA, xaxs = "i", yaxt = "n")
    counts <- integer(16)
    for (i in 1:16) {
      idx <- c(6 * 1:12, 73, 75, 78, 83)
      idx2 <- c(8.9, 23, 37.1, 51.2, 65.3, 79.4, 93.5, 
                107.6, 121.7, 135.8, 149.9, 164, 172.2, 175.5, 
                182, 191)
      if (i == 1) {
        counts[i] <- sum(catalog[1:idx[1], 1])
      }
      else {
        counts[i] <- sum(catalog[(idx[i - 1] + 1):idx[i], 
                                 1])
      }
      text(idx2[i], ymax * 0.6, labels = counts[i], cex = 0.68, 
           adj = 1, xpd = NA)
    }
  }
  else if (attributes(catalog)$catalog.type == "counts.signature") {
    ymax <- ifelse(max(catalog[, 1]) * 1.3 > 1, 1, max(catalog[,1]) * 1.3)*extraY
    bp <- barplot(catalog[, 1], ylim = c(0, ymax), axes = FALSE, 
                  xaxt = "n", lwd = 3, space = 1.35, border = NA, col = cols, 
                  xpd = NA, xaxs = "i", yaxt = "n")
  }
  rect(xleft = bp[1] - 1.5, 0, xright = bp[num.classes] + 1, 
       ymax, col = NA, border = "grey60", lwd = 0.5, xpd = NA)
  segments(bp[1] - 1.5, seq(0, ymax, ymax/4), bp[num.classes] + 
             1, seq(0, ymax, ymax/4), col = "grey60", lwd = 0.5, xpd = NA)
  maj.class.names <- c("1bp deletion", "1bp insertion", ">1bp deletions at repeats\n(Deletion length)", 
                       ">1bp insertions at repeats\n(Insertion length)", "Deletions with microhomology\n(Deletion length)")
  x.left <- bp[c(seq(0, 66, 6), 72, 73, 75, 78) + 1] - 0.5
  x.right <- bp[c(seq(6, 72, 6), 73, 75, 78, 83)] + 0.5
  class.pos <- c((x.left[seq(1, 4, 2)] + x.right[seq(2, 5, 
                                                     2)])/2, (x.left[c(6, 10)] + x.right[c(8, 12)] - 12)/2, 
                 (x.left[13] + x.right[length(x.left)])/2)
  category.lab <- c(rep(c("C", "T"), 2), rep(c("2", "3", "4", 
                                               "5+"), 3))
  category.col <- c(rep(c("black", "white"), 2), rep(c("black", 
                                                       "black", "black", "white"), 3))
  rect(xleft = x.left, ymax * 1.02, xright = x.right, ymax * 
         1.11, col = indel.class.col, border = NA, xpd = NA)
  text((x.left + x.right)/2, ymax * 1.06, labels = category.lab, 
       cex = 0.65, col = category.col, xpd = NA)
  text(class.pos, ymax * 1.27, labels = maj.class.names, cex = 0.75, 
       xpd = NA)
  text(1.5, ymax * 7/8, labels = colnames(catalog), adj = 0, 
       cex = 0.8, font = 2)
  y.axis.values <- seq(0, ymax, ymax/4)
  if (attributes(catalog)$catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
    text(-9, ymax/2, labels = "counts proportion", srt = 90, 
         xpd = NA, cex = 0.8)
  }
  else {
    y.axis.labels <- y.axis.values
    text(-9, ymax/2, labels = "counts", srt = 90, xpd = NA, 
         cex = 0.8)
  }
  text(0, y.axis.values, labels = y.axis.labels, las = 1, adj = 1, 
       xpd = NA, cex = 0.75)
  mut.type <- c(rep(c("1", "2", "3", "4", "5", "6+"), 2), rep(c("0", 
                                                                "1", "2", "3", "4", "5+"), 2), rep(c("1", "2", "3", "4", 
                                                                                                     "5", "6+"), 4), rep(c("0", "1", "2", "3", "4", "5+"), 
                                                                                                                         4), "1", "1", "2", "1", "2", "3", "1", "2", "3", "4", 
                "5+")
  bottom.pos <- c((x.left[1] + x.right[2])/2, (x.left[3] + 
                                                 x.right[4])/2, class.pos[3:length(class.pos)])
  bottom.lab <- c("Homopolymer length", "Homopolymer length", 
                  "Number of repeat units", "Number of repeat units", "Microhomology length")
  rect(xleft = x.left, -ymax * 0.09, xright = x.right, -ymax * 
         0.01, col = indel.class.col, border = NA, xpd = NA)
  text(bp, -ymax * 0.15, labels = mut.type, cex = 0.65, xpd = NA)
  text(bottom.pos, -ymax * 0.27, labels = bottom.lab, cex = 0.75, 
       xpd = NA)
  return(list(plot.success = TRUE, plot.object = bp))
}
