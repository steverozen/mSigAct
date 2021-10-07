#' @title Read an exposure matrix from a file
#'
#' @param file CSV file containing an exposure matrix.
#'
#' @param check.names Passed to \code{read.csv}.
#' \strong{IMPORTANT}: If \code{TRUE} this will replace the double
#' colon in identifiers of the form <tumor_type>::<sample_id>
#' with two periods (i.e. <tumor_type>..<sample_id>.
#' If \code{check.names} is true, generate a warning
#' if double colons were present.
#'
#' @return Matrix of exposures.
#'
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "mSigAct")
#' exposure <- ReadExposure(file)
ReadExposure <- function(file, check.names = FALSE) {
  if (check.names) {
    headers <- read.csv(file, nrow = 1, header = FALSE, stringsAsFactors = FALSE)
    double.colon <- grep("::", unlist(headers)[-1], fixed = TRUE)
    if (length(double.colon) > 0) {
      warning(":: in sample ID replaced by ..; suggest calling with check.names = FALSE")
    }
  }
  retval <- read.csv(file, row.names = 1, check.names = check.names)
  if (any(duplicated(colnames(retval)))) {
    stop("There is duplicated column name in the input file")
  }
  return(data.matrix(retval))
}

#' @title Write an exposure matrix to a file
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param file File to which to write the exposure matrix (as a CSV file).
#'
#' @param row.names Either a logical value indicating whether the row names of
#'   \code{exposure} are to be written along with \code{exposure}, or a
#'   character vector of row names to be written.
#'
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "mSigAct")
#' exposure <- ReadExposure(file)
#' WriteExposure(exposure, file = file.path(tempdir(), "synthetic.exposure.csv"))
WriteExposure <- function(exposure, file, row.names = TRUE) {
  old.digits <- getOption("digits")
  options(digits = 22)
  write.csv(exposure, file, row.names = row.names)
  on.exit(options(digits = old.digits))
}

#' Sort columns of an exposure matrix from largest to smallest (or vice versa)
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param decreasing If \code{TRUE}, sort from largest to smallest.
#'
#' @return The original \code{exposure} with columns sorted.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "mSigAct")
#' exposure <- ReadExposure(file)
#' exposure.sorted <- SortExposure(exposure)
SortExposure <- function(exposure, decreasing = TRUE) {
  retval <- exposure[, order(colSums(exposure), decreasing = decreasing),
                     drop = FALSE]
  return(retval)
}

#' Plot a single exposure plot
#'
#' @param exposure Exposures as a numerical \code{matrix} (or \code{data.frame})
#'   with signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs. If you want
#'   \code{exposure} sorted from largest to smallest, use
#'   \code{\link{SortExposure}}. Do not use column names that start with
#'   multiple underscores. The exposures will often be mutation counts, but
#'   could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proportions rather than counts.
#'
#' @param xlim,ylim Limits for the x and y axis. If \code{NULL}(default), the
#'   function tries to do something reasonable.
#'
#' @param legend.x,legend.y The x and y co-ordinates to be used to position the
#'   legend.
#'
#' @param cex.legend A numerical value giving the amount by which legend
#'   plotting text and symbols should be magnified relative to the default.
#'   
#' @param cex.yaxis A numerical value giving the amount by which y axis values
#'   should be magnified relative to the default.
#'   
#' @param cex.xaxis A numerical value giving the amount by which x axis values
#'   should be magnified relative to the default. If
#'   \code{NULL}(default), the function tries to do something reasonable.
#'   
#' @param plot.sample.names Whether to plot sample names below the x axis.
#'   Default is TRUE.
#'   
#' @param yaxis.labels User defined y axis labels to be plotted. If
#'   \code{NULL}(default), the function tries to do something reasonable.
#'   
#' @param ... Other arguments passed to \code{\link[graphics]{barplot}}. If
#'   \code{ylab} is not included, it defaults to a value depending on
#'   \code{plot.proportion}. If \code{col} is not supplied the function tries to
#'   do something reasonable.
#'
#' @import graphics
#'
#' @return An \strong{invisible} list whose first element is a logic value
#'   indicating whether the plot is successful. The second element is a numeric
#'   vector giving the coordinates of all the bar midpoints drawn, useful for
#'   adding to the graph.
#'
#' @keywords internal
PlotExposureInternal <-
  function(exposure, # This is actually the exposure "counts"
           plot.proportion   = FALSE,
           xlim              = NULL,
           ylim              = NULL,
           legend.x          = NULL,
           legend.y          = NULL,
           cex.legend        = 0.9,
           cex.yaxis         = 1,
           cex.xaxis         = NULL,
           plot.sample.names = TRUE,
           yaxis.labels      = NULL,
           ...
  ) {
    exposure <- as.matrix(exposure) # In case it is a data frame
    num.sigs  <- nrow(exposure)
    num.samples <- ncol(exposure)
    args <- list(...)

    if (is.null(args$col)) {
      if (num.sigs <= 8) {
        args$col <-
          c("red", "black", "grey", "yellow", "blue", "brown", "green4", "skyblue")
      } else {
        # Lots of signatures; use shading lines to differentiate
        args$col <- grDevices::rainbow(num.sigs, alpha = 1)
      }
    }

    barplot.col <- character(num.sigs)
    barplot.col <-
      suppressWarnings(replace(barplot.col, 1:num.sigs, args$col))

    barplot.density <- barplot.angle <- numeric(num.sigs)

    if (num.sigs <= 12) {
      # Specify the density of shading lines, in lines per inch, for the bars or
      # bar components. Will repeat as needed, non-positive values of density
      # inhibit the drawing of shading lines.
      barplot.density <- rep(-1, num.sigs)

      # Sepcify the slope of shading lines, given as an angle in degrees
      # (counter-clockwise), for the bars or bar components. Will repeat as
      # needed.
      barplot.angle <- rep(0, num.sigs)
    } else {
      # Lots of signatures; use shading lines to differentiate
      barplot.density <-
        suppressWarnings(replace(barplot.density, 1:num.sigs,
                                 c(-1, 50, 50, 50, 50)))
      barplot.angle <-
        suppressWarnings(replace(barplot.angle, 1:num.sigs,
                                 c(0, 45, 135, 45, 135)))
    }
    # For legend, we need to put colour, density and angle in reversed order to make
    # sure this matches order used in barplot.
    legend.col <- rev(barplot.col)
    legend.density <- rev(barplot.density)
    legend.angle <- rev(barplot.angle)

    if (plot.proportion) {
      # Matrix divided by vector goes column-wise, not row-wise, so transpose twice
      plot.what <- t(t(exposure)/colSums(exposure))

      # If some columns in exposure matrix are "padded", i.e. all the values
      # are 0, then the division above will make NaNs. We try to replace all
      # NaN with 0
      plot.what[is.na(plot.what)] <- 0
    } else {
      plot.what <- exposure
    }

    # Set the limits for x axis, leave extra space for plotting legend
    if (is.null(xlim)) {
      xmax <- round(ncol(plot.what) * 1.7)
      xlim <- c(0, xmax)
    } else {
      xmax <- xlim[2]
    }

    if (is.null(ylim)) {
      ymax <- round(max(colSums(plot.what)))
      ylim <- c(0, ymax)
    } else {
      ymax <- ylim[2]
    }
    
    # Check whether user specifies argument for barplot
    ylab <- args$ylab
    args$ylab <- NULL
    if (is.null(ylab)) {
      ylab <- ifelse(plot.proportion,
                     "Proportion of mutations",
                     "Number of mutations")
    }
    
    space <- args$space
    args$space <- NULL
    if (is.null(space)) {
      space <- xmax * 0.01
    }
    
    density <- args$density
    args$density <- NULL
    if (is.null(density)) {
      density <- barplot.density
    }
    
    # Ignore column names, we'll plot them separately to make them fit.
    mp <- do.call(
      barplot,
      args = c(list(height    = plot.what,
                    ylab      = ylab,
                    # The amount of space left before each bar
                    space     = space,
                    xaxs      = "i", # No extra spacing at each end of x axis
                    xaxt      = "n", # Do not plot the X axis
                    yaxt      = "n", # Do not plot the Y axis
                    density   = density,
                    angle     = barplot.angle,
                    border    = "white",
                    xlim      = xlim,
                    ylim      = ylim),
               args))

    # Get locations for y axis annotations
    y.axis.values <- seq(0, ymax, ymax/4)
    
    if(is.null(yaxis.labels)) {
      if (plot.proportion) {
        y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
      } else {
        y.axis.labels <- round(y.axis.values)
      }
    } else {
      y.axis.labels <- yaxis.labels
    }

    # Draw y axis and labels
    Axis(side = 2, at = y.axis.values, las = 1, labels = FALSE)
    text(x = -xmax* 0.03, y = y.axis.values, labels = y.axis.labels,
         las = 1, adj = 1, xpd = NA, cex = cex.yaxis)

    # Setting the default parameters for legend plotting
    if (is.null(legend.x)) {
      legend.x <- xmax * 0.9
    }
    if (is.null(legend.y)) {
      legend.y <- ymax
    }

    legend(x         = legend.x,
           y         = legend.y,
           legend    = rev(row.names(exposure)),
           density   = density,
           angle     = barplot.angle,
           xpd       = NA,
           fill      = legend.col,
           x.intersp = 0.3,
           y.intersp = 1,
           bty       = "n",
           border    = "white",
           cex       = cex.legend,
           title     = "Signature")
    
    # Now add sample names, rotated to hopefully fit,
    if (isTRUE(plot.sample.names)) {
      if (length(mp) < 50) {
        size.adj <- 0.75
      } else if (length(mp) < 80) {
        size.adj <- 0.65
      } else if (length(mp) < 100) {
        size.adj <- 0.4
      } else if (length(mp) < 120) {
        size.adj <- 0.4
      } else if (length(mp) < 150) {
        size.adj <- 0.3
      } else {
        size.adj <- 0.3
      }
      cnames <- colnames(exposure)
      cnames <- sub("_____.*", "", cnames)
      if (is.null(cex.xaxis)) {
        mtext(cnames, side = 1, line = 0.38, at = mp, las = 2, cex = size.adj)
      } else {
        mtext(cnames, side = 1, line = 0.38, at = mp, las = 2, cex = cex.xaxis)
      }
    }

    invisible(list(plot.success = TRUE, mp.coordinates = mp))
  }

#' Plot exposures in multiple plots each with a manageable number of samples
#'
#' @inheritParams PlotExposureInternal
#'
#' @param samples.per.line Number of samples to show in each plot.
#'
#' @param xlim,ylim Limits for the x and y axis. If \code{NULL}(default), the
#'   function tries to do something reasonable.
#'
#' @param legend.x,legend.y The x and y co-ordinates to be used to position the
#'   legend.
#'
#' @return An \strong{invisible} list whose first element is a logic value
#'   indicating whether the plot is successful. The second element is a numeric
#'   vector giving the coordinates of all the bar midpoints drawn, useful for
#'   adding to the graph.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "mSigAct")
#' exposure <- ReadExposure(file)
#' PlotExposure(exposure[, 1:30])
PlotExposure <- function(exposure,
                         samples.per.line   = 30,
                         plot.proportion    = FALSE,
                         xlim               = NULL,
                         ylim               = NULL,
                         legend.x           = NULL,
                         legend.y           = NULL,
                         cex.legend         = 0.9,
                         cex.yaxis          = 1,
                         cex.xaxis          = NULL,
                         plot.sample.names  = TRUE,
                         yaxis.labels       = NULL,
                         ...
) {
  n.sample <- ncol(exposure)
  num.ranges <- n.sample %/% samples.per.line
  size.of.last.range <- n.sample %% samples.per.line
  
  if (num.ranges > 0) {
    if (size.of.last.range > 0) {
      padding.len <- samples.per.line - size.of.last.range
      padding <- matrix(0,nrow = nrow(exposure), ncol = padding.len)
      # The column names starting with lots of underscore
      # will not be plotted in the final output.
      colnames(padding) <- paste("_____", 1:ncol(padding), sep = "_")
      exposure <- cbind(exposure, padding)
      starts <- 0:num.ranges * samples.per.line + 1
    } else {
      starts <- 0:(num.ranges - 1) *samples.per.line + 1
    }
    ends   <- starts + samples.per.line - 1
  } else {
    starts <- 1
    ends <- n.sample
  }
  
  for (i in 1:length(starts)) {
    list <- PlotExposureInternal(exposure[ , starts[i]:ends[i], drop = FALSE],
                                 plot.proportion   = plot.proportion,
                                 xlim              = xlim,
                                 ylim              = ylim,
                                 legend.x          = legend.x,
                                 legend.y          = legend.y,
                                 cex.legend        = cex.legend,
                                 cex.yaxis         = cex.yaxis,
                                 cex.xaxis         = cex.xaxis,
                                 plot.sample.names = plot.sample.names,
                                 yaxis.labels      = yaxis.labels,
                                 ...               = ...)
  }
  invisible(list(plot.success = TRUE, mp.coordinates = list$mp.coordinates))
}

#' Plot exposures in multiple plots each with a manageable number of samples to PDF
#'
#' @inheritParams PlotExposureInternal
#'
#' @param file The name of the PDF file to be produced.
#'
#' @param mfrow A vector of the form \code{c(nr, nc)}.
#' Subsequent figures will be drawn in an \code{nr}-by-\code{nc}
#' array on the device by rows.
#'
#' @param mar A numerical vector of the form \code{c(bottom,
#' left, top, right)} which gives the number of lines of margin to be
#' specified on the four sides of the plot.
#'
#' @param oma A vector of the form \code{c(bottom, left, top,
#' right)} giving the size of the outer margins in lines of text.
#'
#' @param samples.per.line Number of samples to show in each plot.
#'
#' @param xlim,ylim Limits for the x and y axis. If \code{NULL}(default), the
#'   function tries to do something reasonable.
#'
#' @param legend.x,legend.y The x and y co-ordinates to be used to position the
#'   legend.
#'   
#' @param width,height The width and height of the graphics region in inches. 
#'
#' @return An \strong{invisible} list whose first element is a logic value
#'   indicating whether the plot is successful. The second element is a numeric
#'   vector giving the coordinates of all the bar midpoints drawn, useful for
#'   adding to the graph.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata",
#'                     "synthetic.exposure.csv",
#'                     package = "mSigAct")
#' exposure <- ReadExposure(file)
#' PlotExposureToPdf(exposure, file = file.path(tempdir(), "exposure.pdf"))
PlotExposureToPdf <- function(exposure,
                              file,
                              mfrow             = c(2, 1),
                              mar               = c(6, 4, 3, 2),
                              oma               = c(3, 2, 0, 2),
                              samples.per.line  = 30,
                              plot.proportion   = FALSE,
                              xlim              = NULL,
                              ylim              = NULL,
                              legend.x          = NULL,
                              legend.y          = NULL,
                              cex.legend        = 0.9,
                              cex.yaxis         = 1,
                              cex.xaxis         = NULL,
                              plot.sample.names = TRUE,
                              yaxis.labels      = NULL,
                              width             = 8.2677,
                              height            = 11.6929,
                              ...
) {
  # Setting the width and length for A4 size plotting
  grDevices::pdf(file, width = width, height = height, onefile = TRUE)

  opar <- par(mfrow = mfrow, mar = mar, oma = oma)
  on.exit(par(opar))

  list <- PlotExposure(exposure          = exposure,
                       samples.per.line  = samples.per.line,
                       plot.proportion   = plot.proportion,
                       xlim              = xlim,
                       ylim              = ylim,
                       legend.x          = legend.x,
                       legend.y          = legend.y,
                       cex.legend        = cex.legend,
                       cex.yaxis         = cex.yaxis,
                       cex.xaxis         = cex.xaxis,
                       plot.sample.names = plot.sample.names,
                       yaxis.labels      = yaxis.labels,
                       ...               = ...)

  grDevices::dev.off()
  invisible(list(plot.success = TRUE, mp.coordinates = list$mp.coordinates))
}
