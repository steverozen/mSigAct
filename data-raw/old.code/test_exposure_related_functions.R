context("Exposure related functions")

test_that("Read and write exposure", {
  file <- "testdata/tiny.exposure.csv"
  x <- ReadExposure(file)
  tfile <- tempfile()
  WriteExposure(x, tfile)
  reread.x <- ReadExposure(tfile)
  expect_equal(x, reread.x)

  x <- ReadExposure(file, check.names = FALSE)
  WriteExposure(x, tfile)
  x2 <- ReadExposure(tfile, check.names = FALSE)
  expect_equal(x, x2)

  file2 <- "testdata/tiny.exposure.dup.csv"
  expect_error(x <- ReadExposure(file2, check = FALSE))
})

test_that("PlotExposureInternal function", {
  exposure <- ReadExposure("testdata/Liver-HCC.exposure.csv")
  par(mar = c(6.9, 3.9, 1, 1))
  out <-
    PlotExposureInternal(SortExposure(exposure),
                         xlim = c(0, ncol(exposure) * 2.5),
                         cex.legend = 0.5, main = "test", cex.main = 0.8,
                         cex.xaxis = 0.6)
  expect_equal(out$plot.success, TRUE)

  # Only plot the first 30 samples
  out1 <- PlotExposureInternal(SortExposure(exposure[, 1:30]),
                               plot.proportion = TRUE, cex.legend = 0.5,
                               main = "test1", cex.main = 0.8)
  expect_equal(out1$plot.success, TRUE)

  # Only plot one sample
  out2 <- PlotExposureInternal(SortExposure(exposure[, 1, drop = FALSE]),
                               xlim = c(0, 5), legend.x = 1.5,
                               cex.legend = 0.4)
  expect_equal(out2$plot.success, TRUE)

  # Only plot selected samples
  out3 <-
    PlotExposureInternal(SortExposure(exposure[, 3:6]), cex.legend = 0.5)
  expect_equal(out3$plot.success, TRUE)

  # We can change the position of legend relative to the bar plot by changing
  # xlim and legend.x
  out4 <-
    PlotExposureInternal(SortExposure(exposure[ , 3:6]), xlim = c(0, 25),
                         cex.legend = 0.5, legend.x = 5.5)
  expect_equal(out4$plot.success, TRUE)

  # Plot exposure proportions rather than counts
  out5 <-
    PlotExposureInternal(SortExposure(exposure[ , 3:6]), xlim = c(0, 25),
                         plot.proportion = TRUE, col = c("red", "blue"),
                         cex.legend = 0.5, legend.x = 5.5)
  expect_equal(out5$plot.success, TRUE)
  
  # Change y axis labels
  out6 <-
    PlotExposureInternal(sqrt(SortExposure(exposure)),
                         xlim = c(0, ncol(exposure) * 2.5),
                         cex.legend = 0.5, main = "test", cex.main = 0.8,
                         yaxis.labels = round(seq(0, max(colSums(exposure)), 
                                                  max(colSums(exposure))/4)))
  expect_equal(out6$plot.success, TRUE)
  
  graphics.off()
})

test_that("PlotExposure function", {
  exposure <- ReadExposure("testdata/Liver-HCC.exposure.csv")

  old.par <-
    par(mfcol = c(2, 1), mar = c(2, 3.9, 3.9, 2), oma = c(4, 0, 0, 0))
  on.exit(par(old.par))
  out <- PlotExposure(exposure = SortExposure(exposure[, 1:43]),
                      main = "test", cex.main = 0.8, cex.legend = 0.3,
                      cex.xaxis = 0.6)
  expect_equal(out$plot.success, TRUE)

  # Only plot the first 30 samples
  par(old.par)
  out1 <- PlotExposure(exposure = SortExposure(exposure[, 1:10]),  # Test a trick edge case
                       main = "test1", cex.main = 0.8, cex.legend = 0.45,
                       plot.proportion = TRUE)
  expect_equal(out1$plot.success, TRUE)

  # Only plot one sample
  out2 <- PlotExposure(exposure = SortExposure(exposure[, 1, drop = FALSE]),
                       samples.per.line = 4, cex.legend = 0.45,
                       legend.x = 1.5)
  expect_equal(out2$plot.success, TRUE)

  # Only plot selected samples
  out3 <- PlotExposure(exposure = SortExposure(exposure[, 3:6]),
                       samples.per.line = 4, cex.legend = 0.45,
                       legend.x = 4.5)
  expect_equal(out3$plot.success, TRUE)

  # Plot exposure proportions rather than counts
  par(mfcol = c(2, 1), mar = c(2, 3.9, 1.9, 0), oma = c(2, 0, 0, 0))
  out4 <- PlotExposure(exposure = SortExposure(exposure[, 1:43 ]),
                       plot.proportion = TRUE, cex.legend = 0.2)
  expect_equal(out4$plot.success, TRUE)

  par(old.par)

  # Only plot selected samples with selected signatures
  out5 <- PlotExposure(exposure = SortExposure(exposure[1:4, 3:6]),
                       samples.per.line = 4, plot.proportion = TRUE,
                       col = c("red", "blue"), cex.legend = 0.45,
                       legend.x = 4.5)
  expect_equal(out5$plot.success, TRUE)

  out6 <- PlotExposure(exposure = SortExposure(exposure[, 3:6]),
                       samples.per.line = 4, plot.proportion = TRUE,
                       col = c("red", "blue"), cex.legend = 0.45,
                       legend.x = 4.5)
  expect_equal(out6$plot.success, TRUE)
  
  
  # Change y axis labels
  out7 <-
    PlotExposureInternal(sqrt(SortExposure(exposure)),
                         xlim = c(0, ncol(exposure) * 2.5),
                         cex.legend = 0.5, main = "test", cex.main = 0.8,
                         yaxis.labels = round(seq(0, max(colSums(exposure)), 
                                                  max(colSums(exposure))/4)))
  expect_equal(out7$plot.success, TRUE)
  graphics.off()
})

test_that("PlotExposureToPdf function", {
  exposure <- ReadExposure("testdata/Liver-HCC.exposure.csv")

  out <- PlotExposureToPdf(exposure, file = file.path(tempdir(), "test.pdf"))
  expect_equal(out$plot.success, TRUE)
  
  # Plot the mutations in sqrt scale
  out0 <- PlotExposureToPdf(exposure = sqrt(SortExposure(exposure)), 
                            file = file.path(tempdir(), "test0.pdf"),
                            ylab = "Number of mutations (square root)",
                            cex.xaxis = 0.9)
  expect_equal(out0$plot.success, TRUE)
  
  # Plot the mutations in sqrt scale and change y axis labels
  out0.1 <- PlotExposureToPdf(exposure = sqrt(SortExposure(exposure)), 
                            file = file.path(tempdir(), "test0.1.pdf"),
                            ylab = "Number of mutations (square root)",
                            ylim = c(0, max(colSums(sqrt(exposure)))),
                            yaxis.labels = round(seq(0, max(colSums(exposure)), 
                                                     max(colSums(exposure))/4)))
  expect_equal(out0.1$plot.success, TRUE)

  # Set ylim for the plots
  out1 <- PlotExposureToPdf(exposure = SortExposure(exposure),
                            file = file.path(tempdir(), "test1.pdf"),
                            ylim = c(0, max(colSums(exposure))))
  expect_equal(out1$plot.success, TRUE)

  # Only plot the first 30 samples
  out2 <- PlotExposureToPdf(exposure = SortExposure(exposure[, 1:30]),
                            file = file.path(tempdir(), "test2.pdf"))  # Test a trick edge case
  expect_equal(out2$plot.success, TRUE)

  # Only plot selected samples, we can change mfrow, cex.legend and legend.x
  # to make the plot look nicer.
  out3 <- PlotExposureToPdf(exposure = SortExposure(exposure[, 1:10]),
                            file = file.path(tempdir(), "test3.pdf"),
                            mfrow = c(2, 2), cex.legend = 1.2,
                            legend.x = 5.5,
                            samples.per.line = 5)
  expect_equal(out3$plot.success, TRUE)

  # Only plot selected samples with selected signatures
  out4 <- PlotExposureToPdf(exposure = exposure[1:5, 1:15],
                            file = file.path(tempdir(), "test4.pdf"),
                            mfrow = c(3, 3), cex.legend = 1.38,
                            legend.x = 5.5, ylim = c(0, 3800),
                            legend.y = 3300,
                            samples.per.line = 5)
  expect_equal(out4$plot.success, TRUE)

  # Plot exposure proportions rather than counts
  out5 <- PlotExposureToPdf(exposure = SortExposure(exposure),
                    file = file.path(tempdir(), "test5.pdf"),
                    plot.proportion = TRUE)
  expect_equal(out5$plot.success, TRUE)

  # Only plot selected samples with selected signatures
  out6 <- PlotExposureToPdf(exposure = exposure[1:5, 1:15],
                            file = file.path(tempdir(), "test6.pdf"),
                            mfrow = c(3, 3), cex.legend = 1.38,
                            legend.x = 5.5,
                            col = c("red", "blue", "orange", "green", "purple"),
                            samples.per.line = 5, plot.proportion = TRUE)
  expect_equal(out6$plot.success, TRUE)
  dev.off()
  unlink(file.path(tempdir(), paste0("test", c("", 1:6), ".pdf")))
  unlink("Rplots.pdf")
})
