context("ShowSigActivity")

test_that("ShowSigActivity for only one SBS96 catalog", {
  spect <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$SBS96))
  exposure <- PCAWG7::exposure$PCAWG$SBS96[, index, drop = FALSE]
  sigs <- PCAWG7::COSMIC.v3.1$signature$genome$SBS96
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs)
  expect_equal(length(retval), 1)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[2],
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for multiple SBS96 catalog", {
  spect <- PCAWG7::spectra$PCAWG$SBS96[, 1:2, drop = FALSE]
  spect.name <- colnames(spect)
  indices <- 
    sapply(spect.name, FUN = grep, x = colnames(PCAWG7::exposure$PCAWG$SBS96))
  exposure <- PCAWG7::exposure$PCAWG$SBS96[, indices, drop = FALSE]
  sigs <- PCAWG7::COSMIC.v3.1$signature$genome$SBS96
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs)
  expect_equal(length(retval), 2)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[2],
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for only one DBS78 catalog", {
  spect <- PCAWG7::spectra$PCAWG$DBS78[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$DBS78))
  exposure <- PCAWG7::exposure$PCAWG$DBS78[, index, drop = FALSE]
  sigs <- PCAWG7::COSMIC.v3.1$signature$genome$DBS78
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs)
  expect_equal(length(retval), 1)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[2],
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for multiple DBS78 catalog", {
  spect <- PCAWG7::spectra$PCAWG$DBS78[, 1:2, drop = FALSE]
  spect.name <- colnames(spect)
  indices <- 
    sapply(spect.name, FUN = grep, x = colnames(PCAWG7::exposure$PCAWG$DBS78))
  exposure <- PCAWG7::exposure$PCAWG$DBS78[, indices, drop = FALSE]
  sigs <- PCAWG7::COSMIC.v3.1$signature$genome$DBS78
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs)
  expect_equal(length(retval), 2)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[2],
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for only one ID catalog", {
  spect <- PCAWG7::spectra$PCAWG$ID[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$ID))
  exposure <- PCAWG7::exposure$PCAWG$ID[, index, drop = FALSE]
  sigs <- PCAWG7::COSMIC.v3.1$signature$genome$ID
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs)
  expect_equal(length(retval), 1)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[2],
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for multiple ID catalog", {
  spect <- PCAWG7::spectra$PCAWG$ID[, 1:2, drop = FALSE]
  spect.name <- colnames(spect)
  indices <- 
    sapply(spect.name, FUN = grep, x = colnames(PCAWG7::exposure$PCAWG$ID))
  exposure <- PCAWG7::exposure$PCAWG$ID[, indices, drop = FALSE]
  sigs <- PCAWG7::COSMIC.v3.1$signature$genome$ID
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs)
  expect_equal(length(retval), 2)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval, 
                  output.dir = output.dirs[2],
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
  unlink("Rplots.pdf")
})