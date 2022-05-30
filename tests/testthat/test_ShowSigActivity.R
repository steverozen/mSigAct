context("ShowSigActivity")

test_that("ShowSigActivity for only one SBS96 catalog", {
  spect <- PCAWG7::spectra$PCAWG$SBS96[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$SBS96))
  exposure <- PCAWG7::exposure$PCAWG$SBS96[, index, drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 1)

  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
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
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
  sigs.prop <- ExposureProportions(mutation.type = "SBS96",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 2)

  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for only one SBS192 catalog", {
  spect <- PCAWG7::spectra$PCAWG$SBS192[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$SBS96))
  exposure <- PCAWG7::exposure$PCAWG$SBS96[, index, drop = FALSE]
  rownames(exposure) <- cosmicsig::SBS96_ID_to_SBS192_ID(rownames(exposure))
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS192
  exposure2 <- ReassignmentQP(spectra = spect, 
                              exposure = exposure, 
                              sigs = sigs)
  sigs.prop <- ExposureProportions(mutation.type = "SBS192",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure2, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 1)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for multiple SBS192 catalog", {
  spect <- PCAWG7::spectra$PCAWG$SBS192[, 1:2, drop = FALSE]
  spect.name <- colnames(spect)
  indices <-
    sapply(spect.name, FUN = grep, x = colnames(PCAWG7::exposure$PCAWG$SBS96))
  exposure <- PCAWG7::exposure$PCAWG$SBS96[, indices, drop = FALSE]
  rownames(exposure) <- cosmicsig::SBS96_ID_to_SBS192_ID(rownames(exposure))
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS192
  exposure2 <- ReassignmentQP(spectra = spect, 
                              exposure = exposure, 
                              sigs = sigs)
  sigs.prop <- ExposureProportions(mutation.type = "SBS192",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure2, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 2)
  
  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for only one DBS78 catalog", {
  spect <- PCAWG7::spectra$PCAWG$DBS78[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$DBS78))
  exposure <- PCAWG7::exposure$PCAWG$DBS78[, index, drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
  sigs.prop <- ExposureProportions(mutation.type = "DBS78",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 1)

  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  base.filename = "Biliary-AdenoCA",
                  output.dir = output.dirs[1])
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
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
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$DBS78
  sigs.prop <- ExposureProportions(mutation.type = "DBS78",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 2)

  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
})

test_that("ShowSigActivity for only one ID catalog", {
  spect <- PCAWG7::spectra$PCAWG$ID[, 1, drop = FALSE]
  spect.name <- colnames(spect)
  index <- grep(spect.name, x = colnames(PCAWG7::exposure$PCAWG$ID))
  exposure <- PCAWG7::exposure$PCAWG$ID[, index, drop = FALSE]
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  sigs.prop <- ExposureProportions(mutation.type = "ID",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 1)

  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
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
  sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
  sigs.prop <- ExposureProportions(mutation.type = "ID",
                                   cancer.type = "Biliary-AdenoCA")
  retval <- AddSigActivity(spect = spect, exposure = exposure, sigs = sigs,
                           sigs.presence.prop = sigs.prop)
  expect_equal(length(retval), 2)

  output.dirs <- file.path(tempdir(), paste0("test", 1:2))
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[1],
                  base.filename = "Biliary-AdenoCA")
  ShowSigActivity(list.of.sig.activity = retval,
                  output.dir = output.dirs[2],
                  base.filename = "Biliary-AdenoCA",
                  plot.all.samples.in.one.pdf = FALSE,
                  plot.exposure.proportion = TRUE)
  sapply(output.dirs, FUN = unlink, recursive = TRUE)
  unlink("Rplots.pdf")
})
