oo <- PCAWG7::spectra$other.genome$SBS96
ii <- grep("ML_31_T_01", colnames(oo))
mel <- oo[ , ii, drop = F]
sgs <- PCAWG7::signature$genome$SBS96
universe <- ExposureProportions("SBS96", "Skin-Melanoma")
universe <- universe[-(which(names(universe) == "SBS58"))]
new.sgs <- sgs[ , names(universe)]

res <- MAPAssignActivity(mel, new.sgs,
                         sigs.presence.prop = universe,
                         output.dir = "tmp2")


test.sgs <- rownames(res$proposed.assignment)

ss <- sgs[ , test.sgs]
qp.exp <- OptimizeExposureQP(spectrum = mel, signatures = ss)
map.exp <- res$proposed.assignment
View(cbind(qp.exp, map.exp))

qp.recon <- ReconstructSpectrum(sigs = ss, exp = qp.exp)
map.recon <- res$proposed.reconstruction
cossim(qp.recon, mel)
cossim(map.recon, mel)


qp.ll <- ObjFnBinomMaxLHRound(exp = qp.exp, spectrum = mel, sigs = ss, nbinom.size = 5)
qp.ll

map.ll <- ObjFnBinomMaxLHRound(exp = map.exp, spectrum = mel, sigs = ss, nbinom.size = 5)

### Next time to get per-channel likelihoods call
# stats::dnbinom(
# x = spectrum, mu = expected.counts, size = nbinom.size, log =TRUE)

class(qp.all.channel) <- class(qp.all.channel)[-1]
class(per.channel.map) <- class(per.channel.map)[-1]



exp.times.sigs <- function(exp.vec, sig.matrix) {
  out <- sig.matrix
  for (ii in 1:length(exp.vec)) {
    out[ ,ii ] <- exp.vec[ii] * sig.matrix[ , ii]
  }
  return(out)
}

qp.partial <- exp.times.sigs(qp.exp, ss)
map.partial <- exp.times.sigs(map.exp, ss)

class(qp.partial) <- class(qp.partial)[-1]
class(map.partial) <- class(map.partial)[-1]

comp <- cbind(qp.ll = qp.all.channel[ , 1],
              map.ll = per.channel.map[ ,1],
              diff.ll = (qp.all.channel[ , 1] - per.channel.map[, 1]),
              spect = mel[  , 1],
              qp.recon = qp.recon[ , 1],
              map.recon = map.recon[ , 1],
              qp.partial,
              map.partial)


dmultinom(mel[ , 1], prob = map.recon[ , 1], log = T)
