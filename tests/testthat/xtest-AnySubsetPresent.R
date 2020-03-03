if (FALSE) {
  any.retval <- TestAny1("SBS17a", 1)
}

if (FALSE) {
  spt.retval <- TestSignaturePresenceTest("SBS17a", 1)
  stopifnot(all.equal(eso.17a$`Eso-AdenoCA::SP111062`$chisq.p, 0.1019716, tolerance = 1e-7))
  
}
