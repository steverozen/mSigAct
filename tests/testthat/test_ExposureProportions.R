test_that("ExposureProportions 1", {
  rr <- ExposureProportions(mutation.type = "SBS96", "Lung-AdenoCA")
  expect_equal(rr,
               c(SBS1 = 0.947368421052632, SBS2 = 0.578947368421053, 
                 SBS3 = 0.0789473684210526,  SBS4 = 0.552631578947368, 
                 SBS5 = 0.921052631578947, SBS9 = 0.0263157894736842,  
                 SBS13 = 0.578947368421053, SBS17a = 0.0789473684210526, 
                 SBS17b = 0.0789473684210526,  SBS18 = 0.210526315789474, 
                 SBS28 = 0.0263157894736842, SBS40 = 0.315789473684211 ))

})

test_that("ExposureProportions 2", {
  rr <- ExposureProportions(
    mutation.type = "SBS96", 
    "Lung-AdenoCA",
    must.include = "SBS6")
  expect_equal(rr,
               c(SBS1 = 0.947368421052632, SBS2 = 0.578947368421053, 
                 SBS3 = 0.0789473684210526,  SBS4 = 0.552631578947368, 
                 SBS5 = 0.921052631578947, SBS9 = 0.0263157894736842,  
                 SBS13 = 0.578947368421053, SBS17a = 0.0789473684210526, 
                 SBS17b = 0.0789473684210526,  SBS18 = 0.210526315789474, 
                 SBS28 = 0.0263157894736842, SBS40 = 0.315789473684211,
                 SBS6  = 0.1))
  
})

test_that("ExposureProportions 3", {
  rr <- ExposureProportions(
    mutation.type = "SBS96", 
    "Lung-AdenoCA",
    must.include = "SBS6",
    must.include.prop = 0.3)
  expect_equal(rr,
               c(SBS1 = 0.947368421052632, SBS2 = 0.578947368421053, 
                 SBS3 = 0.0789473684210526,  SBS4 = 0.552631578947368, 
                 SBS5 = 0.921052631578947, SBS9 = 0.0263157894736842,  
                 SBS13 = 0.578947368421053, SBS17a = 0.0789473684210526, 
                 SBS17b = 0.0789473684210526,  SBS18 = 0.210526315789474, 
                 SBS28 = 0.0263157894736842, SBS40 = 0.315789473684211,
                 SBS6  = 0.3))
  
})

test_that("ExposureProportions 3", {
  rr <- ExposureProportions(
    mutation.type = "SBS192", 
    "Lung-AdenoCA",
    must.include = c("SBS6", "SBS29"),
    must.include.prop = 0.3)
  expect_equal(rr,
               c(SBS1 = 0.947368421052632, SBS2 = 0.578947368421053, 
                 SBS3 = 0.0789473684210526,  SBS4 = 0.552631578947368, 
                 SBS5 = 0.921052631578947, SBS9 = 0.0263157894736842,  
                 SBS13 = 0.578947368421053, SBS17a = 0.0789473684210526, 
                 SBS17b = 0.0789473684210526,  SBS18 = 0.210526315789474, 
                 SBS28 = 0.0263157894736842, SBS40 = 0.315789473684211,
                 SBS6  = 0.3, SBS29 = 0.3))
  
})

test_that("ExposureProportions 4", {
  rr <- ExposureProportions(
    mutation.type = "SBS192", 
    "Lung-AdenoCA",
    must.include = c("SBS1", "SBS6", "SBS29"),
    must.include.prop = 0.3)
  expect_equal(rr,
               c(SBS1 = 0.947368421052632, SBS2 = 0.578947368421053, 
                 SBS3 = 0.0789473684210526,  SBS4 = 0.552631578947368, 
                 SBS5 = 0.921052631578947, SBS9 = 0.0263157894736842,  
                 SBS13 = 0.578947368421053, SBS17a = 0.0789473684210526, 
                 SBS17b = 0.0789473684210526,  SBS18 = 0.210526315789474, 
                 SBS28 = 0.0263157894736842, SBS40 = 0.315789473684211,
                 SBS6  = 0.3, SBS29 = 0.3))
})

test_that("ExposureProportions 5", {
  rr <- ExposureProportions(mutation.type = "SBS192", 
                            cancer.type   = "Liver-HCC",
                            all.sigs      = PCAWG7::signature$genome$SBS192)
  expect_equal(rr,
               c(c(SBS1 = 0.653374233128834, SBS4 = 0.269938650306748, SBS5 = 1, 
                   SBS6 = 0.00920245398773006, 
                   SBS9 = 0.00613496932515337, 
                   SBS12 = 0.607361963190184, SBS14 = 0.00920245398773006, 
                   SBS16 = 0.119631901840491, SBS17a = 0.0122699386503067, 
                   SBS17b = 0.0184049079754601, SBS18 = 0.00306748466257669, 
                   SBS19 = 0.0184049079754601, SBS22 = 0.0521472392638037, 
                   SBS24 = 0.00306748466257669, SBS26 = 0.00613496932515337, 
                   SBS28 = 0.00306748466257669, SBS30 = 0.00613496932515337, 
                   SBS31 = 0.00306748466257669, SBS35 = 0.0306748466257669,  
                   SBS40 = 0.119631901840491, SBS53 = 0.00306748466257669, 
                   SBS54 = 0.00613496932515337, SBS56 = 0.00306748466257669) ))
  
})

test_that("ExposureProportions 6", {
  rr <- ExposureProportions(mutation.type = "SBS192", 
                            cancer.type   = "Unknown",
                            all.sigs      = PCAWG7::signature$genome$SBS192)
  expect_equal(rr, numeric(0))
})

test_that("ExposureProportions 7", {
    rr <- ExposureProportions(mutation.type = "SBS192", 
                              cancer.type   = "Unknown",
                              all.sigs      = PCAWG7::signature$genome$SBS192,
                              must.include  = "FOO")
  expect_equal(rr, numeric(0), check.attributes = FALSE)
})

test_that("ExposureProportions 8", {
  rr <- ExposureProportions(mutation.type = "SBS192", 
                            cancer.type   = "Unknown",
                            all.sigs      = PCAWG7::signature$genome$SBS192,
                            must.include  = "SBS9")
  expect_equal(rr, c(SBS9 = 0.1))
})


