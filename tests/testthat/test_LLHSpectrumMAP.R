context("LLHSpectrumMAP")

test_that("MAPAssignActivity for ID Catalog", {
  catalog <- ICAMS::ReadCatalog(file = "testdata/PCAWG7-Prost-AdenoCA-ten-samples.csv")
  sample.index <- 1
  catID <- catalog[, sample.index, drop = FALSE]
  ID.sigs <- ICAMS::ReadCatalog(file = "testdata/COSMIC-v3-genome-ID-sigs.csv")
  mutation.type <- "ID"
  cancer.type <- "Prost-AdenoCA"
  sigs.prop <- ExposureProportions(mutation.type = mutation.type,
                                   cancer.type = cancer.type)
  sigs <- ID.sigs[, names(sigs.prop), drop = FALSE]
  
  expected.counts <- c(10.473178824347, 8.8355441775243, 5.95720789740363,  
                       3.00307169694877, 1.46064236540156, 1.05919407783509, 
                       11.1645059582945,  9.06512953022888, 10.4947514740634, 
                       9.36910156382, 6.55963438116599,  12.917766549936, 
                       0.13767780662152, 0.175373832589374, 0.310391779088189,  
                       0.264672648437922, 0.320776356031523, 1.77872534449583, 
                       1.51951577508639,  2.72961552882838, 3.17249860489714, 
                       2.28518843406716, 2.06004293871377,  54.7082205129126, 
                       2.66129698055793, 2.50278837030748, 0.749946929222408, 
                       0.205710377708171, 0.142371857801692, 1.61711382470492, 
                       1.57394947119105,  1.47516352222981, 0.234256824468529, 
                       0.120286924958027, 0.437861967833275,  0.626936302611306, 
                       0.988189346983376, 0.414376060912531, 0.0416418056306033,  
                       0.292100188326646, 0.34192646103419, 0.12577885483542, 
                       6.13710437727828,  0.449842312584933, 0.0856488982219786, 
                       0.0243964009683029, 0.0134145712934847,  0.0621721734524383, 
                       0.949717507958856, 2.00633758443542, 0.95147116433723,  
                       0.154124866365715, 0.211045628971386, 1.36957591893657, 
                       0.544662507842371,  1.25707904053465, 0.223138469280451, 
                       0.0180847475807297, 0.0131884542917751,  0.109338544666039, 
                       0.413921682702962, 1.36290286609563, 0.0902353192869763,  
                       0.0159621477793368, 0.0574437812676172, 0.337864025494402, 
                       1.02274060604814,  4.3875368377243, 0.164601726942931, 
                       0.0800501936759003, 0.0709279160827568,  0.0959952615434431, 
                       2.52109834229767, 1.27534039166227, 1.1701591751858,  
                       0.806013974569732, 0.560267728798564, 0.339141895066014, 
                       6.23322800688953,  3.06579321790787, 1.46975218542343, 
                       0.538193616225794, 0.966361879353789 )
  
  model <- c("ID1", "ID2", "ID3", "ID5", "ID8", "ID9", "ID10")
  loglh.MAP <- LLHSpectrumMAP(spectrum = catID,
                              expected.counts = expected.counts,
                              nbinom.size = 5,
                              model = model,
                              likelihood.dist = "neg.binom",
                              sigs.presence.prop = sigs.prop)
  expect_equal(loglh.MAP, -124.6515, tolerance = 1e-3)
})