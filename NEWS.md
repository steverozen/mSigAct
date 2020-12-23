# 2.0.8.9003

* Wrapped MAPAssignActivity1 in a tryCatch and improved some error messages

# 2.0.8.9002

* Can now specify proportions of additional signatures (signatures
  not previously observed in a given cancer type) to ExposureProportions.

* ExposureProportions can no accept cancer.type = "Unknown" or "None"

# 2.0.8.9001

* Changed return value from OptimizeExposureQPBootstrap (not backward compatible)

# 2.0.8

* Cleaned up output of AnySigSubsetPresent
* Added RareSignatures
* Added OptimizeExposureQPBootstrap

# 2.0.7.9001

* Removed the eval_f and eval_g_ineq arguments everywhere; they are now
  in the options list.

# 2.0.6.9006

* Updated the test for AnySigSubsetPresent and added notes on number of sigs in each tumor type

# 2.0.6.9005

* Changed the name of callback function of MAPAssignActivity1 to `progress.monitor`. Added code to update the progress bar when exiting 
the for loop

* Exported function `PossibleArtifacts` and `OptimizeExposureQP`

# 2.0.6.9004

* More performance tests

# 2.0.6.9003

* More performance tests

# 2.0.6.9002

* Added callback fn to MAPAssignActivity1

# 2.0.6.9001

* Added additional elements to return list from MAPAssignActivity1

# 2.0.6

* Created a stable branch, added sigs.prop on OneMapAssignTest and PCAWGMAPTest

# 2.0.5.9011

* Corrected minor bugs in testing MAPAssignActivity

# mSigAct 2.0.5.9010

* Added MAPAssignActivity1 and tests

# mSigAct 2.0.5.9000

* Minor code cleanup

# mSigAct 2.0.5

* Lots new tracing information in SparseAssignActivity1 (controlled
  by the $trace element in argument m.opts; for lots of info
  set this to e.g. 11)
  
* New defaults to keep from spawning too many children from
  SparseAssignActivity1.

# mSigAct 2.0.4.9000

* Multiple changes, mostly backward compatible

# mSigAct 2.0.0.9018

* Had previously failed to export OptimizeExposure.

# mSigAct 2.0.0.9017

* Aliased one.lh.and.exp as OptimizeExposure and made
  OptimizeExposure public.

# mSigAct 2.0.0.9016

* Improved documentation for AnySigSubsetPresent; 
  considering deprecating ObjFnBinomMaxLHMustRound, but need to do more
  tests first -- is it possible that rounding improves attribution?

# mSigAct 2.0.0.9014

* Added AnySigSubsetPresent() and tests.
* Parallelized SparseAssignActivity1().

# mSigAct 2.0.0.0010

* Upgraded functions to test for signature presence and associated tests.

# mSigAct 2.0.0.0009

* Substantially improved separation of background signature, include some tutorial
  examples.

# mSigAct 2.0.0.0005

* Create ObjFnBinomMaxLH2() to provide an option to round the reconstruction only if no mutation
  types have 0 counts. (In this case log likelihoods of observed data may be -inf, and the 
  numerical search might not advance.)

# mSigAct 2.0.0.0004

* Adding background subtraction code.

# mSigAct 2.0.0 

* First version as a package. Not all functionality from the alpha version has been incorporated.

# mSigAct v1.2-alpha

* Version of the code used in Ng et al., 2017, "Aristolochic acids and their derivatives are widely 
  implicated in liver cancers in Taiwan and throughout Asia", Science Translational Medicine 2017
  https://doi.org/10.1126/scitranslmed.aan6446. The supplementary information for
  the paper has additional code that used v1.2-alpha.
