# mSigAct 2.0.5

* Lots new tracing information in SparseAssign (controlled
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
