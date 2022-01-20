# 2.1.3.9008
* Fixed a bug in internal function `MAPAssignActivityInternal` when there are 
more signatures in argument `sigs.presence.prop` than that in `sigs`.

# 2.1.3.9007
* Fixed a bug in internal function `NumFromId` when there are no numerical parts from
the signature ids.

# 2.1.3.9006
* Fixed a bug in internal function `TestAltSolutions` when sorting the all
tested solutions using sparse assignment.

# 2.1.3.9005
* Updated message information when `use.sparse.assign = TRUE` in functions
`MAPAssignActivity1` and `MAPAssignActivity`.

# 2.1.3.9004
* Added a new internal function `MergeListOfExposures`.

# 2.1.3.9003
* Added a new internal function `CalculateDistance`.

# 2.1.3.9002
* Added a new internal function `RemoveZeroActivitySig`.

# 2.1.3.9001
* Added a new internal function `RemoveZeroMutationSample` for `AddSigActivity`
to remove samples which have zero counts either in spectra or exposure matrix.

# 2.1.3.9000
* Fixed a bug in parsing argument `drop.low.mut.samples` in function `MAPAssignActivity` and `SparseAssignActivity`.

# 2.1.3
* Added new exported functions `ReadExposure`, `WriteExposure`, `SortExposure`, 
`PlotExposure`, `PlotExposureToPdf` and necessary tests and test data.

* Added new internal function `NumFromId` and removed `ICAMSxtra` from the dependency packages.

* Enabled argument `target.sig.index` to also take the signature id as input in
function `SignaturePresenceTest`. Changed argument `m.opts` to use return value
from `DefaultManyOpts` by default. Changed the default value for argument `mc.cores` to 2.

* Added `inst/extdata` folder to be used in running the examples.

* Added `ipc` as one of the Suggests packages in DESCRIPTION.

# 2.1.2
* Enabled argument `target.sig.index` to also take the signature id as input in
function `SignaturePresenceTest1`. Changed argument `m.opts` to use return value
from `DefaultManyOpts` by default.

* Updated exported function `RareSignatures` to include SBS91 and SBS94
as rare signatures.

* Added documentation of return value for internal function `SignaturePresenceTest1`.
Changed the names of first two returned elements: `with` -> `loglh.with` and
`without` -> `loglh.without`. Removed two NULL elements `everything.else.with` and 
`everything.else.without`.

* Added documentation of return value for exported function `SignaturePresenceTest`.

* Moved several functions `LLHSpectrumNegBinomDebug`, `ObjFnBinomMaxLHMustRound`,
`ObjFnBinomMaxLHNoRoundOK` and `OptimizeExposureQPBootstrap` to data-raw old.code
folder.

* Made several functions `AddSigActivity`, `g_ineq_for_ObjFnBinomMaxLH2`, 
`g_ineq_for_ObjFnMultinomMaxLH`, `LLHSpectrumMAP`, `LLHSpectrumNegBinom`,
`MAPAssignActivity1`, `ObjFnBinomMaxLHRound`, `ObjFnMultinomMaxLH`,
`OneMAPAssignTest `, `OptimizeExposure`, `OptimizeExposureQP`, `PCAWGMAPTest`,
`ShowSigActivity`, `SignaturePresenceTest1` **internal** functions.

* Removed `rlang` from Imports field in DESCRIPTION file.

* Removed argument `max.presence.proportion` from internal functions `OneMAPAssignTest`
and `PCAWGMAPTest`.

# 2.1.1.9016
* Fixed a bug in function `ReconstructSpectrum` to allow reconstructing multiple
spectra

* Added examples to the exported functions

# 2.1.1.9015
* Updated function `SparseAssignActivity` to use `MAPAssignActivity` with 
argument `use.sparse.assign` set to be TRUE

# 2.1.1.9014
* Changed name of returned element from `time.for.MAP.assign` to
`time.for.assignment` in function `MAPAssignActivity`

# 2.1.1.9013
* Fixed a bug in function `AddSigActivity` when some samples in the spectra have
zero mutations.

# 2.1.1.9012
* Added new arguments `likelihood.dist` and `use.sparse.assign` in functions `AddSigActivity1` 
and `AddSigActivity1`.

# 2.1.1.9011
* Added new internal function `DropLowMutationSamples`

* Added new argument `drop.low.mut.samples` for functions `MAPAssignActivity`
and `MAPAssignActivity1`. Samples with SBS total mutations less than 100, DBS or
ID total mutations less than 25 will be excluded from the analysis


# 2.1.1.9010
* Used Benjamini-Hochberg false discovery rate to adjust p values to determine the `alt.solutions`
that are statistically as good as the `proposed.assignment` from `MAPAssignActivity` and `MAPAssignActivity1`

# 2.1.1.9009
* Updated functions `MAPAssignActivity` and `MAPAssignActivity1` to return extra element `alt.solutions` that are statistically as good as the `proposed.assignment` that can
plausibly reconstruct the spectra

# 2.1.1.9008
* Removed columns `sig.indices` and `removed.sig.names` in the returned all tested table
from function `MAPAssignActivity` and `MAPAssignActivity1`

* Added new internal functions `GetAllTestedTables` and `GetTimeForMAPAssign`

* Changed element `results.details` to `all.tested` in the returned value of
function `MAPAssignActivity`

* Added new element `time.for.MAP.assign` in the returned value of
function `MAPAssignActivity`

# 2.1.1.9007

* Added a new argument `use.sparse.assign` in functions `MAPAssignActivity`,
`MAPAssignActivity1` and `MAPAssignActivityInternal`

# 2.1.1.9006

* Changed the column name from `neg.log.likelihood` to `log.likelihood` in the
output of functions `DistanceMeasures` and `DistanceMeasuresSparse`

# 2.1.1.9005

* Added example, new arguments and return values to function
`SparseAssignActivity`. See its documentation for more details

* Changed the default tracing behavior to only show the information of sample
name, df and number of subsets to remove

# 2.1.1.9004

* Added a new argument `likelihood.dist` to functions `DefaultManyOpts`,
`LLHSpectrumMAP` and `DistanceMeasures` so that user can specify the probability
distribution to calculate the likelihood

* Added a new element `likelihood.dist` in the returned value of function
`DefaultManyOpts`

# 2.1.1.9003

* Updated documentation for the return values of `DefaultManyOpts` to include three additional 
elements `global_eval_f`, `local_eval_f` and `local_eval_g_ineq`

# 2.1.1.9002

* Added new argument `signatures` to internal function `DistanceMeasures` to enable calculating
distance measures using quadratic programming assignment

# 2.1.1.9001

* Updated function `ShowSigActivity` to create the `output.dir` recursively if it
does not exist

# 2.1.1.9000

* Updated function `P.of.M` to give more informative error message when there are signatures used in the model but don't have presence
proportions.

# 2.1.1

* Added back argument `max.subsets` in function `MAPAssignActivity1` and `MAPAssignActivity`

* Added an additional argument `base.filename` in function `ShowSigActivity`

# 2.1.0.9009

* Updated function `MAPAssignActivity1` to exclude signatures in the
`proposed.assignment` if the mutations get assigned to the signatures is zero
after rounding

# 2.1.0.9008

* Updated long running tests by setting `Sys.setenv("MSIGACT_TEST_LENGTH" = "long")`

# 2.1.0.9007

* Adapted to the new version of PCAWG7
* Added dependency on ICAMSxtra version

# 2.1.0.9006
* Round the `proposed.assignment` and `proposed.reconstruction` from output running
`MAPAssignActivity1 `, `MAPAssignActivity`.

* Added new internal function `AddAttributes`

# 2.1.0.9005
* Added examples to exported function
`MAPAssignActivity1 `, `MAPAssignActivity` and `ExposureProportions`

* Added new exported function `CancerTypes`

* Added new exported 
# 2.1.0.9004
* Decreased the number of returned elements from exported functions
`MAPAssignActivity1 ` and `MAPAssignActivity` and renamed them to make them more
informative

# 2.1.0.9003
* Removed parameter `max.subsets` from exported functions
`MAPAssignActivity1 ` and `MAPAssignActivity`

# 2.1.0.9002
* Removed parameter `max.presence.proportion` from exported functions
`MAPAssignActivity1 ` and `MAPAssignActivity`

* Fixed a bug in internal function `RunMAPOnOneSample` to replace `::` to `.` in the file 
path (otherwise error will occur on Windows OS)

# 2.1.0.9001
* Fixed a bug in internal function `PlotSigActivityToPdf` to replace `::` to `.` in the file 
path otherwise `grDevices::pdf` will throw an error

# 2.1.0.9000
* Created new exported function `AddSigActivity` to add contributing signature
activity information for multiple spectra

* Created new exported function `ShowSigActivity` to show signature activity
from the output generated by `AddSigActivity`

# 2.0.10.9003

* Added two internal functions `GetDistanceInfo` and `GetExposureInfo` to
retrieve information from output returned by `MAPAssignActivity`

# 2.0.10.9002

* Added signature etiologies information to the reconstructed spectrum plot when calling `MAPAssignActivity`

# 2.0.10.9001

* Added `LLHSpectrumMAP` and tests

# 2.0.10

* Added `MAPAssignActivity` and tests

# 2.0.9.9000

* Added sample information to show to the user when calling `MAPAssignActivity1`

# 2.0.8.9005

* Round the non integers in spectrum to integers and give a warning when calling
`MAPAssignActivity1`

# 2.0.8.9004

* Fixed typo in documentation for argument `must.include` in function `ExposureProportions`

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
