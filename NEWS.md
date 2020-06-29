# statespacer 0.2.0 (in development)

## Breaking changes

* API change by making use of S3 classes and methods:

  * Replaced `StateSpaceFit()` for `statespacer()`.

  * Replaced `StateSpaceEval(...)` for `statespacer(..., fit = FALSE)`.

  * `statespacer()` returns a list object of class `statespacer`.

  * Replaced `StateSpaceForecast()` for the S3 method `predict.statespacer()`.

## Performance improvements

* Major:

  * Calculation of loglikelihood now done using c++. Major performance improvement as the loglikelihood potentially gets calculated a lot of times by the optimisation procedure.

* Minor:
  
  * Making use of crossprod and tcrossprod.

  * Improved efficiency of the computation of standard_errors.
  
  * y_temp outside of LogLikelihood function.

## Extra functionality

* Printing progress now optional using `verbose = TRUE`.

* Computation of standard_errors now optional using `standard_errors = TRUE`.

# statespacer 0.1.0

* Initial release to CRAN!
