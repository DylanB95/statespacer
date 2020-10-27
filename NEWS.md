# statespacer 0.3.9000

## Extra functionality

* Introduces the simulation smoother for drawing random samples conditional on the observed data. See `?SimSmoother` for details.

## Performance improvements

* Extraction of components now done a bit more efficiently.

## Bug fixes

* Fixed incrementation of initialisation_steps during initialisation when the dependent variable was missing.

# statespacer 0.3.0

## Performance improvements

* Kalman Filter and Smoother now fully written in c++.

## Extra functionality

* Computation of diagnostics now optional using `diagnostics = TRUE`.
* Added full system matrices to the statespacer object.
* Added predicted, filtered, and forecasted `P_star` and `P_inf`.

## Bug fixes

* Fixed bug in predicting using explanatory variables.
* Fixed edge case in adding explanatory variables. The bug occurred when there was only 1 explanatory variable and 1 dependent variable supplied, causing the Z system matrix to be a vector instead of an array.

# statespacer 0.2.1

## Bug fixes

* Patch for macOS. On macOS, assignment of nested lists is handled a bit differently than on other platforms. For instance, if `first <- list()`, then assigning `first$second[[["third"]] <- 1` returns on Windows (and other platforms) a list named `first`, that contains another list named `second`, that contains a named element `third` equal to an unnamed length 1 numerical vector. On macOS though, it returns a list named `first`, that contains a named element `second` equal to a named (`third`) length 1 numerical vector. So `second` is a list on the other platforms, while being a named numerical vector on macOS. This caused a bug on macOS while computing standard errors.

# statespacer 0.2.0

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

  * Easier to find proper initial values, reducing time spent on trial and error by the user.

## Extra functionality

* Printing progress now optional using `verbose = TRUE`.

* Computation of standard_errors now optional using `standard_errors = TRUE`.

# statespacer 0.1.0

* Initial release to CRAN!
