#' statespacer: A package for state space modelling in R.
#'
#' The statespacer package provides functions that make estimating models in
#' State Space form a breeze. This package implements state-of-the-art
#' algorithms developed by various time series practitioners such as J. Durbin
#' and S.J. Koopman. Details about the algorithms can be found in their book,
#' "Time Series Analysis by State Space Methods".
#'
#' @section State Space Components:
#' This package supports numerous state space components:
#' * The Local Level
#' * The Local Level + Slope
#' * Smoothing Splines
#' * Trigonometric Seasonality, BSM
#' * Cycles
#' * Business Cycles
#' * Explanatory Variables
#' * Explanatory Variables with time-varying coefficients
#' * Explanatory Variables in the Local Level
#' * Explanatory Variables in the Local Level + Slope
#' * ARIMA
#' * SARIMA
#'
#' These components can be used for both univariate, as multivariate models.
#' The components can be combined in order to get more extensive models.
#' Moreover, the user can control the format of the variance - covariance
#' matrices of each of the components. This way, one could specify the
#' components to be deterministic instead of stochastic. In the multivariate
#' case, one could set some correlations to be 0 if desired. Also, the rank of
#' the variance - covariance matrix could be set such that commonalities in the
#' components are estimated (common levels, common slopes, etc.).
#'
#' @section Fitting Procedure:
#' The package employs a univariate treatment, and an exact initialisation for
#' diffuse elements, to estimate the state parameters.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#'
#' @references
#' \insertRef{durbin2012time}{statespacer}
#'
#' @docType package
#' @keywords internal
"_PACKAGE"
