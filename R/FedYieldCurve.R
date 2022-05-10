#' Federal Reserve Interest Rates
#'
#' A dataset containing the interest rates of the Federal Reserve, from January
#' 1982 up to April 2022. The interest rates are market yields on United States
#' Treasury securities with constant maturity (CMT). The maturities contained
#' in this dataset are the 3 months, 6 months, 1 year, 2 years, 3 years,
#' 5 years, 7 years, and 10 years maturities. Each interest rate is quoted on
#' investment basis, and are reported monthly.
#'
#' @format A data frame with 484 rows and 9 variables:
#' \describe{
#'   \item{Month}{The month of the quoted interest rates, format is yyyy-mm-dd}
#'   \item{M3}{Market yield on U.S. Treasury securities at 3-month constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{M6}{Market yield on U.S. Treasury securities at 6-month constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{Y1}{Market yield on U.S. Treasury securities at 1-year constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{Y2}{Market yield on U.S. Treasury securities at 2-year constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{Y3}{Market yield on U.S. Treasury securities at 3-year constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{Y5}{Market yield on U.S. Treasury securities at 5-year constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{Y7}{Market yield on U.S. Treasury securities at 7-year constant
#'     maturity, quoted on investment basis, in percent per year.}
#'   \item{Y10}{Market yield on U.S. Treasury securities at 10-year constant
#'     maturity, quoted on investment basis, in percent per year.}
#' }
#' @source \url{https://www.federalreserve.gov/datadownload/Build.aspx?rel=H15}
"FedYieldCurve"
