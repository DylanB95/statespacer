#' Construct a Valid Variance - Covariance Matrix
#'
#' Constructs a valid variance - covariance matrix by using the Cholesky LDL
#' decomposition.
#'
#' @param param Vector containing the parameters used to construct the
#'   variance - covariance matrix.
#' @param format Matrix representing the format for the Loading matrix L
#'   and Diagonal matrix D. The lower triangular part of the format is used
#'   as the format for the Loading matrix L. The diagonal of the format is
#'   used as the format for the Diagonal matrix D. Must be a matrix.
#' @param decompositions Boolean indicating whether the loading and diagonal
#'   matrix of the Cholesky decomposition, and the correlation matrix and
#'   standard deviations should be returned.
#'
#' @details
#' `format` is used to specify which elements of the loading and diagonal
#' matrix should be non-zero. The elements of `param` are then distributed
#' along the non-zero elements of the loading and diagonal matrix.
#' The parameters for the diagonal matrix are transformed using `exp(2 * x)`.
#'
#' @return
#' A valid variance - covariance matrix.
#' If `decompositions = TRUE` then it returns a list containing:
#' * `cov_mat`: The variance - covariance matrix.
#' * `loading_matrix`: The loading matrix of the Cholesky decomposition.
#' * `diagonal_matrix`: The diagonal matrix of the Cholesky decomposition.
#' * `correlation_matrix`: Matrix containing the correlations.
#' * `stdev_matrix`: Matrix containing the standard deviations on the diagonal.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#'
#' @examples
#' format <- diag(1, 2, 2)
#' format[2, 1] <- 1
#' Cholesky(param = c(2, 4, 1), format = format, decompositions = TRUE)
#' @export
Cholesky <- function(param = NULL, format = NULL, decompositions = TRUE) {

  # Check if at least one of param and format is specified
  if (is.null(param) & is.null(format)) {
    stop(
      paste(
        "Both `param` and `format` were not specified.",
        "Must specify at least one of them."
      ),
      call. = FALSE
    )
  }

  # Number of parameters that are specified
  param <- param[!is.na(param)]
  n_par <- length(param)

  # If no format is specified
  if (is.null(format)) {

    # Calculate the dimension using the number of parameters that are specified
    dimension <- (-1 + sqrt(1 + 8 * n_par)) / 2

    # if calculated dimension is not an integer, return an error message
    if ((dimension %% 1) != 0) {
      stop(
        paste(
          "Number of parameters supplied result in non-integer dimensions",
          "of the variance - covariance matrix,",
          "specify a correct number of parameters.",
          "Hint: The dimension is calculated as (-1 + sqrt(1 + 8 * n_par))/2."
        ),
        call. = FALSE
      )
    }

    # LDL Cholesky algorithm for covariance matrix
    #   L matrix: Putting ones on the diagonal
    #             and the parameters on the lower triangle of the matrix
    #   D matrix: Diagonal matrix containing the magnitudes
    chol_L <- diag(1, dimension, dimension)
    chol_L[lower.tri(chol_L)] <- param[(dimension + 1):n_par]
    chol_D <- diag(exp(2 * param[1:dimension]), dimension, dimension)

    # If format is specified
  } else {

    # Dimension of format
    format_dim <- dim(format)

    # Number of columns must not exceed the number of rows
    if (format_dim[1] < format_dim[2]) {
      stop(
        paste(
          "Number of columns of `format` must be less than",
          "or equal to the number of rows."
        ),
        call. = FALSE
      )
    }

    # Only the lower triangle (including diagonal) of the format is needed,
    # because the LDL decomposition uses a lower triangular Loading matrix
    format[upper.tri(format)] <- 0

    # Non zero elements of the format
    format_n0 <- format != 0

    # Non zero elements of the lower triangular part of the format
    format_n0_lt <- format_n0 & lower.tri(format)

    # Non zero elements of the diagonal of the format
    format_n0_diag <- diag(format_n0)

    # Number of parameters required for the matrix (lower + diagonal)
    lower <- sum(format_n0_lt)
    diagonal <- sum(format_n0_diag)

    # Check if magnitudes that are set to 0, have non-zero coefficients
    # specified in the loading matrix L
    if (diagonal < sum(colSums(format_n0) != 0)) {
      stop(
        paste(
          "The specified format is not valid.",
          "Columns of which the diagonal element is zero,",
          "should be set to zero."
        ),
        call. = FALSE
      )
    }

    # If format is a matrix of zeroes, return the format itself
    if (diagonal == 0) {
      return(format)
    }

    # Check if more parameters are specified than needed
    if ((lower + diagonal) < n_par) {
      stop("Too many parameters supplied.", call. = FALSE)
    }

    # Check if not enough parameters are specified
    if ((lower + diagonal) > n_par) {
      stop("Not enough parameters supplied.", call. = FALSE)
    }

    # Initialising L matrix with specified dimensions
    chol_L <- diag(1, format_dim[1], format_dim[2])

    # Putting the parameters into the right places according to format
    if (lower > 0) {
      chol_L[format_n0_lt] <- param[(diagonal + 1):(diagonal + lower)]
    }

    # Constructing D matrix
    # Note: exp(-Inf) = 0. Using a magnitude equal to zero,
    #       where the diagonal entry of the format equals zero
    param_diag <- rep(-Inf, format_dim[2])
    param_diag[format_n0_diag] <- param[1:diagonal]
    chol_D <- diag(exp(2 * param_diag), format_dim[2], format_dim[2])
  }

  # LDL Cholesky algorithm, resulting in a valid Variance - Covariance matrix
  cov_mat <- chol_L %*% chol_D %*% t(chol_L)

  # Returning the result
  if (decompositions) {

    # Diagonal matrix containing the standard deviations
    stdev_matrix <- diag(sqrt(diag(cov_mat)), dim(cov_mat)[1], dim(cov_mat)[2])

    # Inverse of stdev_matrix
    stdev_inv <- stdev_matrix
    diag(stdev_inv)[diag(stdev_inv) > 0] <-
      1 / diag(stdev_inv)[diag(stdev_inv) > 0]

    # Correlation matrix
    correlation_matrix <- stdev_inv %*% cov_mat %*% stdev_inv

    result <- list(
      cov_mat = cov_mat,
      loading_matrix = chol_L,
      diagonal_matrix = chol_D,
      correlation_matrix = correlation_matrix,
      stdev_matrix = stdev_matrix
    )
    return(result)
  } else {
    return(cov_mat)
  }
}
