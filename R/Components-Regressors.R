#' Construct the System Matrices of a Explanatory Variables Component
#'
#' Constructs the system matrices of a explanatory variables component.
#'
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceFit
#' @inheritParams StateSpaceEval
#' @inheritParams LocalLevel
#' @inheritParams Cholesky
#'
#' @return
#' A list containing the system matrices.
#'
#' @noRd
AddVar <- function(p = 1,
                   addvar_list,
                   fixed_part = TRUE,
                   update_part = TRUE,
                   param = NULL,
                   format_addvar = diag(0, 1, 1),
                   decompositions = TRUE) {

  # Check if the addvar list has p elements
  if (length(addvar_list) != p) {
    stop(
      paste(
        "The number of elements in `addvar_list` must",
        "be equal to the number of dependent variables."
      )
    )
  }

  # Variable that is used to check if Q should be a matrix of 0s
  diag_addvar <- sum(diag(format_addvar) != 0)

  # Number of coefficients
  k <- sum(sapply(addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[2]}}))

  # Initialising the list to return
  result <- list()

  if (fixed_part) {

    # Number of observations
    N <- max(sapply(addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[1]}}))

    # Z = a p x k x N array
    Ztemp <- do.call(BlockMatrix, addvar_list)
    index <- 1:N
    Ztemp2 <- NULL
    null_mat <- matrix(0, N, k)
    for (i in seq_along(addvar_list)) {
      if (is.null(addvar_list[[i]])) {
        Ztemp2 <- rbind(Ztemp2, null_mat)
      } else {
        Ztemp2 <- rbind(Ztemp2, Ztemp[index,])
        index <- index + N
      }
    }
    result$Z <- array(
      sapply(
        seq(N),
        function(i) {Ztemp2[seq(i, p * N, N),, drop = FALSE]},
        simplify = "array"
      ), dim = c(p, k, N)
    )

    # T matrix = a k x k diagonal matrix, containing ones on the diagonal
    result$Tmat <- diag(1, k, k)

    # R matrix = a k x k diagonal matrix, containing ones on the diagonal
    result$R <- diag(1, k, k)

    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, k, 1)
    result$P_inf <- diag(1, k, k)
    result$P_star <- matrix(0, k, k)

    # Check if Q depends on parameters
    if (diag_addvar == 0) {
      result$Q <- matrix(0, k, k)
    }
  }

  if (update_part & diag_addvar > 0) {

    # Check whether the number of rows of format_addvar is valid
    if (dim(format_addvar)[1] != k) {
      stop(
        paste(
          "The number of rows of `format_addvar` for the explanatory variables",
          "must be equal to the number of explanatory variables."
        )
      )
    }

    # Using Cholesky function to get a valid variance - covariance matrix
    # for the Q matrix
    Q <- Cholesky(
      param = param,
      format = format_addvar,
      decompositions = decompositions
    )

    # Check what to return
    if (decompositions) {
      result$Q <- Q$cov_mat
      result$loading_matrix <- Q$loading_matrix
      result$diagonal_matrix <- Q$diagonal_matrix
      result$correlation_matrix <- Q$correlation_matrix
      result$stdev_matrix <- Q$stdev_matrix
    } else {
      result$Q <- Q
    }
  }

  return(result)
}

#' Construct the System Matrices of a Explanatory Variables
#' in the Level Component
#'
#' Constructs the system matrices of a explanatory variables
#' in the level component.
#'
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceFit
#' @inheritParams StateSpaceEval
#' @inheritParams LocalLevel
#' @inheritParams Cholesky
#'
#' @return
#' A list containing the system matrices.
#'
#' @noRd
LevelAddVar <- function(p = 1,
                        exclude_level = NULL,
                        level_addvar_list,
                        fixed_part = TRUE,
                        update_part = TRUE,
                        param = rep(1, p),
                        format_level = diag(1, p, p),
                        format_level_addvar = diag(0, 1, 1),
                        decompositions = TRUE) {

  # Check if the level_addvar list has p elements
  if (length(level_addvar_list) != p) {
    stop(
      paste(
        "The number of elements in `level_addvar_list`",
        "must be equal to the number of dependent variables."
      )
    )
  }

  # The number of dependent variables that should get a local level
  n_level <- p - length(exclude_level)

  # Variable that is used to check if Q_level should be a matrix of 0s
  diag_level <- sum(diag(format_level) != 0)

  # Variable that is used to check if Q_level_addvar should be a matrix of 0s
  diag_level_addvar <- sum(diag(format_level_addvar) != 0)

  # Initialising the list to return
  result <- list()

  if (fixed_part) {

    # Number of coefficients
    k <- sum(sapply(level_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[2]}}))

    # Number of observations
    N <- max(sapply(level_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[1]}}))

    # Z matrix = matrix of p rows and n_level + k columns
    #   First n_level columns containing one 1, 0s elsewhere
    #   Last k columns containing zeroes
    Z <- cbind(diag(1, p, p), matrix(0, p, k))
    if (n_level < p) {
      Z <- Z[, -exclude_level, drop = FALSE]
    }
    result$Z <- Z

    # Tmat = a (n_level + k) x (n_level + k) x N array
    Ttemp <- do.call(BlockMatrix, level_addvar_list)
    index <- 1:N
    Ttemp2 <- NULL
    null_mat <- matrix(0, N, k)
    for (i in seq_along(level_addvar_list)) {
      if (is.null(level_addvar_list[[i]])) {
        if (!(i %in% exclude_level)) {
          Ttemp2 <- rbind(Ttemp2, null_mat)
        }
      } else {
        if (i %in% exclude_level) {
          stop(
            paste0(
              "The ", i, "th dependent variable was specified to be excluded ",
              "from getting a local level, ",
              "but did have explanatory variables specified for the local level."
            )
          )
        } else {
          Ttemp2 <- rbind(Ttemp2, Ttemp[index,])
          index <- index + N
        }
      }
    }
    result$Tmat <- sapply(
      seq(N),
      function(i) {rbind(
        cbind(
          diag(1, n_level, n_level),
          Ttemp2[seq(i, n_level * N, N),, drop = FALSE]
        ),
        cbind(
          matrix(0, k, n_level),
          diag(1, k, k)
        )
      )},
      simplify = "array"
    )

    # R matrix = a (n_level + k) x (n_level + k) diagonal matrix containing 1s
    result$R <- diag(1, n_level + k, n_level + k)

    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, n_level + k, 1)
    result$P_inf <- diag(1, n_level + k, n_level + k)
    result$P_star <- matrix(0, n_level + k, n_level + k)

    # Check if Q_level depends on parameters
    if (diag_level == 0) {
      result$Q_level <- matrix(0, n_level, n_level)
    }

    # Check if Q_level_addvar depends on parameters
    if (diag_level_addvar == 0) {
      result$Q_level_addvar <- matrix(0, k, k)
    }
  }

  if (update_part & (diag_level + diag_level_addvar) > 0) {

    if (diag_level > 0) {

      # Check whether the number of rows of format_level is valid
      if (dim(format_level)[1] != n_level) {
        stop(
          paste(
            "The number of rows of `format_level` for the level +",
            "explanatory variables in the level component must be equal",
            "to the number of dependent variables minus the number of",
            "excluded dependent variables."
          )
        )
      }

      # Which parameters should be used for the Q matrix corresponding to the level
      split <- sum(format_level != 0 & lower.tri(format_level, diag = TRUE))

      # Using Cholesky function to get a valid variance - covariance matrix
      # for the Q matrix for the level
      Q_level <- Cholesky(
        param = param[1:split],
        format = format_level,
        decompositions = decompositions
      )

      # Check what to return
      if (decompositions) {
        result$Q_level <- Q_level$cov_mat
        result$loading_matrix_level <- Q_level$loading_matrix
        result$diagonal_matrix_level <- Q_level$diagonal_matrix
        result$correlation_matrix_level <- Q_level$correlation_matrix
        result$stdev_matrix_level <- Q_level$stdev_matrix
      } else {
        result$Q_level <- Q_level
      }

    } else {
      split <- 0
    }

    if (diag_level_addvar > 0) {

      # Check whether the number of rows of format_level_addvar is valid
      if (dim(format_level_addvar)[1] != k) {
        stop(
          paste(
            "The number of rows of `format_level_addvar` for the level +",
            "explanatory variables in the level component",
            "must be equal to the number of explanatory variables."
          )
        )
      }

      # Using Cholesky function to get a valid variance - covariance matrix
      # for the Q matrix for the level_addvar
      Q_level_addvar <- Cholesky(
        param = param[(split + 1) : length(param)],
        format = format_level_addvar,
        decompositions = decompositions
      )

      # Check what to return
      if (decompositions) {
        result$Q_level_addvar <- Q_level_addvar$cov_mat
        result$loading_matrix_level_addvar <- Q_level_addvar$loading_matrix
        result$diagonal_matrix_level_addvar <- Q_level_addvar$diagonal_matrix
        result$correlation_matrix_level_addvar <- Q_level_addvar$correlation_matrix
        result$stdev_matrix_level_addvar <- Q_level_addvar$stdev_matrix
      } else {
        result$Q_level_addvar <- Q_level_addvar
      }
    }
  }

  return(result)
}

#' Construct the System Matrices of a Explanatory Variables in the
#' Level + Slope Component
#'
#' Constructs the system matrices of a explanatory variables in the
#' level + slope component.
#'
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceFit
#' @inheritParams StateSpaceEval
#' @inheritParams LocalLevel
#' @inheritParams Cholesky
#'
#' @return
#' A list containing the system matrices.
#'
#' @noRd
SlopeAddVar <- function(p = 1,
                        exclude_level = NULL,
                        exclude_slope = NULL,
                        level_addvar_list,
                        fixed_part = TRUE,
                        update_part = TRUE,
                        param = rep(1, 2 * p),
                        format_level = diag(1, p, p),
                        format_slope = diag(1, p, p),
                        format_level_addvar = diag(0, 1, 1),
                        decompositions = TRUE) {

  # Check if the level_addvar list has p elements
  if (length(level_addvar_list) != p) {
    stop(
      paste(
        "The number of elements in `level_addvar_list` must be equal",
        "to the number of dependent variables."
      )
    )
  }

  # The number of dependent variables that should get a local level
  n_level <- p - length(exclude_level)

  # The number of local levels that should get a slope
  n_slope <- p - length(exclude_slope)

  # Variable that is used to check if Q_level should be a matrix of 0s
  diag_level <- sum(diag(format_level) != 0)

  # Variable that is used to check if Q_slope should be a matrix of 0s
  diag_slope <- sum(diag(format_slope) != 0)

  # Variable that is used to check if Q_level_addvar should be a matrix of 0s
  diag_level_addvar <- sum(diag(format_level_addvar) != 0)

  # Initialising the list to return
  result <- list()

  if (fixed_part) {

    # Number of coefficients
    k <- sum(
      sapply(
        level_addvar_list,
        function(X) {if (is.null(X)) {0} else {dim(X)[2]}}
      )
    )

    # Number of observations
    N <- max(
      sapply(
        level_addvar_list,
        function(X) {if (is.null(X)) {0} else {dim(X)[1]}}
      )
    )

    # Z matrix = matrix of p rows and n_level + n_slope + k columns
    #   First n_level columns containing one 1, 0s elsewhere
    #   Second n_slope columns containing zeroes
    #   Last k columns containing zeroes
    Z <- cbind(diag(1, p, p), matrix(0, p, n_slope), matrix(0, p, k))
    if (n_level < p) {
      Z <- Z[, -exclude_level, drop = FALSE]
    }
    result$Z <- Z

    # Tmat = a (n_level + n_slope + k) x (n_level + n_slope + k) x N array
    Ttemp <- do.call(BlockMatrix, level_addvar_list)
    index <- 1:N
    Ttemp2 <- NULL
    null_mat <- matrix(0, N, k)
    for (i in seq_along(level_addvar_list)) {
      if (is.null(level_addvar_list[[i]])) {
        if (!(i %in% exclude_level)) {
          Ttemp2 <- rbind(Ttemp2, null_mat)
        }
      } else {
        if (i %in% exclude_level) {
          stop(
            paste0(
              "The ", i, "th dependent variable was specified to be excluded ",
              "from getting a local level, ",
              "but did have explanatory variables specified for the local level."
            )
          )
        } else {
          Ttemp2 <- rbind(Ttemp2, Ttemp[index,])
          index <- index + N
        }
      }
    }
    Ttemp3 <- rbind(
      cbind(
        diag(1, p, p),
        diag(1, p, p)
      ),
      cbind(
        matrix(0, p, p),
        diag(1, p, p)
      )
    )
    if ((n_level + n_slope) < (2 * p)) {
      exclude <- c(exclude_level, p + exclude_slope)
      Ttemp3 <- Ttemp3[-exclude, -exclude, drop = FALSE]
    }
    result$Tmat <- sapply(
      seq(N),
      function(i) {rbind(
        cbind(
          Ttemp3,
          rbind(
            Ttemp2[seq(i, n_level * N, N),, drop = FALSE],
            matrix(0, n_slope, k)
          )
        ),
        cbind(
          matrix(0, k, n_level + n_slope),
          diag(1, k, k)
        )
      )},
      simplify = "array"
    )

    # R matrix = a (n_level + n_slope + k) x (n_level + n_slope + k)
    #            diagonal matrix containing 1s
    result$R <- diag(1, n_level + n_slope + k, n_level + n_slope + k)

    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, n_level + n_slope + k, 1)
    result$P_inf <- diag(1, n_level + n_slope + k, n_level + n_slope + k)
    result$P_star <- matrix(0, n_level + n_slope + k, n_level + n_slope + k)

    # Check if Q_level depends on parameters
    if (diag_level == 0) {
      result$Q_level <- matrix(0, n_level, n_level)
    }

    # Check if Q_slope depends on parameters
    if (diag_slope == 0) {
      result$Q_slope <- matrix(0, n_slope, n_slope)
    }

    # Check if Q_level_addvar depends on parameters
    if (diag_level_addvar == 0) {
      result$Q_level_addvar <- matrix(0, k, k)
    }
  }

  if (update_part & (diag_level + diag_slope + diag_level_addvar) > 0) {

    if (diag_level > 0) {

      # Check whether the number of rows of format_level is valid
      if (dim(format_level)[1] != n_level) {
        stop(
          paste(
            "The number of rows of `format_level` for the level",
            "+ slope + explanatory variables in the level component must be",
            "equal to the number of dependent variables minus",
            "the number of excluded dependent variables."
          )
        )
      }

      # Which parameters should be used for the Q matrix corresponding to the level
      split <- sum(format_level != 0 & lower.tri(format_level, diag = TRUE))

      # Using Cholesky function to get a valid variance - covariance matrix
      # for the Q matrix for the level
      Q_level <- Cholesky(
        param = param[1:split],
        format = format_level,
        decompositions = decompositions
      )

      # Check what to return
      if (decompositions) {
        result$Q_level <- Q_level$cov_mat
        result$loading_matrix_level <- Q_level$loading_matrix
        result$diagonal_matrix_level <- Q_level$diagonal_matrix
        result$correlation_matrix_level <- Q_level$correlation_matrix
        result$stdev_matrix_level <- Q_level$stdev_matrix
      } else {
        result$Q_level <- Q_level
      }

    } else {
      split <- 0
    }

    if (diag_slope > 0) {

      # Check whether the number of rows of format_slope is valid
      if (dim(format_slope)[1] != n_slope) {
        stop(
          paste(
            "The number of rows of `format_slope` for the level",
            "+ slope + explanatory variables in the level component must be",
            "equal to the number of local levels minus",
            "the number of excluded local levels."
          )
        )
      }

      # Which parameters should be used for the Q matrix corresponding to the slope
      split2 <- split + sum(format_slope != 0 & lower.tri(format_slope, diag = TRUE))

      # Using Cholesky function to get a valid variance - covariance matrix
      # for the Q matrix for the slope
      Q_slope <- Cholesky(
        param = param[(split + 1):split2],
        format = format_slope,
        decompositions = decompositions
      )

      # Check what to return
      if (decompositions) {
        result$Q_slope <- Q_slope$cov_mat
        result$loading_matrix_slope <- Q_slope$loading_matrix
        result$diagonal_matrix_slope <- Q_slope$diagonal_matrix
        result$correlation_matrix_slope <- Q_slope$correlation_matrix
        result$stdev_matrix_slope <- Q_slope$stdev_matrix
      } else {
        result$Q_slope <- Q_slope
      }

    } else {
      split2 <- 0
    }

    if (diag_level_addvar > 0) {

      # Check whether the number of rows of format_level_addvar is valid
      if (dim(format_level_addvar)[1] != k) {
        stop(
          paste(
            "The number of rows of `format_level_addvar` for the level",
            "+ slope + explanatory variables in the level component must be",
            "equal to the number of explanatory variables."
          )
        )
      }

      # Using Cholesky function to get a valid variance - covariance matrix
      # for the Q matrix for the level_addvar
      Q_level_addvar <- Cholesky(
        param = param[(split + split2 + 1) : length(param)],
        format = format_level_addvar,
        decompositions = decompositions
      )

      # Check what to return
      if (decompositions) {
        result$Q_level_addvar <- Q_level_addvar$cov_mat
        result$loading_matrix_level_addvar <- Q_level_addvar$loading_matrix
        result$diagonal_matrix_level_addvar <- Q_level_addvar$diagonal_matrix
        result$correlation_matrix_level_addvar <- Q_level_addvar$correlation_matrix
        result$stdev_matrix_level_addvar <- Q_level_addvar$stdev_matrix
      } else {
        result$Q_level_addvar <- Q_level_addvar
      }
    }
  }

  return(result)
}
