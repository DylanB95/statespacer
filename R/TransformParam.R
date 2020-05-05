#' Transform the optimisation parameters to the parameters of a State Space Model
#'
#' Transforms the optimisation parameters to the parameters
#' of a State Space Model.
#'
#' @inheritParams StateSpaceFit
#' @inheritParams StateSpaceEval
#' @inheritParams GetSysMat
#'
#' @return
#' A vector containing the transformed parameters.
#'
#' @noRd
TransformParam <- function(param = NULL,
                           p,
                           H_format = NULL,
                           local_level_ind = FALSE,
                           slope_ind = FALSE,
                           BSM_vec = NULL,
                           cycle_ind = FALSE,
                           addvar_list = NULL,
                           level_addvar_list = NULL,
                           arima_list = NULL,
                           sarima_list = NULL,
                           self_spec_list = NULL,
                           exclude_level = NULL,
                           exclude_slope = NULL,
                           exclude_BSM_list = lapply(BSM_vec, FUN = function(x) 0),
                           exclude_cycle_list = list(0),
                           exclude_arima_list = lapply(arima_list, FUN = function(x) 0),
                           exclude_sarima_list = lapply(sarima_list, FUN = function(x) 0),
                           damping_factor_ind = rep(TRUE, length(exclude_cycle_list)),
                           format_level = NULL,
                           format_slope = NULL,
                           format_BSM_list = lapply(BSM_vec, FUN = function(x) NULL),
                           format_cycle_list = lapply(exclude_cycle_list, FUN = function(x) NULL),
                           format_addvar = NULL,
                           format_level_addvar = NULL) {

  # Construct the system matrices
  sys_mat <- GetSysMat(p = p,
                       param = param,
                       update_part = TRUE,
                       add_residuals = FALSE,
                       H_format = H_format,
                       local_level_ind = local_level_ind,
                       slope_ind = slope_ind,
                       BSM_vec = BSM_vec,
                       cycle_ind = cycle_ind,
                       addvar_list = addvar_list,
                       level_addvar_list = level_addvar_list,
                       arima_list = arima_list,
                       sarima_list = sarima_list,
                       self_spec_list = self_spec_list,
                       exclude_level = exclude_level,
                       exclude_slope = exclude_slope,
                       exclude_BSM_list = exclude_BSM_list,
                       exclude_cycle_list = exclude_cycle_list,
                       exclude_arima_list = exclude_arima_list,
                       exclude_sarima_list = exclude_sarima_list,
                       damping_factor_ind = damping_factor_ind,
                       format_level = format_level,
                       format_slope = format_slope,
                       format_BSM_list = format_BSM_list,
                       format_cycle_list = format_cycle_list,
                       format_addvar = format_addvar,
                       format_level_addvar = format_level_addvar
  )

  # Initialise vector of transformed parameters
  result <- c()

  # H
  if (!is.null(sys_mat[["H"]])) {
    for (i in sys_mat[["H"]]) {
      result <- c(result, i)
    }
  }

  # Q
  if (!is.null(sys_mat[["Q"]])) {
    for (i in sys_mat[["Q"]]) {
      result <- c(result, i)
    }
  }

  # Q_loading_matrix
  if (!is.null(sys_mat$Q_loading_matrix)) {
    for (i in sys_mat$Q_loading_matrix) {
      result <- c(result, i)
    }
  }

  # Q_diagonal_matrix
  if (!is.null(sys_mat$Q_diagonal_matrix)) {
    for (i in sys_mat$Q_diagonal_matrix) {
      result <- c(result, i)
    }
  }

  # Q_correlation_matrix
  if (!is.null(sys_mat$Q_correlation_matrix)) {
    for (i in sys_mat$Q_correlation_matrix) {
      result <- c(result, i)
    }
  }

  # Q_stdev_matrix
  if (!is.null(sys_mat$Q_stdev_matrix)) {
    for (i in sys_mat$Q_stdev_matrix) {
      result <- c(result, i)
    }
  }

  # lambda
  if (!is.null(sys_mat$lambda)) {
    for (i in sys_mat$lambda) {
      result <- c(result, i)
    }
  }

  # rho
  if (!is.null(sys_mat$rho)) {
    for (i in sys_mat$rho) {
      result <- c(result, i)
    }
  }

  # AR
  if (!is.null(sys_mat$AR)) {
    for (i in sys_mat$AR) {
      result <- c(result, i)
    }
  }

  # MA
  if (!is.null(sys_mat$MA)) {
    for (i in sys_mat$MA) {
      result <- c(result, i)
    }
  }

  # SAR
  if (!is.null(sys_mat$SAR)) {
    for (j in sys_mat$SAR) {
      for (i in j) {
        result <- c(result, i)
      }
    }
  }

  # SMA
  if (!is.null(sys_mat$SMA)) {
    for (j in sys_mat$SMA) {
      for (i in j) {
        result <- c(result, i)
      }
    }
  }

  # Self Specified transformed parameters
  if (!is.null(sys_mat$self_spec)) {
    result <- c(result, sys_mat$self_spec)
  }

  return(result)
}

#' Add structure to the parameters of a State Space Model
#'
#' Adds structure to the parameters of a State Space Model.
#'
#' @param sys_mat The list as returned by GetSysMat.
#' @inheritParams StateSpaceEval
#'
#' @return
#' A list containing the structured parameters of a State Space Model.
#'
#' @noRd
StructParam <- function(param = NULL,
                        sys_mat) {

  # Initialise list to return
  result <- list()

  # H
  if (!is.null(sys_mat[["H"]])) {
    H <- list()
    for (i in sys_mat[["H"]]) {
      indices <- 1:length(i)
      if (is.matrix(i)) {
        param_mat <- matrix(param[indices], dim(i)[1], dim(i)[2])
      } else {
        param_mat <- array(param[indices], dim = dim(i))
      }
      param <- param[-indices]
      H <- c(H, list(param_mat))
    }
    names(H) <- names(sys_mat[["H"]])
    result[["H"]] <- H
  }

  # Q
  if (!is.null(sys_mat[["Q"]])) {
    Q <- list()
    for (i in sys_mat[["Q"]]) {
      indices <- 1:length(i)
      if (is.matrix(i)) {
        param_mat <- matrix(param[indices], dim(i)[1], dim(i)[2])
      } else {
        param_mat <- array(param[indices], dim = dim(i))
      }
      param <- param[-indices]
      Q <- c(Q, list(param_mat))
    }
    names(Q) <- names(sys_mat[["Q"]])
    result[["Q"]] <- Q
  }

  # Q_loading_matrix
  if (!is.null(sys_mat$Q_loading_matrix)) {
    Q_loading_matrix <- list()
    for (i in sys_mat$Q_loading_matrix) {
      indices <- 1:length(i)
      param_mat <- matrix(param[indices], dim(i)[1], dim(i)[2])
      param <- param[-indices]
      Q_loading_matrix <- c(Q_loading_matrix, list(param_mat))
    }
    names(Q_loading_matrix) <- names(sys_mat$Q_loading_matrix)
    result$Q_loading_matrix <- Q_loading_matrix
  }

  # Q_diagonal_matrix
  if (!is.null(sys_mat$Q_diagonal_matrix)) {
    Q_diagonal_matrix <- list()
    for (i in sys_mat$Q_diagonal_matrix) {
      indices <- 1:length(i)
      param_mat <- matrix(param[indices], dim(i)[1], dim(i)[2])
      param <- param[-indices]
      Q_diagonal_matrix <- c(Q_diagonal_matrix, list(param_mat))
    }
    names(Q_diagonal_matrix) <- names(sys_mat$Q_diagonal_matrix)
    result$Q_diagonal_matrix <- Q_diagonal_matrix
  }

  # Q_correlation_matrix
  if (!is.null(sys_mat$Q_correlation_matrix)) {
    Q_correlation_matrix <- list()
    for (i in sys_mat$Q_correlation_matrix) {
      indices <- 1:length(i)
      param_mat <- matrix(param[indices], dim(i)[1], dim(i)[2])
      param <- param[-indices]
      Q_correlation_matrix <- c(Q_correlation_matrix, list(param_mat))
    }
    names(Q_correlation_matrix) <- names(sys_mat$Q_correlation_matrix)
    result$Q_correlation_matrix <- Q_correlation_matrix
  }

  # Q_stdev_matrix
  if (!is.null(sys_mat$Q_stdev_matrix)) {
    Q_stdev_matrix <- list()
    for (i in sys_mat$Q_stdev_matrix) {
      indices <- 1:length(i)
      param_mat <- matrix(param[indices], dim(i)[1], dim(i)[2])
      param <- param[-indices]
      Q_stdev_matrix <- c(Q_stdev_matrix, list(param_mat))
    }
    names(Q_stdev_matrix) <- names(sys_mat$Q_stdev_matrix)
    result$Q_stdev_matrix <- Q_stdev_matrix
  }

  # lambda
  if (!is.null(sys_mat$lambda)) {
    lambda <- list()
    for (i in sys_mat$lambda) {
      lambda <- c(lambda, list(param[1]))
      param <- param[-1]
    }
    names(lambda) <- names(sys_mat$lambda)
    result$lambda <- lambda
  }

  # rho
  if (!is.null(sys_mat$rho)) {
    rho <- list()
    for (i in sys_mat$rho) {
      rho <- c(rho, list(param[1]))
      param <- param[-1]
    }
    names(rho) <- names(sys_mat$rho)
    result$rho <- rho
  }

  # AR
  if (!is.null(sys_mat$AR)) {
    AR <- list()
    for (i in sys_mat$AR) {
      indices <- 1:length(i)
      if (is.vector(i)) {
        param_ar <- param[indices]
      } else {
        param_ar <- array(param[indices], dim = dim(i))
      }
      param <- param[-indices]
      AR <- c(AR, list(param_ar))
    }
    names(AR) <- names(sys_mat$AR)
    result$AR <- AR
  }

  # MA
  if (!is.null(sys_mat$MA)) {
    MA <- list()
    for (i in sys_mat$MA) {
      indices <- 1:length(i)
      if (is.vector(i)) {
        param_ma <- param[indices]
      } else {
        param_ma <- array(param[indices], dim = dim(i))
      }
      param <- param[-indices]
      MA <- c(MA, list(param_ma))
    }
    names(MA) <- names(sys_mat$MA)
    result$MA <- MA
  }

  # SAR
  if (!is.null(sys_mat$SAR)) {
    SAR <- list()
    for (j in sys_mat$SAR) {
      temp <- list()
      for (i in j) {
        indices <- 1:length(i)
        if (is.vector(i)) {
          param_sar <- param[indices]
        } else {
          param_sar <- array(param[indices], dim = dim(i))
        }
        param <- param[-indices]
        temp <- c(temp, list(param_sar))
      }
      names(temp) <- names(j)
      SAR <- c(SAR, list(temp))
    }
    names(SAR) <- names(sys_mat$SAR)
    result$SAR <- SAR
  }

  # SMA
  if (!is.null(sys_mat$SMA)) {
    SMA <- list()
    for (j in sys_mat$SMA) {
      temp <- list()
      for (i in j) {
        indices <- 1:length(i)
        if (is.vector(i)) {
          param_sma <- param[indices]
        } else {
          param_sma <- array(param[indices], dim = dim(i))
        }
        param <- param[-indices]
        temp <- c(temp, list(param_sma))
      }
      names(temp) <- names(j)
      SMA <- c(SMA, list(temp))
    }
    names(SMA) <- names(sys_mat$SMA)
    result$SMA <- SMA
  }

  # Self Specified transformed parameters
  if (!is.null(sys_mat$self_spec)) {
    indices <- 1:length(sys_mat$self_spec)
    result$self_spec <- param[indices]
    param <- param[-indices]
  }

  # Check if all parameters have been used
  if (length(param) != 0) {
    stop("Not all parameters have been used in assigning standard errors!")
  }

  return(result)
}
