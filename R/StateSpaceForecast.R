#' State Space Model Forecasting
#'
#' Produces forecasts using a fitted State Space Model.
#'
#' @param fit A list containing the specifications of the State Space Model, as
#'   returned by \code{\link{StateSpaceFit}} or \code{\link{StateSpaceEval}}.
#' @param addvar_list_fc A list containing the explanatory variables for each
#'   of the dependent variables. The list should contain p (number of dependent
#'   variables) elements. Each element of the list should be a
#'   `forecast_period` x k_p matrix, with k_p being the number of explanatory
#'   variables for the pth dependent variable. If no explanatory variables
#'   should be added for one of the dependent variables, then set the
#'   corresponding element to `NULL`.
#' @param level_addvar_list_fc A list containing the explanatory variables
#'   for each of the dependent variables. The list should contain p
#'   (number of dependent variables) elements. Each element of the list should
#'   be a `forecast_period` x k_p matrix, with k_p being the number of
#'   explanatory variables for the pth dependent variable. If no explanatory
#'   variables should be added for one of the dependent variables, then set
#'   the corresponding element to `NULL`.
#' @param self_spec_list_fc A list containing the specification of the self
#'   specified component. Does not have to be specified if it does not differ
#'   from `self_spec_list` as passed on to `StateSpaceFit()` or
#'   `StateSpaceEval()`. If some system matrices are time-varying then you
#'   should specify this argument. See `StateSpaceFit()` for details about the
#'   format that must be followed for this argument.
#' @param forecast_period Number of time steps to forecast ahead.
#'
#' @return
#' A list containing the forecasts and corresponding uncertainties.
#' In addition, it returns the components of the forecasts, as specified
#' by the State Space model.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#' @references
#' \insertRef{durbin2012time}{statespacer}
#'
#' @examples
#' # Fits a local level model for the Nile data
#' library(datasets)
#' y <- matrix(Nile)
#' fit <- StateSpaceFit(initial = 1, y = y / 100, local_level_ind = TRUE)
#'
#' # Obtain forecasts for 10 steps ahead using the fitted model
#' fc <- StateSpaceForecast(fit, forecast_period = 10)
#'
#' # Plot the forecasts
#' plot(1:10, fc$y_fc * 100, type = 'l')
#'
#' @export
StateSpaceForecast <- function(fit,
                               addvar_list_fc = NULL,
                               level_addvar_list_fc = NULL,
                               self_spec_list_fc = NULL,
                               forecast_period = 1) {

  # Check if specification of addvar_list_fc is in line with the fit object
  if (!is.null(fit$function_call$addvar_list) & is.null(addvar_list_fc)) {
    stop("`addvar_list_fc` must be specified for the forecasting period.")
  }
  if (is.null(fit$function_call$addvar_list) & !is.null(addvar_list_fc)) {
    stop(
      paste(
        "`addvar_list_fc` was specified for the forecasting period,",
        "while explanatory variables were not incorporated in the model."
      )
    )
  }

  # Check if specification of level_addvar_list_fc is in line with the fit object
  if (!is.null(fit$function_call$level_addvar_list) &
      is.null(level_addvar_list_fc)) {
    stop("`level_addvar_list_fc` must be specified for the forecasting period.")
  }
  if (is.null(fit$function_call$level_addvar_list) &
      !is.null(level_addvar_list_fc)) {
    stop(
      paste(
        "`level_addvar_list_fc` was specified for the forecasting period,",
        "while explanatory variables in the level were not",
        "incorporated in the model."
      )
    )
  }

  # Check whether self_spec_list_fc should be used
  if (!is.null(fit$function_call$self_spec_list)) {
    if (!is.null(self_spec_list_fc)) {
      self_spec_list <- self_spec_list_fc
    } else {
      self_spec_list <- fit$function_call$self_spec_list
    }
  } else {
    self_spec_list <- NULL
    if (!is.null(self_spec_list_fc)) {
      stop(
        paste(
          "`self_spec_list_fc` was specified for the forecasting period,",
          "while a self specified component was not incorporated in the model."
        )
      )
    }
  }

  # Initialising list to return
  result <- list()

  # Number of observations
  N <- dim(fit$function_call$y)[1]

  # Number of dependent variables
  p <- dim(fit$function_call$y)[2]

  # Number of state parameters
  m <- length(fit$predicted$a_fc)

  # Initialising matrices and arrays for storing forecasted values
  # and corresponding variances
  y_fc <- matrix(0, forecast_period, p) # N_fc x p
  a_fc <- matrix(0, forecast_period, m) # N_fc x m
  P_fc <- array(0, dim = c(m, m, forecast_period)) # m x m x N_fc
  Fmat_fc <- array(0, dim = c(p, p, forecast_period)) # p x p x N_fc

  # Obtaining forecast for one step ahead to initialise the forecasting sequence
  a_fc[1,] <- fit$predicted$a_fc
  P_fc[,,1] <- fit$predicted$P_fc

  # Parameters used
  if (!is.null(fit$optim$par)) {
    param <- fit$optim$par
  } else {
    param <- fit$function_call$param
  }

  # Construct the system matrices
  sys_mat <- GetSysMat(p = p,
                       param = param,
                       update_part = TRUE,
                       add_residuals = FALSE,
                       H_format = fit$function_call$H_format,
                       local_level_ind = fit$function_call$local_level_ind,
                       slope_ind = fit$function_call$slope_ind,
                       BSM_vec = fit$function_call$BSM_vec,
                       cycle_ind = fit$function_call$cycle_ind,
                       addvar_list = addvar_list_fc,
                       level_addvar_list = level_addvar_list_fc,
                       arima_list = fit$function_call$arima_list,
                       sarima_list = fit$function_call$sarima_list,
                       self_spec_list = self_spec_list,
                       exclude_level = fit$function_call$exclude_level,
                       exclude_slope = fit$function_call$exclude_slope,
                       exclude_BSM_list = fit$function_call$exclude_BSM_list,
                       exclude_cycle_list = fit$function_call$exclude_cycle_list,
                       exclude_arima_list = fit$function_call$exclude_arima_list,
                       exclude_sarima_list = fit$function_call$exclude_sarima_list,
                       damping_factor_ind = fit$function_call$damping_factor_ind,
                       format_level = fit$function_call$format_level,
                       format_slope = fit$function_call$format_slope,
                       format_BSM_list = fit$function_call$format_BSM_list,
                       format_cycle_list = fit$function_call$format_cycle_list,
                       format_addvar = fit$function_call$format_addvar,
                       format_level_addvar = fit$function_call$format_level_addvar
  )

  # Z system matrices augmented with zeroes
  Z_padded <- sys_mat$Z_padded

  # Complete system matrices
  Z_full <- sys_mat$Z_kal
  T_full <- sys_mat$T_kal
  R_full <- sys_mat$R_kal
  Q_full <- sys_mat$Q_kal
  H_full <- sys_mat$H$H

  # Forecasting for t = 1 to forecast_period
  for (i in 1:forecast_period) {

    # System matrices for timepoint
    if (is.matrix(Z_full)) {
      Z_fc <- Z_full
    } else {
      Z_fc <- matrix(Z_full[,,i], nrow = p)
    }
    if (is.matrix(T_full)) {
      T_fc <- T_full
    } else {
      T_fc <- as.matrix(T_full[,,i])
    }
    if (is.matrix(R_full)) {
      R_fc <- R_full
    } else {
      R_fc <- matrix(R_full[,,i], nrow = m)
    }
    if (is.matrix(Q_full)) {
      Q_fc <- Q_full
    } else {
      Q_fc <- as.matrix(Q_full[,,i])
    }
    if (is.matrix(H_full)) {
      H_fc <- H_full
    } else {
      H_fc <- as.matrix(H_full[,,i])
    }

    # Forecast of y and corresponding uncertainty
    y_fc[i,] <- Z_fc %*% matrix(a_fc[i,])
    Fmat_fc[,,i] <- Z_fc %*% as.matrix(P_fc[,,i]) %*% t(Z_fc) + H_fc

    # Forecast of next state and corresponding uncertainty
    if (i < forecast_period) {
      a_fc[i + 1,] <- T_fc %*% matrix(a_fc[i,])
      P_fc[,,i + 1] <- T_fc %*% as.matrix(P_fc[,,i]) %*% t(T_fc) +
        R_fc %*% Q_fc %*% t(R_fc)
    }
  }

  #### Adjusting dimensions of Z matrices of components ####
  ##-- and adding predicted components of the model ------##

  # Local Level
  if (fit$function_call$local_level_ind & !fit$function_call$slope_ind &
      is.null(level_addvar_list_fc)) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    result$level <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      result$level[i,] <- tempZ %*% matrix(a_fc[i,])
    }
    Z_padded$level <- tempZ
  }

  # Local Level + Slope
  if (fit$function_call$slope_ind & is.null(level_addvar_list_fc)) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    result$level <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      result$level[i,] <- tempZ %*% matrix(a_fc[i,])
    }
    Z_padded$level <- tempZ
  }

  # BSM
  if (length(fit$function_call$BSM_vec) > 0) {
    for (s in fit$function_call$BSM_vec) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0('BSM', s)]])] <- Z_padded[[paste0('BSM', s)]]
      result[[paste0('BSM', s)]] <- matrix(0, forecast_period, p)
      for (i in 1:forecast_period) {
        result[[paste0('BSM', s)]][i,] <- tempZ %*% matrix(a_fc[i,])
      }
      Z_padded[[paste0('BSM', s)]] <- tempZ
    }
  }

  # Explanatory variables
  if (!is.null(addvar_list_fc)) {
    tempZ <- array(0, dim = c(p, m, forecast_period))
    result$addvar <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      tempZ[,,i][1:length(Z_padded$addvar[,,i])] <- Z_padded$addvar[,,i]
      result$addvar[i,] <- matrix(tempZ[,,i], nrow = p) %*% matrix(a_fc[i,])
    }
    Z_padded$addvar[,,i] <- tempZ
  }

  # level_addvar
  if (!is.null(level_addvar_list_fc) & !fit$function_call$slope_ind) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    result$level <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      result$level[i,] <- tempZ %*% matrix(a_fc[i,])
    }
    Z_padded$level <- tempZ
  }

  # slope_addvar
  if (!is.null(level_addvar_list_fc) & fit$function_call$slope_ind) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    result$level <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      result$level[i,] <- tempZ %*% matrix(a_fc[i,])
    }
    Z_padded$level <- tempZ
  }

  # Cycle
  if (fit$function_call$cycle_ind) {
    for (j in seq_along(fit$function_call$format_cycle_list)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0('Cycle', j)]])] <- Z_padded[[paste0('Cycle', j)]]
      result[[paste0('Cycle', j)]] <- matrix(0, forecast_period, p)
      for (i in 1:forecast_period) {
        result[[paste0('Cycle', j)]][i,] <- tempZ %*% matrix(a_fc[i,])
      }
      Z_padded[[paste0('Cycle', j)]] <- tempZ
    }
  }

  # ARIMA
  if (!is.null(fit$function_call$arima_list)) {
    for (j in seq_along(fit$function_call$arima_list)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0('ARIMA', j)]])] <- Z_padded[[paste0('ARIMA', j)]]
      result[[paste0('ARIMA', j)]] <- matrix(0, forecast_period, p)
      for (i in 1:forecast_period) {
        result[[paste0('ARIMA', j)]][i,] <- tempZ %*% matrix(a_fc[i,])
      }
      Z_padded[[paste0('ARIMA', j)]] <- tempZ
    }
  }

  # SARIMA
  if (!is.null(fit$function_call$sarima_list)) {
    for (j in seq_along(fit$function_call$sarima_list)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0('SARIMA', j)]])] <- Z_padded[[paste0('SARIMA', j)]]
      result[[paste0('SARIMA', j)]] <- matrix(0, forecast_period, p)
      for (i in 1:forecast_period) {
        result[[paste0('SARIMA', j)]][i,] <- tempZ %*% matrix(a_fc[i,])
      }
      Z_padded[[paste0('SARIMA', j)]] <- tempZ
    }
  }

  # Self Specified
  if (!is.null(fit$function_call$self_spec_list)) {
    if (is.matrix(Z_padded$self_spec)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded$self_spec)] <- Z_padded$self_spec
    } else {
      tempZ <- array(0, dim = c(p, m, forecast_period))
    }
    result$self_spec <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      if (is.matrix(Z_padded$self_spec)) {
        result$self_spec[i,] <- tempZ %*% matrix(a_fc[i,])
      } else {
        tempZ[,,i][1:length(Z_padded$self_spec[,,i])] <- Z_padded$self_spec[,,i]
        result$self_spec[i,] <- matrix(tempZ[,,i], nrow = p) %*% matrix(a_fc[i,])
      }
    }
    Z_padded$self_spec <- tempZ
  }

  result <- c(
    list(y_fc = y_fc, a_fc = a_fc, Fmat_fc = Fmat_fc, P_fc = P_fc),
    result, list(Z_padded = Z_padded)
  )
  return(result)
}
