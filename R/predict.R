#' State Space Model Forecasting
#'
#' Produces forecasts and out of sample simulations using a fitted State Space Model.
#'
#' @param object A statespacer object as returned by \code{\link{statespacer}}.
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
#'   from `self_spec_list` as passed on to \code{\link{statespacer}}. If some
#'   system matrices are time-varying then you should specify this argument.
#'   See \code{\link{statespacer}} for details about the format that must be
#'   followed for this argument.
#' @param forecast_period Number of time steps to forecast ahead.
#' @param nsim Number of simulations to generate over the forecast period.
#' @param ... Arguments passed on to the `predict` generic. Should not be used!
#'
#' @return
#' A list containing the forecasts and corresponding uncertainties.
#' In addition, it returns the components of the forecasts, as specified
#' by the State Space model.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#'
#' @references
#' \insertRef{durbin2012time}{statespacer}
#'
#' @examples
#' # Fit a SARIMA model on the AirPassengers data
#' library(datasets)
#' Data <- matrix(log(AirPassengers))
#' sarima_list <- list(list(s = c(12, 1), ar = c(0, 0), i = c(1, 1), ma = c(1, 1)))
#' fit <- statespacer(y = Data, 
#'                    H_format = matrix(0), 
#'                    sarima_list = sarima_list, 
#'                    initial = c(0.5*log(var(diff(Data))), 0, 0))
#'
#' # Obtain forecasts for 100 steps ahead using the fitted model
#' fc <- predict(fit, forecast_period = 100, nsim = 10)
#'
#' # Plot the forecasts and one of the simulation paths
#' plot(fc$y_fc, type = 'l')
#' lines(fc$sim$y[, 1, 1], type = 'p')
#' @export
predict.statespacer <- function(object,
                                addvar_list_fc = NULL,
                                level_addvar_list_fc = NULL,
                                self_spec_list_fc = NULL,
                                forecast_period = 1,
                                nsim = 0,
                                ...) {

  # Check if specification of addvar_list_fc is in line with the object
  if (!is.null(object$function_call$addvar_list) && is.null(addvar_list_fc)) {
    stop(
      "`addvar_list_fc` must be specified for the forecasting period.",
      call. = FALSE
    )
  }
  if (is.null(object$function_call$addvar_list) && !is.null(addvar_list_fc)) {
    stop(
      paste(
        "`addvar_list_fc` was specified for the forecasting period,",
        "while explanatory variables were not incorporated in the model."
      ),
      call. = FALSE
    )
  }

  # Check if specification of level_addvar_list_fc is in line with the object
  if (!is.null(object$function_call$level_addvar_list) &&
    is.null(level_addvar_list_fc)) {
    stop(
      "`level_addvar_list_fc` must be specified for the forecasting period.",
      call. = FALSE
    )
  }
  if (is.null(object$function_call$level_addvar_list) &&
    !is.null(level_addvar_list_fc)) {
    stop(
      paste(
        "`level_addvar_list_fc` was specified for the forecasting period,",
        "while explanatory variables in the level were not",
        "incorporated in the model."
      ),
      call. = FALSE
    )
  }

  # Check whether self_spec_list_fc should be used
  if (!is.null(object$function_call$self_spec_list)) {
    if (!is.null(self_spec_list_fc)) {
      self_spec_list <- self_spec_list_fc
    } else {
      self_spec_list <- object$function_call$self_spec_list
    }
  } else {
    self_spec_list <- NULL
    if (!is.null(self_spec_list_fc)) {
      stop(
        paste(
          "`self_spec_list_fc` was specified for the forecasting period,",
          "while a self specified component was not incorporated in the model."
        ),
        call. = FALSE
      )
    }
  }

  # Initialising list to return
  result <- list()
  
  # Number of observations
  N <- dim(object$function_call$y)[[1]]

  # Number of dependent variables
  p <- dim(object$function_call$y)[[2]]

  # Number of state parameters
  m <- length(object$predicted$a_fc)

  # Initialising matrices and arrays for storing forecasted values
  # and corresponding variances
  y_fc <- matrix(0, forecast_period, p) # N_fc x p
  a_fc <- matrix(0, forecast_period, m) # N_fc x m
  P_fc <- array(0, dim = c(m, m, forecast_period)) # m x m x N_fc
  Fmat_fc <- array(0, dim = c(p, p, forecast_period)) # p x p x N_fc
  
  # Initialising objects for storing simulations
  if (nsim > 0) {
    r <- dim(object$system_matrices$R$full)[[2]]
    sim_list <- list()
    y_sim <- array(0, dim = c(forecast_period, p, nsim))
    eta_sim <- array(0, dim = c(forecast_period, r, nsim))
    a_sim <- array(0, dim = c(forecast_period, m, nsim))
  }
  
  # Current last filled column indices of a and eta
  eta_index <- 0
  a_index <- 0

  # Obtaining forecast for one step ahead to initialise the forecasting sequence
  a_fc[1, ] <- object$predicted$a_fc
  P_fc[, , 1] <- object$predicted$P_fc

  # Parameters used
  if (!is.null(object$optim$par)) {
    param <- object$optim$par
  } else {
    param <- object$function_call$param
  }

  # Construct the system matrices
  sys_mat <- GetSysMat(
    p = p,
    param = param,
    update_part = TRUE,
    add_residuals = FALSE,
    H_format = object$function_call$H_format,
    local_level_ind = object$function_call$local_level_ind,
    slope_ind = object$function_call$slope_ind,
    BSM_vec = object$function_call$BSM_vec,
    cycle_ind = object$function_call$cycle_ind,
    addvar_list = addvar_list_fc,
    level_addvar_list = level_addvar_list_fc,
    arima_list = object$function_call$arima_list,
    sarima_list = object$function_call$sarima_list,
    self_spec_list = self_spec_list,
    exclude_level = object$function_call$exclude_level,
    exclude_slope = object$function_call$exclude_slope,
    exclude_BSM_list = object$function_call$exclude_BSM_list,
    exclude_cycle_list = object$function_call$exclude_cycle_list,
    exclude_arima_list = object$function_call$exclude_arima_list,
    exclude_sarima_list = object$function_call$exclude_sarima_list,
    damping_factor_ind = object$function_call$damping_factor_ind,
    format_level = object$function_call$format_level,
    format_slope = object$function_call$format_slope,
    format_BSM_list = object$function_call$format_BSM_list,
    format_cycle_list = object$function_call$format_cycle_list,
    format_addvar = object$function_call$format_addvar,
    format_level_addvar = object$function_call$format_level_addvar
  )

  # Z system matrices augmented with zeroes
  Z_padded <- sys_mat$Z_padded

  # Complete system matrices
  Z_full <- sys_mat$Z_kal
  T_full <- sys_mat$T_kal
  R_full <- sys_mat$R_kal
  Q_full <- sys_mat$Q_kal
  H_full <- sys_mat[["H"]][["H"]]

  # Forecasting for t = 1 to forecast_period
  for (i in 1:forecast_period) {

    # System matrices for time point
    if (is.matrix(Z_full)) {
      Z_fc <- Z_full
    } else {
      Z_fc <- matrix(Z_full[, , i], nrow = p)
    }
    if (is.matrix(T_full)) {
      T_fc <- T_full
    } else {
      T_fc <- as.matrix(T_full[, , i])
    }
    if (is.matrix(R_full)) {
      R_fc <- R_full
    } else {
      R_fc <- matrix(R_full[, , i], nrow = m)
    }
    if (is.matrix(Q_full)) {
      Q_fc <- Q_full
    } else {
      Q_fc <- as.matrix(Q_full[, , i])
    }
    if (is.matrix(H_full)) {
      H_fc <- H_full
    } else {
      H_fc <- as.matrix(H_full[, , i])
    }

    # Forecast of y and corresponding uncertainty
    y_fc[i, ] <- Z_fc %*% matrix(a_fc[i, ])
    Fmat_fc[, , i] <- tcrossprod(Z_fc %*% as.matrix(P_fc[, , i]), Z_fc) + H_fc

    # Forecast of next state and corresponding uncertainty
    if (i < forecast_period) {
      a_fc[i + 1, ] <- T_fc %*% matrix(a_fc[i, ])
      P_fc[, , i + 1] <- tcrossprod(T_fc %*% as.matrix(P_fc[, , i]), T_fc) +
        tcrossprod(R_fc %*% Q_fc, R_fc)
    }
  }

  #### Adjusting dimensions of Z matrices of components ####
  ## -- and adding predicted components of the model ------##

  # Local Level
  if (object$function_call$local_level_ind && !object$function_call$slope_ind &&
    is.null(level_addvar_list_fc)) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    result$level <- a_fc %*% tempZ_t
    Z_padded$level <- tempZ
    
    # Generate simulations
    if (nsim > 0) {
      
      # Model components
      Z <- object$system_matrices$Z$level
      Tmat <- object$system_matrices$T$level
      R <- object$system_matrices$R$level
      Q <- object$system_matrices$Q$level
      dim(Z) <- c(dim(Z), 1)
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)
      
      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num
      
      # Forecasted state and uncertainty
      a1 <- object$predicted$a_fc[a_indices]
      P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
      
      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = forecast_period,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P1,
        draw_initial = TRUE,
        eta_only = FALSE,
        transposed_state = FALSE
      )
      sim_list$level <- sim$y
      y_sim <- y_sim + sim$y
      eta_sim[ , eta_indices, ] <- sim$eta
      a_sim[ , a_indices, ] <- sim$a
    }
  }

  # Local Level + Slope
  if (object$function_call$slope_ind && is.null(level_addvar_list_fc)) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    result$level <- a_fc %*% tempZ_t
    Z_padded$level <- tempZ
    
    # Generate simulations
    if (nsim > 0) {
      
      # Model components
      Z <- object$system_matrices$Z$level
      Tmat <- object$system_matrices$T$level
      R <- object$system_matrices$R$level
      Q <- BlockMatrix(
        object$system_matrices$Q$level,
        object$system_matrices$Q$slope
      )
      dim(Z) <- c(dim(Z), 1)
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)
      
      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num
      
      # Forecasted state and uncertainty
      a1 <- object$predicted$a_fc[a_indices]
      P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
      
      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = forecast_period,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P1,
        draw_initial = TRUE,
        eta_only = FALSE,
        transposed_state = FALSE
      )
      sim_list$level <- sim$y
      y_sim <- y_sim + sim$y
      eta_sim[ , eta_indices, ] <- sim$eta
      a_sim[ , a_indices, ] <- sim$a
    }
  }

  # BSM
  if (length(object$function_call$BSM_vec) > 0) {
    for (s in object$function_call$BSM_vec) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0("BSM", s)]])] <- Z_padded[[paste0("BSM", s)]]
      tempZ_t <- t(tempZ)
      result[[paste0("BSM", s)]] <- a_fc %*% tempZ_t
      Z_padded[[paste0("BSM", s)]] <- tempZ
      
      # Generate simulations
      if (nsim > 0) {
        
        # Model components
        Z <- object$system_matrices$Z[[paste0("BSM", s)]]
        Tmat <- object$system_matrices$T[[paste0("BSM", s)]]
        R <- object$system_matrices$R[[paste0("BSM", s)]]
        Q <- object$system_matrices$Q[[paste0("BSM", s)]]
        dim(Z) <- c(dim(Z), 1)
        dim(Tmat) <- c(dim(Tmat), 1)
        dim(R) <- c(dim(R), 1)
        dim(Q) <- c(dim(Q), 1)
        
        # Indices of eta and a components
        eta_num <- dim(R)[[2]]
        a_num <- dim(Z)[[2]]
        eta_indices <- eta_index + 1:eta_num
        a_indices <- a_index + 1:a_num
        eta_index <- eta_index + eta_num
        a_index <- a_index + a_num
        
        # Forecasted state and uncertainty
        a1 <- object$predicted$a_fc[a_indices]
        P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
        
        # Obtain random samples of component
        sim <- SimulateC(
          nsim = nsim,
          repeat_Q = s - 1,
          N = forecast_period,
          a = a1,
          Z = Z,
          T = Tmat,
          R = R,
          Q = Q,
          P_star = P1,
          draw_initial = TRUE,
          eta_only = FALSE,
          transposed_state = FALSE
        )
        sim_list[[paste0("BSM", s)]] <- sim$y
        y_sim <- y_sim + sim$y
        eta_sim[ , eta_indices, ] <- sim$eta
        a_sim[ , a_indices, ] <- sim$a
      }
    }
  }

  # Explanatory variables
  if (!is.null(addvar_list_fc)) {
    tempZ <- array(0, dim = c(p, m, forecast_period))
    result$addvar <- matrix(0, forecast_period, p)
    for (i in 1:forecast_period) {
      tempZ[, , i][1:length(Z_padded$addvar[, , i])] <- Z_padded$addvar[, , i]
      result$addvar[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_fc[i, ])
    }
    Z_padded$addvar <- tempZ
    
    # Generate simulations
    if (nsim > 0) {
      
      # Model components
      Z <- sys_mat$Z$addvar
      Tmat <- object$system_matrices$T$addvar
      R <- object$system_matrices$R$addvar
      Q <- object$system_matrices$Q$addvar
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)
      
      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num
      
      # Forecasted state and uncertainty
      a1 <- object$predicted$a_fc[a_indices]
      P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
      
      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = forecast_period,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P1,
        draw_initial = TRUE,
        eta_only = FALSE,
        transposed_state = FALSE
      )
      sim_list$addvar <- sim$y
      y_sim <- y_sim + sim$y
      eta_sim[ , eta_indices, ] <- sim$eta
      a_sim[ , a_indices, ] <- sim$a
    }
  }

  # level_addvar
  if (!is.null(level_addvar_list_fc) && !object$function_call$slope_ind) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    result$level <- a_fc %*% tempZ_t
    Z_padded$level <- tempZ
    
    # Generate simulations
    if (nsim > 0) {
      
      # Model components
      Z <- object$system_matrices$Z$level
      Tmat <- sys_mat$Tmat$level
      R <- object$system_matrices$R$level
      Q <- BlockMatrix(
        object$system_matrices$Q$level,
        object$system_matrices$Q$level_addvar
      )
      dim(Z) <- c(dim(Z), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)
      
      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num
      
      # Forecasted state and uncertainty
      a1 <- object$predicted$a_fc[a_indices]
      P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
      
      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = forecast_period,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P1,
        draw_initial = TRUE,
        eta_only = FALSE,
        transposed_state = FALSE
      )
      sim_list$level <- sim$y
      y_sim <- y_sim + sim$y
      eta_sim[ , eta_indices, ] <- sim$eta
      a_sim[ , a_indices, ] <- sim$a
    }
  }

  # slope_addvar
  if (!is.null(level_addvar_list_fc) && object$function_call$slope_ind) {
    tempZ <- matrix(0, p, m)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    result$level <- a_fc %*% tempZ_t
    Z_padded$level <- tempZ
    
    # Generate simulations
    if (nsim > 0) {
      
      # Model components
      Z <- object$system_matrices$Z$level
      Tmat <- sys_mat$Tmat$level
      R <- object$system_matrices$R$level
      Q <- BlockMatrix(
        object$system_matrices$Q$level,
        object$system_matrices$Q$slope,
        object$system_matrices$Q$level_addvar
      )
      dim(Z) <- c(dim(Z), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)
      
      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num
      
      # Forecasted state and uncertainty
      a1 <- object$predicted$a_fc[a_indices]
      P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
      
      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = forecast_period,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P1,
        draw_initial = TRUE,
        eta_only = FALSE,
        transposed_state = FALSE
      )
      sim_list$level <- sim$y
      y_sim <- y_sim + sim$y
      eta_sim[ , eta_indices, ] <- sim$eta
      a_sim[ , a_indices, ] <- sim$a
    }
  }

  # Cycle
  if (object$function_call$cycle_ind) {
    for (j in seq_along(object$function_call$format_cycle_list)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0("Cycle", j)]])] <- Z_padded[[paste0("Cycle", j)]]
      tempZ_t <- t(tempZ)
      result[[paste0("Cycle", j)]] <- a_fc %*% tempZ_t
      Z_padded[[paste0("Cycle", j)]] <- tempZ
      
      # Generate simulations
      if (nsim > 0) {
        
        # Model components
        Z <- object$system_matrices$Z[[paste0("Cycle", j)]]
        Tmat <- object$system_matrices$T[[paste0("Cycle", j)]]
        R <- object$system_matrices$R[[paste0("Cycle", j)]]
        Q <- object$system_matrices$Q[[paste0("Cycle", j)]]
        dim(Z) <- c(dim(Z), 1)
        dim(Tmat) <- c(dim(Tmat), 1)
        dim(R) <- c(dim(R), 1)
        dim(Q) <- c(dim(Q), 1)

        # Indices of eta and a components
        eta_num <- dim(R)[[2]]
        a_num <- dim(Z)[[2]]
        eta_indices <- eta_index + 1:eta_num
        a_indices <- a_index + 1:a_num
        eta_index <- eta_index + eta_num
        a_index <- a_index + a_num
        
        # Forecasted state and uncertainty
        a1 <- object$predicted$a_fc[a_indices]
        P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]

        # Obtain random samples of component
        sim <- SimulateC(
          nsim = nsim,
          repeat_Q = 2,
          N = forecast_period,
          a = a1,
          Z = Z,
          T = Tmat,
          R = R,
          Q = Q,
          P_star = P1,
          draw_initial = TRUE,
          eta_only = FALSE,
          transposed_state = FALSE
        )
        sim_list[[paste0("Cycle", j)]] <- sim$y
        y_sim <- y_sim + sim$y
        eta_sim[ , eta_indices, ] <- sim$eta
        a_sim[ , a_indices, ] <- sim$a
      }
    }
  }

  # ARIMA
  if (!is.null(object$function_call$arima_list)) {
    for (j in seq_along(object$function_call$arima_list)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0("ARIMA", j)]])] <- Z_padded[[paste0("ARIMA", j)]]
      tempZ_t <- t(tempZ)
      result[[paste0("ARIMA", j)]] <- a_fc %*% tempZ_t
      Z_padded[[paste0("ARIMA", j)]] <- tempZ
      
      # Generate simulations
      if (nsim > 0) {
        
        # Model components
        Z <- object$system_matrices$Z[[paste0("ARIMA", j)]]
        Tmat <- object$system_matrices$T[[paste0("ARIMA", j)]]
        R <- object$system_matrices$R[[paste0("ARIMA", j)]]
        Q <- object$system_matrices$Q[[paste0("ARIMA", j)]]
        dim(Z) <- c(dim(Z), 1)
        dim(Tmat) <- c(dim(Tmat), 1)
        dim(R) <- c(dim(R), 1)
        dim(Q) <- c(dim(Q), 1)
        
        # Indices of eta and a components
        eta_num <- dim(R)[[2]]
        a_num <- dim(Z)[[2]]
        eta_indices <- eta_index + 1:eta_num
        a_indices <- a_index + 1:a_num
        eta_index <- eta_index + eta_num
        a_index <- a_index + a_num
        
        # Forecasted state and uncertainty
        a1 <- object$predicted$a_fc[a_indices]
        P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
        
        # Obtain random samples of component
        sim <- SimulateC(
          nsim = nsim,
          repeat_Q = 1,
          N = forecast_period,
          a = a1,
          Z = Z,
          T = Tmat,
          R = R,
          Q = Q,
          P_star = P1,
          draw_initial = TRUE,
          eta_only = FALSE,
          transposed_state = FALSE
        )
        sim_list[[paste0("ARIMA", j)]] <- sim$y
        y_sim <- y_sim + sim$y
        eta_sim[ , eta_indices, ] <- sim$eta
        a_sim[ , a_indices, ] <- sim$a
      }
    }
  }

  # SARIMA
  if (!is.null(object$function_call$sarima_list)) {
    for (j in seq_along(object$function_call$sarima_list)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded[[paste0("SARIMA", j)]])] <- Z_padded[[paste0("SARIMA", j)]]
      tempZ_t <- t(tempZ)
      result[[paste0("SARIMA", j)]] <- a_fc %*% tempZ_t
      Z_padded[[paste0("SARIMA", j)]] <- tempZ
      
      # Generate simulations
      if (nsim > 0) {
        
        # Model components
        Z <- object$system_matrices$Z[[paste0("SARIMA", j)]]
        Tmat <- object$system_matrices$T[[paste0("SARIMA", j)]]
        R <- object$system_matrices$R[[paste0("SARIMA", j)]]
        Q <- object$system_matrices$Q[[paste0("SARIMA", j)]]
        dim(Z) <- c(dim(Z), 1)
        dim(Tmat) <- c(dim(Tmat), 1)
        dim(R) <- c(dim(R), 1)
        dim(Q) <- c(dim(Q), 1)
        
        # Indices of eta and a components
        eta_num <- dim(R)[[2]]
        a_num <- dim(Z)[[2]]
        eta_indices <- eta_index + 1:eta_num
        a_indices <- a_index + 1:a_num
        eta_index <- eta_index + eta_num
        a_index <- a_index + a_num
        
        # Forecasted state and uncertainty
        a1 <- object$predicted$a_fc[a_indices]
        P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
        
        # Obtain random samples of component
        sim <- SimulateC(
          nsim = nsim,
          repeat_Q = 1,
          N = forecast_period,
          a = a1,
          Z = Z,
          T = Tmat,
          R = R,
          Q = Q,
          P_star = P1,
          draw_initial = TRUE,
          eta_only = FALSE,
          transposed_state = FALSE
        )
        sim_list[[paste0("SARIMA", j)]] <- sim$y
        y_sim <- y_sim + sim$y
        eta_sim[ , eta_indices, ] <- sim$eta
        a_sim[ , a_indices, ] <- sim$a
      }
    }
  }

  # Self Specified
  if (!is.null(object$function_call$self_spec_list)) {
    if (is.matrix(Z_padded$self_spec)) {
      tempZ <- matrix(0, p, m)
      tempZ[1:length(Z_padded$self_spec)] <- Z_padded$self_spec
      tempZ_t <- t(tempZ)
      result$self_spec <- a_fc %*% tempZ_t
    } else {
      tempZ <- array(0, dim = c(p, m, forecast_period))
      result$self_spec <- matrix(0, forecast_period, p)
      for (i in 1:forecast_period) {
        tempZ[, , i][1:length(Z_padded$self_spec[, , i])] <- Z_padded$self_spec[, , i]
        result$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_fc[i, ])
      }
    }
    Z_padded$self_spec <- tempZ
    
    # Generate simulations
    if (nsim > 0) {
      
      # Model components
      Z <- sys_mat$Z$self_spec
      Tmat <- sys_mat$Tmat$self_spec
      R <- sys_mat$R$self_spec
      Q <- sys_mat$Q$self_spec
      if (is.matrix(Z)) {
        dim(Z) <- c(dim(Z), 1)
      }
      if (is.matrix(Tmat)) {
        dim(Tmat) <- c(dim(Tmat), 1)
      }
      if (is.matrix(R)) {
        dim(R) <- c(dim(R), 1)
      }
      if (is.matrix(Q)) {
        dim(Q) <- c(dim(Q), 1)
      }
      
      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num
      
      # Forecasted state and uncertainty
      a1 <- object$predicted$a_fc[a_indices]
      P1 <- object$predicted$P_fc[a_indices, a_indices, drop = FALSE]
      
      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = forecast_period,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P1,
        draw_initial = TRUE,
        eta_only = FALSE,
        transposed_state = FALSE
      )
      sim_list$self_spec <- sim$y
      y_sim <- y_sim + sim$y
      eta_sim[ , eta_indices, ] <- sim$eta
      a_sim[ , a_indices, ] <- sim$a
    }
  }

  if (nsim > 0) {
    ### Generate Simulated Residuals ###
    
    # Model components
    Q <- sys_mat$H$H
    if (is.matrix(Q)) {
      dim(Q) <- c(dim(Q), 1)
    }
    
    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = forecast_period,
      a = matrix(0),
      Z = array(0, dim = c(1, 1, 1)),
      T = array(0, dim = c(1, 1, 1)),
      R = array(0, dim = c(1, 1, 1)),
      Q = Q,
      P_star = matrix(0),
      draw_initial = FALSE,
      eta_only = TRUE,
      transposed_state = FALSE
    )
    sim_list$epsilon <- sim$eta
    y_sim <- y_sim + sim$eta
  }
  
  result <- c(
    list(y_fc = y_fc, a_fc = a_fc, Fmat_fc = Fmat_fc, P_fc = P_fc),
    result, list(Z_padded = Z_padded)
  )
  
  # Store simulated objects
  if (nsim > 0) {
    sim_list$y <- y_sim
    sim_list$a <- a_sim
    sim_list$eta <- eta_sim
    result$sim <- sim_list
  }

  return(result)
}
