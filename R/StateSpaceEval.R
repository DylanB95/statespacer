#' State Space Model Evaluation at Specified Parameters
#'
#' Evaluates the specified State Space model at the parameters
#' specified by the user.
#'
#' @param param Parameters used to construct the system matrices.
#' @inheritParams statespacer
#'
#' @return
#' A list containing:
#' * `function_call`: A list containing the input to the function.
#' * `system_matrices`: A list containing the system matrices of
#'   the State Space model.
#' * `predicted`: A list containing the predicted components of
#'   the State Space model.
#' * `filtered`: A list containing the filtered components of
#'   the State Space model.
#' * `smoothed`: A list containing the smoothed components of
#'   the State Space model.
#' * `diagnostics`: A list containing items useful for diagnostical tests.
#'
#' @noRd
StateSpaceEval <- function(param,
                           y,
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
                           exclude_BSM_list = lapply(BSM_vec, function(x) 0),
                           exclude_cycle_list = list(0),
                           exclude_arima_list = lapply(arima_list, function(x) 0),
                           exclude_sarima_list = lapply(sarima_list, function(x) 0),
                           damping_factor_ind = rep(TRUE, length(exclude_cycle_list)),
                           format_level = NULL,
                           format_slope = NULL,
                           format_BSM_list = lapply(BSM_vec, function(x) NULL),
                           format_cycle_list = lapply(exclude_cycle_list, function(x) NULL),
                           format_addvar = NULL,
                           format_level_addvar = NULL,
                           diagnostics = TRUE) {

  # Initialising lists to return
  filtered <- list()
  predicted <- list()
  smoothed <- list()
  system_matrices <- list()
  diagnostics_ind <- diagnostics # Storing diagnostics indicator
  diagnostics <- list()

  # N = Number of observations
  N <- dim(y)[[1]]

  # p = Number of dependent variables
  p <- dim(y)[[2]]

  # Construct the system matrices
  sys_mat <- GetSysMat(
    p = p,
    param = param,
    update_part = TRUE,
    add_residuals = TRUE,
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

  # Adjust input arguments
  H_format <- sys_mat$function_call$H_format
  local_level_ind <- sys_mat$function_call$local_level_ind
  slope_ind <- sys_mat$function_call$slope_ind
  BSM_vec <- sys_mat$function_call$BSM_vec
  cycle_ind <- sys_mat$function_call$cycle_ind
  addvar_list <- sys_mat$function_call$addvar_list
  level_addvar_list <- sys_mat$function_call$level_addvar_list
  arima_list <- sys_mat$function_call$arima_list
  sarima_list <- sys_mat$function_call$sarima_list
  self_spec_list <- sys_mat$function_call$self_spec_list
  exclude_level <- sys_mat$function_call$exclude_level
  exclude_slope <- sys_mat$function_call$exclude_slope
  exclude_BSM_list <- sys_mat$function_call$exclude_BSM_list
  exclude_cycle_list <- sys_mat$function_call$exclude_cycle_list
  exclude_arima_list <- sys_mat$function_call$exclude_arima_list
  exclude_sarima_list <- sys_mat$function_call$exclude_sarima_list
  damping_factor_ind <- sys_mat$function_call$damping_factor_ind
  format_level <- sys_mat$function_call$format_level
  format_slope <- sys_mat$function_call$format_slope
  format_BSM_list <- sys_mat$function_call$format_BSM_list
  format_cycle_list <- sys_mat$function_call$format_cycle_list
  format_addvar <- sys_mat$function_call$format_addvar
  format_level_addvar <- sys_mat$function_call$format_level_addvar

  # Z system matrices augmented with zeroes
  Z_padded <- sys_mat$Z_padded

  # Complete system matrices
  Z_kal <- sys_mat$Z_kal
  T_kal <- sys_mat$T_kal
  R_kal <- sys_mat$R_kal
  Q_kal <- sys_mat$Q_kal
  a <- sys_mat$a_kal
  P_inf <- sys_mat$P_inf_kal
  P_star <- sys_mat$P_star_kal

  # NA in y
  y_isna <- is.na(y)

  # Augment system matrices to arrays
  if (is.matrix(Z_kal)) {
    dim(Z_kal) <- c(dim(Z_kal), 1)
  }
  if (is.matrix(T_kal)) {
    dim(T_kal) <- c(dim(T_kal), 1)
  }
  if (is.matrix(R_kal)) {
    dim(R_kal) <- c(dim(R_kal), 1)
  }
  if (is.matrix(Q_kal)) {
    dim(Q_kal) <- c(dim(Q_kal), 1)
  }

  # Number of state parameters
  m <- dim(a)[[1]]

  # Number of diffuse elements
  diffuse_num <- sum(P_inf)

  # Kalman Filter and Smoother
  kalman <- KalmanC(
    y = y,
    y_isna = y_isna,
    a = a,
    P_inf = P_inf,
    P_star = P_star,
    Z = Z_kal,
    T = T_kal,
    R = R_kal,
    Q = Q_kal,
    diagnostics = diagnostics_ind
  )

  # Assign objects computed by KalmanC
  initialisation_steps <- kalman$nested$initialisation_steps
  loglik <- kalman$nested$loglik
  a_pred <- kalman$nested$a_pred
  a_fil <- kalman$nested$a_fil
  a_smooth <- kalman$nested$a_smooth
  P_pred <- kalman$nested$P_pred
  P_fil <- kalman$nested$P_fil
  V <- kalman$nested$V
  P_inf_pred <- kalman$nested$P_inf_pred
  P_star_pred <- kalman$P_star_pred
  P_inf_fil <- kalman$P_inf_fil
  P_star_fil <- kalman$P_star_fil
  yfit <- kalman$yfit
  v <- kalman$v
  eta <- kalman$eta
  Fmat <- kalman$Fmat
  eta_var <- kalman$eta_var
  a_fc <- kalman$a_fc
  P_inf_fc <- kalman$P_inf_fc
  P_star_fc <- kalman$P_star_fc
  P_fc <- kalman$P_fc
  r_vec <- kalman$r_vec
  Nmat <- kalman$Nmat
  if (diagnostics_ind) {
    v_norm <- kalman$v_norm
    e <- kalman$e
    D <- kalman$D
    Tstat_observation <- kalman$Tstat_observation
    Tstat_state <- kalman$Tstat_state
  }

  # Diagnostics
  if (diagnostics_ind) {

    # Calculating test statistics based on the normalised prediction errors
    obs_vec <- colSums(!is.na(v_norm))
    m1 <- colMeans(v_norm, na.rm = TRUE)
    v_norm_centered <- t(t(v_norm) - m1)
    m2 <- colMeans(v_norm_centered^2, na.rm = TRUE)
    m3 <- colMeans(v_norm_centered^3, na.rm = TRUE)
    m4 <- colMeans(v_norm_centered^4, na.rm = TRUE)
    Sstat <- m3 / sqrt(m2^3)
    Kstat <- m4 / m2^2
    Nstat <- obs_vec * (Sstat^2 / 6 + (Kstat - 3)^2 / 24)

    # Number of observations after initialisation
    obs <- min(obs_vec)

    # Correlogram and Box-Ljung statistic
    if (obs > 0) {
      correlogram <- matrix(0, floor(obs / 2), p)
      Box_Ljung <- matrix(0, floor(obs / 2), p)
      BL_running <- 0
      for (i in 1:floor(obs / 2)) {
        correlogram[i, ] <- colSums(
          v_norm_centered[(i + 1):N, , drop = FALSE] *
            v_norm_centered[1:(N - i), , drop = FALSE],
          na.rm = TRUE
        ) / (obs_vec * m2)
        BL_running <- BL_running +
          obs_vec * (obs_vec + 2) * correlogram[i, ]^2 / (obs_vec - i)
        Box_Ljung[i, ] <- BL_running
      }

      # Heteroscedasticity test
      Heteroscedasticity <- matrix(0, floor(obs / 3), p)
      group1 <- 0
      group2 <- 0
      v_2 <- v_norm[stats::complete.cases(v_norm), , drop = FALSE]
      for (i in 1:floor(obs / 3)) {
        group1 <- group1 + v_2[i, ]^2
        group2 <- group2 + v_2[dim(v_2)[[1]] + 1 - i, ]^2
        Heteroscedasticity[i, ] <- group2 / group1
      }
    }
  }

  # Smoothed observation disturbance and corresponding variance
  epsilon <- a_smooth[, sys_mat$residuals_state, drop = FALSE]
  epsilon_var <- V[sys_mat$residuals_state, sys_mat$residuals_state, , drop = FALSE]

  # Removing residuals
  a_pred <- a_pred[, -sys_mat$residuals_state, drop = FALSE]
  a_fil <- a_fil[, -sys_mat$residuals_state, drop = FALSE]
  a_smooth <- a_smooth[, -sys_mat$residuals_state, drop = FALSE]
  P_pred <- P_pred[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_fil <- P_fil[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  V <- V[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_inf_pred <- P_inf_pred[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_star_pred <- P_star_pred[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_inf_fil <- P_inf_fil[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_star_fil <- P_star_fil[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  a_fc <- a_fc[-sys_mat$residuals_state, drop = FALSE]
  P_fc <- P_fc[-sys_mat$residuals_state, -sys_mat$residuals_state, drop = FALSE]
  P_inf_fc <- P_inf_fc[-sys_mat$residuals_state, -sys_mat$residuals_state, drop = FALSE]
  P_star_fc <- P_star_fc[-sys_mat$residuals_state, -sys_mat$residuals_state, drop = FALSE]
  r_vec <- r_vec[-1, -sys_mat$residuals_state, drop = FALSE]
  Nmat <- Nmat[-sys_mat$residuals_state, -sys_mat$residuals_state, -1, drop = FALSE]
  eta <- eta[, -sys_mat$residuals_state, drop = FALSE]
  eta_var <- eta_var[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  if (diagnostics_ind) {
    Tstat_state <- Tstat_state[-1, -sys_mat$residuals_state, drop = FALSE]
  }

  # Storing components
  predicted$yfit <- yfit
  predicted$v <- v
  predicted$Fmat <- Fmat
  predicted$a <- a_pred
  predicted$P <- P_pred
  predicted$P_inf <- P_inf_pred
  predicted$P_star <- P_star_pred
  predicted$a_fc <- a_fc
  predicted$P_fc <- P_fc
  predicted$P_inf_fc <- P_inf_fc
  predicted$P_star_fc <- P_star_fc
  filtered$a <- a_fil
  filtered$P <- P_fil
  filtered$P_inf <- P_inf_fil
  filtered$P_star <- P_star_fil
  smoothed$a <- a_smooth
  smoothed$V <- V
  smoothed$eta <- eta
  smoothed$eta_var <- eta_var
  smoothed$epsilon <- epsilon
  smoothed$epsilon_var <- epsilon_var
  diagnostics$initialisation_steps <- initialisation_steps
  diagnostics$loglik <- sum(loglik, na.rm = TRUE)
  diagnostics$AIC <- 1 / N *
    (-2 * diagnostics$loglik + 2 * (diffuse_num + length(param)))
  diagnostics$BIC <- 1 / N *
    (-2 * diagnostics$loglik + (diffuse_num + length(param)) * log(N))
  diagnostics$r <- r_vec
  diagnostics$N <- Nmat
  if (diagnostics_ind) {
    diagnostics$e <- e
    diagnostics$D <- D
    diagnostics$Tstat_observation <- Tstat_observation
    diagnostics$Tstat_state <- Tstat_state
    diagnostics$v_normalised <- v_norm
    diagnostics$Skewness <- Sstat
    diagnostics$Kurtosis <- Kstat
    diagnostics$Jarque_Bera <- Nstat
    diagnostics$Jarque_Bera_criticalvalue <- stats::qchisq(0.95, df = 2)
    if (obs > 0) {
      diagnostics$correlogram <- correlogram
      diagnostics$Box_Ljung <- Box_Ljung
      diagnostics$Box_Ljung_criticalvalues <- matrix(
        stats::qchisq(0.95, df = 1:floor(obs / 2))
      )
      diagnostics$Heteroscedasticity <- Heteroscedasticity
      diagnostics$Heteroscedasticity_criticalvalues <- cbind(
        stats::qf(0.025, df1 = 1:floor(obs / 3), df2 = 1:floor(obs / 3)),
        stats::qf(0.975, df1 = 1:floor(obs / 3), df2 = 1:floor(obs / 3))
      )
    }
  }

  #### Adjusting dimensions of Z matrices of components ####
  ## -- and adding fitted components of the model --------##

  # Local Level
  if (local_level_ind && !slope_ind && is.null(level_addvar_list)) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    predicted$level <- a_pred %*% tempZ_t
    filtered$level <- a_fil %*% tempZ_t
    smoothed$level <- a_smooth %*% tempZ_t
    Z_padded$level <- tempZ
  }

  # Local Level + Slope
  if (slope_ind && is.null(level_addvar_list)) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    predicted$level <- a_pred %*% tempZ_t
    filtered$level <- a_fil %*% tempZ_t
    smoothed$level <- a_smooth %*% tempZ_t
    Z_padded$level <- tempZ
  }

  # BSM
  if (length(BSM_vec) > 0) {
    for (s in BSM_vec) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("BSM", s)]])] <- Z_padded[[paste0("BSM", s)]]
      tempZ_t <- t(tempZ)
      predicted[[paste0("BSM", s)]] <- a_pred %*% tempZ_t
      filtered[[paste0("BSM", s)]] <- a_fil %*% tempZ_t
      smoothed[[paste0("BSM", s)]] <- a_smooth %*% tempZ_t
      Z_padded[[paste0("BSM", s)]] <- tempZ
    }
  }

  # Explanatory variables
  if (!is.null(addvar_list)) {
    tempZ <- array(0, dim = c(p, m - p, N))
    predicted$addvar <- matrix(0, N, p)
    filtered$addvar <- matrix(0, N, p)
    smoothed$addvar <- matrix(0, N, p)
    for (i in 1:N) {
      tempZ[, , i][1:length(Z_padded$addvar[, , i])] <- Z_padded$addvar[, , i]
      predicted$addvar[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_pred[i, ])
      filtered$addvar[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_fil[i, ])
      smoothed$addvar[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_smooth[i, ])
    }
    Z_padded$addvar <- tempZ
    filtered$addvar_coeff <- a_fil[, sys_mat$addvar_state]
    filtered$addvar_coeff_se <- t(apply(
      P_fil, 3, function(x) sqrt(diag(as.matrix(x))[sys_mat$addvar_state])
    ))
    smoothed$addvar_coeff <- a_smooth[, sys_mat$addvar_state]
    smoothed$addvar_coeff_se <- t(
      apply(V, 3, function(x) sqrt(diag(as.matrix(x))[sys_mat$addvar_state]))
    )
  }

  # level_addvar
  if (!is.null(level_addvar_list) && !slope_ind) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    predicted$level <- a_pred %*% tempZ_t
    filtered$level <- a_fil %*% tempZ_t
    smoothed$level <- a_smooth %*% tempZ_t
    Z_padded$level <- tempZ
    filtered$level_addvar_coeff <- a_fil[, sys_mat$level_addvar_state]
    filtered$level_addvar_coeff_se <- t(apply(
      P_fil, 3, function(x) sqrt(diag(as.matrix(x))[sys_mat$level_addvar_state])
    ))
    smoothed$level_addvar_coeff <- a_smooth[, sys_mat$level_addvar_state]
    smoothed$level_addvar_coeff_se <- t(apply(
      V, 3, function(x) sqrt(diag(as.matrix(x))[sys_mat$level_addvar_state])
    ))
  }

  # slope_addvar
  if (!is.null(level_addvar_list) && slope_ind) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    tempZ_t <- t(tempZ)
    predicted$level <- a_pred %*% tempZ_t
    filtered$level <- a_fil %*% tempZ_t
    smoothed$level <- a_smooth %*% tempZ_t
    Z_padded$level <- tempZ
    filtered$level_addvar_coeff <- a_fil[, sys_mat$level_addvar_state]
    filtered$level_addvar_coeff_se <- t(apply(
      P_fil, 3, function(x) sqrt(diag(as.matrix(x))[sys_mat$level_addvar_state])
    ))
    smoothed$level_addvar_coeff <- a_smooth[, sys_mat$level_addvar_state]
    smoothed$level_addvar_coeff_se <- t(apply(
      V, 3, function(x) sqrt(diag(as.matrix(x))[sys_mat$level_addvar_state])
    ))
  }

  # Cycle
  if (cycle_ind) {
    for (j in seq_along(format_cycle_list)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("Cycle", j)]])] <- Z_padded[[paste0("Cycle", j)]]
      tempZ_t <- t(tempZ)
      predicted[[paste0("Cycle", j)]] <- a_pred %*% tempZ_t
      filtered[[paste0("Cycle", j)]] <- a_fil %*% tempZ_t
      smoothed[[paste0("Cycle", j)]] <- a_smooth %*% tempZ_t
      Z_padded[[paste0("Cycle", j)]] <- tempZ
    }
  }

  # ARIMA
  if (!is.null(arima_list)) {
    for (j in seq_along(arima_list)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("ARIMA", j)]])] <- Z_padded[[paste0("ARIMA", j)]]
      tempZ_t <- t(tempZ)
      predicted[[paste0("ARIMA", j)]] <- a_pred %*% tempZ_t
      filtered[[paste0("ARIMA", j)]] <- a_fil %*% tempZ_t
      smoothed[[paste0("ARIMA", j)]] <- a_smooth %*% tempZ_t
      Z_padded[[paste0("ARIMA", j)]] <- tempZ
    }
  }

  # SARIMA
  if (!is.null(sarima_list)) {
    for (j in seq_along(sarima_list)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("SARIMA", j)]])] <- Z_padded[[paste0("SARIMA", j)]]
      tempZ_t <- t(tempZ)
      predicted[[paste0("SARIMA", j)]] <- a_pred %*% tempZ_t
      filtered[[paste0("SARIMA", j)]] <- a_fil %*% tempZ_t
      smoothed[[paste0("SARIMA", j)]] <- a_smooth %*% tempZ_t
      Z_padded[[paste0("SARIMA", j)]] <- tempZ
    }
  }

  # Self Specified
  if (!is.null(self_spec_list)) {
    if (is.matrix(Z_padded$self_spec)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded$self_spec)] <- Z_padded$self_spec
      tempZ_t <- t(tempZ)
      predicted$self_spec <- a_pred %*% tempZ_t
      filtered$self_spec <- a_fil %*% tempZ_t
      smoothed$self_spec <- a_smooth %*% tempZ_t
    } else {
      tempZ <- array(0, dim = c(p, m - p, N))
      predicted$self_spec <- matrix(0, N, p)
      filtered$self_spec <- matrix(0, N, p)
      smoothed$self_spec <- matrix(0, N, p)
      for (i in 1:N) {
        tempZ[, , i][1:length(Z_padded$self_spec[, , i])] <- Z_padded$self_spec[, , i]
        predicted$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_pred[i, ])
        filtered$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_fil[i, ])
        smoothed$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_smooth[i, ])
      }
    }
    Z_padded$self_spec <- tempZ
  }

  # Filling system_matrices
  system_matrices[["H"]] <- sys_mat[["H"]]
  system_matrices$Z <- sys_mat$Z
  system_matrices$T <- sys_mat$Tmat
  system_matrices$R <- sys_mat$R
  system_matrices[["Q"]] <- sys_mat[["Q"]]
  system_matrices$Q_loading_matrix <- sys_mat$Q_loading_matrix
  system_matrices$Q_diagonal_matrix <- sys_mat$Q_diagonal_matrix
  system_matrices$Q_correlation_matrix <- sys_mat$Q_correlation_matrix
  system_matrices$Q_stdev_matrix <- sys_mat$Q_stdev_matrix
  system_matrices$lambda <- sys_mat$lambda
  system_matrices$rho <- sys_mat$rho
  system_matrices$AR <- sys_mat$AR
  system_matrices$MA <- sys_mat$MA
  system_matrices$SAR <- sys_mat$SAR
  system_matrices$SMA <- sys_mat$SMA
  system_matrices$self_spec <- sys_mat$self_spec
  system_matrices$a1 <- sys_mat$a1
  system_matrices$P_inf <- sys_mat$P_inf
  system_matrices$P_star <- sys_mat$P_star
  system_matrices$Z_padded <- Z_padded
  system_matrices$state_label <- sys_mat$state_label

  # Returning the result
  result <- list()
  result$function_call <- c(
    list(y = y),
    sys_mat$function_call,
    list(fit = FALSE, initial = param, diagnostics = diagnostics_ind)
  )
  result$system_matrices <- system_matrices
  result$predicted <- predicted
  result$filtered <- filtered
  result$smoothed <- smoothed
  result$diagnostics <- diagnostics
  return(result)
}
