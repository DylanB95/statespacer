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
                           loglik_only = FALSE) {

  ##### Initialising lists to return #####
  filtered <- list()
  predicted <- list()
  smoothed <- list()
  system_matrices <- list()
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

  # Uncertainty of initial 'guess' of state vector
  kappa <- 1e7

  # Number of state parameters
  m <- dim(a)[[1]]

  # Number of disturbances in the state equation
  r <- dim(R_kal)[[2]]

  # Number of diffuse elements
  diffuse_num <- sum(P_inf)

  # Check if P_inf is already 0
  initialisation <- !all(abs(P_inf) < 1e-7)

  # Initial P
  P <- kappa * P_inf + P_star

  # Initialise (filtered) state and corresponding variance
  a <- t(matrix(a, m, N * p)) # N*p x m
  P_inf <- array(P_inf, dim = c(m, m, N * p)) # m x m x N*p
  P_inf[, , 2:(N * p)] <- 0
  P_star <- array(P_star, dim = c(m, m, N * p)) # m x m x N*p
  a_pred <- t(matrix(sys_mat$a_kal, m, N)) # N x m
  P_pred <- array(P, dim = c(m, m, N)) # m x m x N
  a_fil <- a_pred # N x m
  P_fil <- P_pred # m x m x N

  ################################# Kalman Filter #################################

  # Initialising residuals and variance
  v <- matrix(0, N, p) # N x p
  Fmat <- array(0, dim = c(p, p, N)) # p x p x N

  # Initialising normalised residuals
  v_norm <- matrix(NA, N, p) # N x p

  # Initialising loglikelihood vector
  loglik <- rep(0, N * p)

  # Initialising fitted values matrix
  yfit <- matrix(0, N, p) # N x p

  ###### Applying Kalman Filter with exact initialisation ######

  # Keep track of the number of initialisation steps
  initialisation_steps <- 0

  # Keep track of which time index and row index to use
  # For selection of system matrices
  t <- 1
  row <- 0

  # Applying Kalman Filter
  for (i in 1:(N * p)) {

    # Updating indices
    row <- row + 1
    if (row == (p + 1)) {
      row <- 1
      t <- t + 1
    }

    # When should a transition to the next timepoint be made?
    if (i %% p == 0) {
      timestep <- TRUE

      # T, R, and Q matrices only needed when a transition to
      # the next timepoint is made
      if (is.matrix(T_kal)) {
        T_input <- T_kal
      } else {
        T_input <- as.matrix(T_kal[, , t])
      }
      if (is.matrix(R_kal)) {
        R_input <- R_kal
      } else {
        R_input <- matrix(R_kal[, , t], nrow = m)
      }
      if (is.matrix(Q_kal)) {
        Q_input <- Q_kal
      } else {
        Q_input <- as.matrix(Q_kal[, , t])
      }
    } else {
      timestep <- FALSE
      T_input <- NULL
      R_input <- NULL
      Q_input <- NULL
    }

    if (is.matrix(Z_kal)) {
      Z_input <- Z_kal[row, , drop = FALSE]
    } else {
      Z_input <- matrix(Z_kal[row, , t], nrow = 1)
    }

    # Apply KalmanEI in initialisation stage, else KalmanUT
    if (initialisation) {

      # Keep track of the number of initialisation steps
      initialisation_steps <- initialisation_steps + 1

      # Calling the Kalman Filter with exact initialisation
      filter_output <- KalmanEI(
        y = y[[t, row]],
        a = matrix(a[i, ]),
        P_inf = as.matrix(P_inf[, , i]),
        P_star = as.matrix(P_star[, , i]),
        Z = Z_input,
        Tmat = T_input,
        R = R_input,
        Q = Q_input,
        timestep = timestep
      )

      # Storing next predicted state and variance used for the next iteration
      if (i < (N * p)) {
        a[i + 1, ] <- filter_output$a
        P_inf[, , i + 1] <- filter_output$P_inf
        P_star[, , i + 1] <- filter_output$P_star

        # Check if P_inf converged to zero
        initialisation <- !all(abs(filter_output$P_inf) < 1e-7)
      }

      # Saving predicted and filtered state and corresponding variance for each timestep
      if (timestep & !loglik_only) {

        # Predicted state and variance
        if (i < (N * p)) {
          a_pred[t + 1, ] <- filter_output$a
          if (initialisation) {
            P_pred[, , t + 1] <- kappa * filter_output$P_inf + filter_output$P_star
          } else {
            P_pred[, , t + 1] <- filter_output$P_star
          }
        } else {

          # Store out of sample forecast separately
          a_fc <- filter_output$a
          if (initialisation) {
            P_fc <- kappa * filter_output$P_inf + filter_output$P_star
          } else {
            P_fc <- filter_output$P_star
          }
        }

        # Filtered state and variance
        a_fil[t, ] <- filter_output$a_fil
        if (initialisation) {
          P_fil[, , t] <- kappa * filter_output$P_inf_fil + filter_output$P_star_fil
        } else {
          P_fil[, , t] <- filter_output$P_star_fil
        }

        if (is.matrix(Z_kal)) {
          Z_full <- Z_kal
        } else {
          Z_full <- matrix(Z_kal[, , t], nrow = p)
        }

        # Storing fitted values
        yfit[t, ] <- Z_full %*% matrix(a_pred[t, ])

        # Storing residuals
        v[t, ] <- y[t, ] - yfit[t, ]

        # Storing Variance - Covariance matrix of fitted values
        Fmat[, , t] <- tcrossprod(Z_full %*% as.matrix(P_pred[, , t]), Z_full)

        if (all(abs(tcrossprod(
          Z_full %*% as.matrix(P_inf[, , i - (p - 1)]),
          Z_full
        )) < 1e-7)) {

          # Inverse of Fmat
          Finv <- solve(Fmat[, , t])

          # Singular Value decomposition to obtain square root of the inverse of Fmat
          svd_Finv <- svd(Finv)

          # Square root of Inverse of Fmat
          Finv_root <- tcrossprod(
            svd_Finv$u %*% sqrt(diag(svd_Finv$d, p, p)),
            svd_Finv$u
          )

          # Normalised prediction error
          v_0na <- v[t, ]
          v_0na[is.na(v_0na)] <- 0
          v_norm[t, ] <- Finv_root %*% matrix(v_0na)
          v_norm[t, ][is.na(v[t, ])] <- NA
        }
      }
    } else {

      # Calling the Kalman Filter
      filter_output <- KalmanUT(
        y = y[[t, row]],
        a = matrix(a[i, ]),
        P = as.matrix(P_star[, , i]),
        Z = Z_input,
        Tmat = T_input,
        R = R_input,
        Q = Q_input,
        timestep = timestep
      )

      # Storing next predicted state and variance used for the next iteration
      if (i < (N * p)) {
        a[i + 1, ] <- filter_output$a
        P_star[, , i + 1] <- filter_output$P
      }

      # Saving predicted and filtered state and corresponding variance
      # for each timestep
      if (timestep & !loglik_only) {

        # Predicted state and variance
        if (i < (N * p)) {
          a_pred[t + 1, ] <- filter_output$a
          P_pred[, , t + 1] <- filter_output$P
        } else {

          # Store out of sample forecast separately
          a_fc <- filter_output$a
          P_fc <- filter_output$P
        }

        # Filtered state and variance
        a_fil[t, ] <- filter_output$a_fil
        P_fil[, , t] <- filter_output$P_fil

        if (is.matrix(Z_kal)) {
          Z_full <- Z_kal
        } else {
          Z_full <- matrix(Z_kal[, , t], nrow = p)
        }

        # Storing fitted values
        yfit[t, ] <- Z_full %*% matrix(a_pred[t, ])

        # Storing residuals
        v[t, ] <- y[t, ] - yfit[t, ]

        # Storing Variance - Covariance matrix of fitted values
        Fmat[, , t] <- tcrossprod(Z_full %*% as.matrix(P_pred[, , t]), Z_full)

        # Inverse of Fmat
        Finv <- solve(Fmat[, , t])

        # Singular Value decomposition to obtain square root of the inverse of Fmat
        svd_Finv <- svd(Finv)

        # Square root of Inverse of Fmat
        Finv_root <- tcrossprod(
          svd_Finv$u %*% sqrt(diag(svd_Finv$d, p, p)),
          svd_Finv$u
        )

        # Normalised prediction error
        v_0na <- v[t, ]
        v_0na[is.na(v_0na)] <- 0
        v_norm[t, ] <- Finv_root %*% matrix(v_0na)
        v_norm[t, ][is.na(v[t, ])] <- NA
      }
    }

    # Store loglikelihood
    loglik[[i]] <- filter_output$loglik
  }

  # Should only the loglikelihood be returned?
  if (loglik_only) {
    return(sum(loglik, na.rm = TRUE))
  }

  #### Diagnostics ####

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
    v_2 <- v_norm[complete.cases(v_norm), , drop = FALSE]
    for (i in 1:floor(obs / 3)) {
      group1 <- group1 + v_2[i, ]^2
      group2 <- group2 + v_2[dim(v_2)[[1]] + 1 - i, ]^2
      Heteroscedasticity[i, ] <- group2 / group1
    }
  }
  ################################ Kalman Smoother ################################

  # Initialise smoothed state and corresponding variance
  a_smooth <- a_fil # N x m
  V <- P_fil # m x m x N

  # Initialise smoothed state disturbance and corresponding variance
  eta <- matrix(0, N, r) # N x r
  eta_var <- array(0, dim = c(r, r, N)) # r x r x N

  # Initialising smoothing error e and corresponding variance D
  # Plus T-statistics for both the observation and state equation
  e <- matrix(0, N, p)
  D <- array(0, dim = c(p, p, N))
  Tstat_observation <- matrix(0, N, p)
  Tstat_state <- matrix(0, N + 1, m)

  # Initialise r and N
  r_UT <- matrix(0, N * p + 1, m)
  N_UT <- array(0, dim = c(m, m, N * p + 1))
  r_vec <- matrix(0, N + 1, m)
  Nmat <- array(0, dim = c(m, m, N + 1))

  # Keep track of which time index and row index to use
  # For selection of system matrices
  t <- N
  row <- p + 1

  # Initialise r_1, N_1, N_2
  r_1 <- matrix(0, m, 1)
  N_1 <- matrix(0, m, m)
  N_2 <- matrix(0, m, m)

  # Normal formulae after the initialisation steps
  for (i in (N * p):1) {

    # Updating indices
    row <- row - 1
    if (row == 0) {
      row <- p
      t <- t - 1
    }

    # When should a transition to the next timepoint be made?
    if (row == 1) {
      timestep <- TRUE

      # T, R, and Q matrices only needed when a transition to the
      # next timepoint is made
      if (is.matrix(T_kal)) {
        T_input <- T_kal
      } else {
        if (t > 1) {
          T_input <- as.matrix(T_kal[, , t - 1])
        }
      }
      if (is.matrix(R_kal)) {
        R_input <- R_kal
      } else {
        if (t > 1) {
          R_input <- matrix(R_kal[, , t - 1], nrow = m)
        }
      }
      if (is.matrix(Q_kal)) {
        Q_input <- Q_kal
      } else {
        if (t > 1) {
          Q_input <- as.matrix(Q_kal[, , t - 1])
        }
      }
    } else {
      timestep <- FALSE
    }

    if (is.matrix(Z_kal)) {
      Z_input <- Z_kal[row, , drop = FALSE]
    } else {
      Z_input <- matrix(Z_kal[row, , t], nrow = 1)
    }

    # For the initialisation steps, other computations are required
    if (i > initialisation_steps) {

      # Check for missing observation
      if (is.na(y[[t, row]])) {
        r_UT[i, ] <- r_UT[i + 1, ]
        N_UT[, , i] <- N_UT[, , i + 1]
      } else {

        # PZ' as in Kalman formulae
        PtZ <- tcrossprod(as.matrix(P_star[, , i]), Z_input)

        # Variance of the current residual/fitted value
        F_scalar <- c(Z_input %*% PtZ)

        # Check if F_scalar is nearly 0
        if (F_scalar < 1e-7) {
          r_UT[i, ] <- r_UT[i + 1, ]
          N_UT[, , i] <- N_UT[, , i + 1]
        } else {

          # Inverse of Fmat
          Finv <- 1 / F_scalar

          # Current residual
          v_UT <- y[[t, row]] - c(Z_input %*% a[i, ])

          # Auxiliary matrices
          K_UT <- PtZ * Finv
          L_UT <- diag(length(Z_input)) - K_UT %*% Z_input

          # r and N
          r_UT[i, ] <- t(Z_input) * Finv * v_UT + crossprod(L_UT, r_UT[i + 1, ])
          N_UT[, , i] <- crossprod(Z_input) * Finv +
            crossprod(L_UT, N_UT[, , i + 1] %*% L_UT)
        }
      }

      # If a transition to the previous timepoint is made,
      # do some additional computations
      if (timestep) {

        # Save r and N for each timepoint and compute smoothed state and variance
        r_vec[t, ] <- r_UT[i, ]
        Nmat[, , t] <- N_UT[, , i]
        a_smooth[t, ] <- a_pred[t, ] + P_pred[, , t] %*% r_vec[t, ]
        V[, , t] <- P_pred[, , t] -
          P_pred[, , t] %*% Nmat[, , t] %*% P_pred[, , t]

        # T-statistic for the state equation
        Tstat_state[t + 1, ] <- r_vec[t + 1, ] /
          sqrt(diag(as.matrix(Nmat[, , t + 1])))

        # Smoothing error e and corresponding variance D
        # Plus T-statistic for the observation equation
        if (is.matrix(Z_kal)) {
          Z_full <- Z_kal
        } else {
          Z_full <- matrix(Z_kal[, , t], nrow = p)
        }
        Finv <- solve(Fmat[, , t])
        K <- T_input %*% P_pred[, , t] %*% crossprod(Z_full, Finv)
        v_0na <- v[t, ]
        v_0na[is.na(v_0na)] <- 0
        e[t, ] <- Finv %*% matrix(v_0na) - crossprod(K, r_vec[t + 1, ])
        e[t, ][is.na(v[t, ])] <- NA
        D[, , t] <- Finv + crossprod(K, Nmat[, , t + 1] %*% K)
        Tstat_observation[t, ] <- e[t, ] / sqrt(diag(as.matrix(D[, , t])))

        # Compute smoothed state disturbance and corresponding variance
        QtR <- tcrossprod(Q_input, R_input)
        eta[t, ] <- QtR %*% r_vec[t + 1, ]
        eta_var[, , t] <- Q_input -
          tcrossprod(QtR %*% Nmat[, , t + 1], QtR)

        # r and N for the next step, not valid/needed for t = 1
        if (t > 1) {
          r_UT[i, ] <- crossprod(T_input, r_UT[i, ])
          N_UT[, , i] <- crossprod(T_input, N_UT[, , i] %*% T_input)
        }
      }
    } else {

      # Check for missing observation
      if (is.na(y[[t, row]])) {
        r_UT[i, ] <- r_UT[i + 1, ]
        N_UT[, , i] <- N_UT[, , i + 1]
      } else {

        # PZ' as in Kalman formulae
        M_inf <- tcrossprod(as.matrix(P_inf[, , i]), Z_input)
        M_star <- tcrossprod(as.matrix(P_star[, , i]), Z_input)

        # Variance matrix of the current residual/fitted value
        F_inf <- c(Z_input %*% M_inf)
        F_star <- c(Z_input %*% M_star)

        # Check if F_inf is nearly 0
        if (F_inf < 1e-7) {

          # Check if F_star is nearly 0
          if (F_star < 1e-7) {
            r_UT[i, ] <- r_UT[i + 1, ]
            N_UT[, , i] <- N_UT[, , i + 1]
          } else {

            # Inverse of Fmat
            F_1 <- 1 / F_star

            # Current residual
            v_UT <- y[[t, row]] - c(Z_input %*% a[i, ])

            # Auxiliary matrices
            K_0 <- M_star * F_1
            L_0 <- diag(length(Z_input)) - K_0 %*% Z_input

            # r and N
            r_UT[i, ] <- t(Z_input) * F_1 * v_UT + crossprod(L_0, r_UT[i + 1, ])
            N_UT[, , i] <- crossprod(Z_input) * F_1 +
              crossprod(L_0, N_UT[, , i + 1] %*% L_0)
            N_1 <- N_1 %*% L_0
          }
        } else {

          # Inverse of Fmat
          F_1 <- 1 / F_inf
          F_2 <- -F_1 * F_star * F_1

          # Current residual
          v_UT <- y[[t, row]] - c(Z_input %*% a[i, ])

          # Auxiliary matrices
          K_0 <- M_inf * F_1
          K_1 <- M_star * F_1 + M_inf * F_2
          L_0 <- diag(length(Z_input)) - K_0 %*% Z_input
          L_1 <- -K_1 %*% Z_input

          # r and N
          r_UT[i, ] <- crossprod(L_0, r_UT[i + 1, ])
          r_1 <- t(Z_input) * F_1 * v_UT +
            crossprod(L_0, r_1) + crossprod(L_1, r_UT[i + 1, ])
          N_UT[, , i] <- crossprod(L_0, N_UT[, , i + 1] %*% L_0)
          tZZ <- crossprod(Z_input)
          N_1 <- tZZ * F_1 + crossprod(L_0, N_1 %*% L_0) +
            crossprod(L_1, N_UT[, , i + 1] %*% L_0) +
            crossprod(L_0, N_UT[, , i + 1] %*% L_1)
          N_2 <- tZZ * F_2 + crossprod(L_0, N_2 %*% L_0) +
            crossprod(L_0, N_1 %*% L_1) +
            crossprod(N_1 %*% L_1, L_0) +
            crossprod(L_1, N_UT[, , i + 1] %*% L_1)
        }
      }

      # If a transition to the previous timepoint is made,
      # do some additional computations
      if (timestep) {

        # Save r and N for each timepoint and compute smoothed state and variance
        r_vec[t, ] <- r_UT[i, ]
        Nmat[, , t] <- N_UT[, , i]
        a_smooth[t, ] <- a_pred[t, ] + P_star[, , i] %*% r_vec[t, ] +
          P_inf[, , i] %*% r_1
        V[, , t] <- P_star[, , i] -
          P_star[, , i] %*% Nmat[, , t] %*% P_star[, , i] -
          t(P_inf[, , i] %*% N_1 %*% P_star[, , i]) -
          P_inf[, , i] %*% N_1 %*% P_star[, , i] -
          P_inf[, , i] %*% N_2 %*% P_inf[, , i]

        # T-statistic for the state equation
        Tstat_state[t + 1, ] <- r_vec[t + 1, ] / sqrt(diag(as.matrix(Nmat[, , t + 1])))

        # Smoothing error e and corresponding variance D
        # Plus T-statistic for the observation equation
        if (is.matrix(Z_kal)) {
          Z_full <- Z_kal
        } else {
          Z_full <- matrix(Z_kal[, , t], nrow = p)
        }
        Finv <- solve(Fmat[, , t])
        K <- T_input %*% P_pred[, , t] %*% crossprod(Z_full, Finv)
        v_0na <- v[t, ]
        v_0na[is.na(v_0na)] <- 0
        e[t, ] <- Finv %*% matrix(v_0na) - crossprod(K, r_vec[t + 1, ])
        e[t, ][is.na(v[t, ])] <- NA
        D[, , t] <- Finv + crossprod(K, Nmat[, , t + 1] %*% K)
        Tstat_observation[t, ] <- e[t, ] / sqrt(diag(as.matrix(D[, , t])))

        # Compute smoothed state disturbance and corresponding variance
        QtR <- tcrossprod(Q_input, R_input)
        eta[t, ] <- QtR %*% r_vec[t + 1, ]
        eta_var[, , t] <- Q_input -
          tcrossprod(QtR %*% Nmat[, , t + 1], QtR)

        # r and N for the next step, not valid/needed for t = 1
        if (t > 1) {
          r_UT[i, ] <- crossprod(T_input, r_UT[i, ])
          r_1 <- crossprod(T_input, r_1)
          N_UT[, , i] <- crossprod(T_input, N_UT[, , i] %*% T_input)
          N_1 <- crossprod(T_input, N_1 %*% T_input)
          N_2 <- crossprod(T_input, N_2 %*% T_input)
        }
      }
    }
  }

  # Smoothed observation disturbance and corresponding variance
  epsilon <- a_smooth[, sys_mat$residuals_state, drop = FALSE]
  epsilon_var <- V[sys_mat$residuals_state, sys_mat$residuals_state, , drop = FALSE]
  ######################################################################################

  ############## Removing residuals from components and storing components #############

  # Removing residuals
  a_pred <- a_pred[, -sys_mat$residuals_state, drop = FALSE]
  P_pred <- P_pred[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  a_fil <- a_fil[, -sys_mat$residuals_state, drop = FALSE]
  P_fil <- P_fil[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_inf <- P_inf[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  P_star <- P_star[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  a_fc <- a_fc[-sys_mat$residuals_state, drop = FALSE]
  P_fc <- P_fc[-sys_mat$residuals_state, -sys_mat$residuals_state, drop = FALSE]
  r_vec <- r_vec[-1, -sys_mat$residuals_state, drop = FALSE]
  Nmat <- Nmat[-sys_mat$residuals_state, -sys_mat$residuals_state, -1, drop = FALSE]
  a_smooth <- a_smooth[, -sys_mat$residuals_state, drop = FALSE]
  V <- V[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  eta <- eta[, -sys_mat$residuals_state, drop = FALSE]
  eta_var <- eta_var[-sys_mat$residuals_state, -sys_mat$residuals_state, , drop = FALSE]
  Tstat_state <- Tstat_state[-1, -sys_mat$residuals_state, drop = FALSE]

  # Storing components
  predicted$yfit <- yfit
  predicted$v <- v
  predicted$Fmat <- Fmat
  predicted$a <- a_pred
  predicted$P <- P_pred
  predicted$P_inf <- P_inf
  predicted$P_star <- P_star
  predicted$a_fc <- a_fc
  predicted$P_fc <- P_fc
  filtered$a <- a_fil
  filtered$P <- P_fil
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
  ######################################################################################

  #### Adjusting dimensions of Z matrices of components ####
  ## -- and adding fitted components of the model ---------##

  # Local Level
  if (local_level_ind & !slope_ind & is.null(level_addvar_list)) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      predicted$level[i, ] <- tempZ %*% matrix(a_pred[i, ])
      filtered$level[i, ] <- tempZ %*% matrix(a_fil[i, ])
      smoothed$level[i, ] <- tempZ %*% matrix(a_smooth[i, ])
    }
    Z_padded$level <- tempZ
  }

  # Local Level + Slope
  if (slope_ind & is.null(level_addvar_list)) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      predicted$level[i, ] <- tempZ %*% matrix(a_pred[i, ])
      filtered$level[i, ] <- tempZ %*% matrix(a_fil[i, ])
      smoothed$level[i, ] <- tempZ %*% matrix(a_smooth[i, ])
    }
    Z_padded$level <- tempZ
  }

  # BSM
  if (length(BSM_vec) > 0) {
    for (s in BSM_vec) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("BSM", s)]])] <- Z_padded[[paste0("BSM", s)]]
      predicted[[paste0("BSM", s)]] <- matrix(0, N, p)
      filtered[[paste0("BSM", s)]] <- matrix(0, N, p)
      smoothed[[paste0("BSM", s)]] <- matrix(0, N, p)
      for (i in 1:N) {
        predicted[[paste0("BSM", s)]][i, ] <- tempZ %*% matrix(a_pred[i, ])
        filtered[[paste0("BSM", s)]][i, ] <- tempZ %*% matrix(a_fil[i, ])
        smoothed[[paste0("BSM", s)]][i, ] <- tempZ %*% matrix(a_smooth[i, ])
      }
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
  if (!is.null(level_addvar_list) & !slope_ind) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      predicted$level[i, ] <- tempZ %*% matrix(a_pred[i, ])
      filtered$level[i, ] <- tempZ %*% matrix(a_fil[i, ])
      smoothed$level[i, ] <- tempZ %*% matrix(a_smooth[i, ])
    }
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
  if (!is.null(level_addvar_list) & slope_ind) {
    tempZ <- matrix(0, p, m - p)
    tempZ[1:length(Z_padded$level)] <- Z_padded$level
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      predicted$level[i, ] <- tempZ %*% matrix(a_pred[i, ])
      filtered$level[i, ] <- tempZ %*% matrix(a_fil[i, ])
      smoothed$level[i, ] <- tempZ %*% matrix(a_smooth[i, ])
    }
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
      predicted[[paste0("Cycle", j)]] <- matrix(0, N, p)
      filtered[[paste0("Cycle", j)]] <- matrix(0, N, p)
      smoothed[[paste0("Cycle", j)]] <- matrix(0, N, p)
      for (i in 1:N) {
        predicted[[paste0("Cycle", j)]][i, ] <- tempZ %*% matrix(a_pred[i, ])
        filtered[[paste0("Cycle", j)]][i, ] <- tempZ %*% matrix(a_fil[i, ])
        smoothed[[paste0("Cycle", j)]][i, ] <- tempZ %*% matrix(a_smooth[i, ])
      }
      Z_padded[[paste0("Cycle", j)]] <- tempZ
    }
  }

  # ARIMA
  if (!is.null(arima_list)) {
    for (j in seq_along(arima_list)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("ARIMA", j)]])] <- Z_padded[[paste0("ARIMA", j)]]
      predicted[[paste0("ARIMA", j)]] <- matrix(0, N, p)
      filtered[[paste0("ARIMA", j)]] <- matrix(0, N, p)
      smoothed[[paste0("ARIMA", j)]] <- matrix(0, N, p)
      for (i in 1:N) {
        predicted[[paste0("ARIMA", j)]][i, ] <- tempZ %*% matrix(a_pred[i, ])
        filtered[[paste0("ARIMA", j)]][i, ] <- tempZ %*% matrix(a_fil[i, ])
        smoothed[[paste0("ARIMA", j)]][i, ] <- tempZ %*% matrix(a_smooth[i, ])
      }
      Z_padded[[paste0("ARIMA", j)]] <- tempZ
    }
  }

  # SARIMA
  if (!is.null(sarima_list)) {
    for (j in seq_along(sarima_list)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded[[paste0("SARIMA", j)]])] <- Z_padded[[paste0("SARIMA", j)]]
      predicted[[paste0("SARIMA", j)]] <- matrix(0, N, p)
      filtered[[paste0("SARIMA", j)]] <- matrix(0, N, p)
      smoothed[[paste0("SARIMA", j)]] <- matrix(0, N, p)
      for (i in 1:N) {
        predicted[[paste0("SARIMA", j)]][i, ] <- tempZ %*% matrix(a_pred[i, ])
        filtered[[paste0("SARIMA", j)]][i, ] <- tempZ %*% matrix(a_fil[i, ])
        smoothed[[paste0("SARIMA", j)]][i, ] <- tempZ %*% matrix(a_smooth[i, ])
      }
      Z_padded[[paste0("SARIMA", j)]] <- tempZ
    }
  }

  # Self Specified
  if (!is.null(self_spec_list)) {
    if (is.matrix(Z_padded$self_spec)) {
      tempZ <- matrix(0, p, m - p)
      tempZ[1:length(Z_padded$self_spec)] <- Z_padded$self_spec
    } else {
      tempZ <- array(0, dim = c(p, m - p, N))
    }
    predicted$self_spec <- matrix(0, N, p)
    filtered$self_spec <- matrix(0, N, p)
    smoothed$self_spec <- matrix(0, N, p)
    for (i in 1:N) {
      if (is.matrix(Z_padded$self_spec)) {
        predicted$self_spec[i, ] <- tempZ %*% matrix(a_pred[i, ])
        filtered$self_spec[i, ] <- tempZ %*% matrix(a_fil[i, ])
        smoothed$self_spec[i, ] <- tempZ %*% matrix(a_smooth[i, ])
      } else {
        tempZ[, , i][1:length(Z_padded$self_spec[, , i])] <- Z_padded$self_spec[, , i]
        predicted$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_pred[i, ])
        filtered$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_fil[i, ])
        smoothed$self_spec[i, ] <- matrix(tempZ[, , i], nrow = p) %*% matrix(a_smooth[i, ])
      }
    }
    Z_padded$self_spec <- tempZ
  }
  ####################################################################################################

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
    list(fit = FALSE, initial = param)
  )
  result$system_matrices <- system_matrices
  result$predicted <- predicted
  result$filtered <- filtered
  result$smoothed <- smoothed
  result$diagnostics <- diagnostics
  return(result)
}
