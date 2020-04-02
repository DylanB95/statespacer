#' State Space Model Evaluation at Specified Parameters
#' 
#' Evaluates the specified State Space model at the parameters specified by the user.
#' 
#' @param param Parameters used to construct the system matrices.
#' @param loglik_only Boolean indicating whether only the loglikelihood should be returned.
#' @inheritParams StateSpaceFit
#' 
#' @return 
#' A list containing:
#' * function_call: A list containing the input to the function.
#' * system_matrices: A list containing the system matrices of the State Space model.
#' * predicted: A list containing the predicted components of the State Space model.
#' * filtered: A list containing the filtered components of the State Space model.
#' * smoothed: A list containing the smoothed components of the State Space model.
#' * diagnostics: A list containing items useful for diagnostical tests.
#' 
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#' @references 
#' \insertRef{durbin2012time}{statespacer}
#' 
#' @examples
#' # Evaluates a local level model for the Nile data at given parameters
#' library(datasets)
#' y <- matrix(Nile)
#' StateSpaceEval(param = c(1, 1), y = y / 100, local_level_ind = TRUE)
#' 
#' @export
#' @importFrom stats qchisq qf
StateSpaceEval <- function(param,
                           y,
                           H_format = NULL,
                           local_level_ind = FALSE,
                           slope_ind = FALSE,
                           BSM_vec = NULL,
                           cycle_ind = FALSE,
                           addvar_list = NULL,
                           level_addvar_list = NULL,
                           slope_addvar_list = NULL,
                           exclude_level = NULL,
                           exclude_slope = NULL,
                           exclude_BSM_list = lapply(BSM_vec, FUN = function(x) {0}),
                           exclude_cycle_list = list(0),
                           damping_factor_ind = rep(TRUE, length(exclude_cycle_list)),
                           format_level = NULL,
                           format_slope = NULL,
                           format_BSM_list = lapply(BSM_vec, FUN = function(x) {NULL}),
                           format_cycle_list = lapply(exclude_cycle_list, FUN = function(x) {NULL}),
                           format_addvar = NULL,
                           format_level_addvar = NULL,
                           loglik_only = FALSE) {

  ##### Initialising lists to return #####
  filtered <- list() # list of items to return that correspond to the Kalman filter
  predicted <- list() # list of items to return that correspond to the Kalman predicter
  smoothed <- list() # list of items to return that correspond to the Kalman smoother
  system_matrices <- list() # list that contains the potentially useful system matrices
  diagnostics <- list() # list that contains the useful items for diagnostic checking
  
  # N = Number of observations
  N <- dim(y)[1]
  
  # p = Number of dependent variables
  p <- dim(y)[2]
  
  # Construct the system matrices
  sys_mat <- GetSysMat(p = p,
                       param = param,
                       update_part = TRUE,
                       H_format = H_format,
                       local_level_ind = local_level_ind,
                       slope_ind = slope_ind,
                       BSM_vec = BSM_vec,
                       cycle_ind = cycle_ind,
                       addvar_list = addvar_list,
                       level_addvar_list = level_addvar_list,
                       slope_addvar_list = slope_addvar_list,
                       exclude_level = exclude_level,
                       exclude_slope = exclude_slope,
                       exclude_BSM_list = exclude_BSM_list,
                       exclude_cycle_list = exclude_cycle_list,
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
  slope_addvar_list <- sys_mat$function_call$slope_addvar_list
  exclude_level <- sys_mat$function_call$exclude_level
  exclude_slope <- sys_mat$function_call$exclude_slope
  exclude_BSM_list <- sys_mat$function_call$exclude_BSM_list
  exclude_cycle_list <- sys_mat$function_call$exclude_cycle_list
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
  H <- sys_mat$H$H
  
  # Dimensions of the system matrices
  Zdim <- length(dim(Z_kal))
  Tdim <- length(dim(T_kal))
  
  # Uncertainty of initial 'guess' of state vector
  kappa <- 1e7
  
  # Number of state parameters
  m <- dim(a)[1]
  
  # Number of disturbances in the state equation
  r <- dim(R_kal)[2]
  
  # Number of diffuse elements
  diffuse_num <- sum(P_inf)
  
  # Check if P_inf is already 0
  Initialisation <- !all(abs(P_inf) < 1e-7)
  
  # Initial P
  P <- kappa * P_inf + P_star
  
  # Initialise (filtered) state and corresponding variance
  a <- t(matrix(a, m, N * p)) # N*p x m
  P_inf <- array(P_inf, dim = c(m, m, N * p)) # m x m x N*p
  P_inf[,,2:(N*p)] <- 0
  P_star <- array(P_star, dim = c(m, m, N * p)) # m x m x N*p
  a_pred <- t(matrix(a, m, N)) # N x m
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
  Initialisation_Steps <- 0
  
  # Keep track of which time index and row index to use
  # For selection of system matrices
  t <- 1
  row <- 0
  
  # Applying Kalman Filter
  for (i in 1:(N*p)) {
    
    # Updating indices
    row <- row + 1
    if (row == (p + 1)) {
      row <- 1
      t <- t + 1
    }
    
    # When should a transition to the next timepoint be made?
    if (i %% p == 0) {
      timestep <- TRUE
      
      # T, R, and Q matrices only needed when a transition to the next timepoint is made
      if (Tdim < 3) {
        T_input <- T_kal
      } else {
        T_input <- as.matrix(T_kal[,,t])
      }
      R_input <- R_kal
      Q_input <- Q_kal
      
    } else {
      timestep <- FALSE
      T_input <- NULL
      R_input <- NULL
      Q_input <- NULL
    }
    
    if (Zdim < 3) {
      Z_input <- matrix(Z_kal[row,], 1, m)
    } else {
      Z_input <- matrix(Z_kal[row,,t], 1, m)
    }
    
    # Apply KalmanEI in initialisation stage, else KalmanUT
    if (Initialisation) {
      
      # Keep track of the number of initialisation steps
      Initialisation_Steps <- Initialisation_Steps + 1
      
      # Calling the Kalman Filter with exact initialisation
      filter_output <- KalmanEI(y = y[t, row],
                                a = as.matrix(a[i,]),
                                P_inf = as.matrix(P_inf[,,i]),
                                P_star = as.matrix(P_star[,,i]), 
                                Z = Z_input, 
                                Tmat = T_input, 
                                R = R_input, 
                                Q = Q_input,
                                timestep = timestep)
      
      # Storing next predicted state and variance used for the next iteration 
      if (i < (N*p)) {
        a[i + 1,] <- filter_output$a
        P_inf[,,i + 1] <- filter_output$P_inf
        P_star[,,i + 1] <- filter_output$P_star
        
        # Check if P_inf converged to zero
        Initialisation <- !all(abs(filter_output$P_inf) < 1e-7)
      }
      
      # Saving predicted and filtered state and corresponding variance for each timestep
      if (timestep & !loglik_only) {
        
        # Predicted state and variance
        if (i < (N*p)) {
          a_pred[t + 1,] <- filter_output$a
          if (Initialisation) {
            P_pred[,,t + 1] <- kappa * filter_output$P_inf + filter_output$P_star
          } else {
            P_pred[,,t + 1] <- filter_output$P_star
          }
          
        } else {
          
          # Store out of sample forecast separately
          a_fc <- filter_output$a
          if (Initialisation) {
            P_fc <- kappa * filter_output$P_inf + filter_output$P_star
          } else {
            P_fc <- filter_output$P_star
          }
        }
        
        # Filtered state and variance
        a_fil[t,] <- filter_output$a_fil
        if (Initialisation) {
          P_fil[,,t] <- kappa * filter_output$P_inf_fil + filter_output$P_star_fil
        } else {
          P_fil[,,t] <- filter_output$P_star_fil
        }
        
        if (Zdim < 3) {
          Z_full <- Z_kal
        } else {
          Z_full <- matrix(Z_kal[,,t], p, m)
        }
        
        # Storing fitted values
        yfit[t,] <- Z_full %*% as.matrix(a_pred[t,])
        
        # Storing residuals
        v[t,] <- y[t,] - yfit[t,]
        
        # Storing Variance - Covariance matrix of fitted values
        Fmat[,,t] <- Z_full %*% as.matrix(P_pred[,,t]) %*% t(Z_full)
      }
      
    } else {
      
      # Calling the Kalman Filter
      filter_output <- KalmanUT(y = y[t, row],
                                a = as.matrix(a[i,]),
                                P = as.matrix(P_star[,,i]), 
                                Z = Z_input, 
                                Tmat = T_input, 
                                R = R_input, 
                                Q = Q_input,
                                timestep = timestep)
      
      # Storing next predicted state and variance used for the next iteration 
      if (i < (N*p)) {
        a[i + 1,] <- filter_output$a
        P_star[,,i + 1] <- filter_output$P
      }
      
      # Saving predicted and filtered state and corresponding variance for each timestep
      if (timestep & !loglik_only) {
        
        # Predicted state and variance
        if (i < (N*p)) {
          a_pred[t + 1,] <- filter_output$a
          P_pred[,,t + 1] <- filter_output$P
          
        } else {
          
          # Store out of sample forecast separately
          a_fc <- filter_output$a
          P_fc <- filter_output$P
        }
        
        # Filtered state and variance
        a_fil[t,] <- filter_output$a_fil
        P_fil[,,t] <- filter_output$P_fil
        
        if (Zdim < 3) {
          Z_full <- Z_kal
        } else {
          Z_full <- matrix(Z_kal[,,t], p, m)
        }
        
        # Storing fitted values
        yfit[t,] <- Z_full %*% as.matrix(a_pred[t,])
        
        # Storing residuals
        v[t,] <- y[t,] - yfit[t,]
        
        # Storing Variance - Covariance matrix of fitted values
        Fmat[,,t] <- Z_full %*% as.matrix(P_pred[,,t]) %*% t(Z_full)
        
        # Inverse of Fmat
        Finv <- solve(Fmat[,,t])
        
        # Singular Value decomposition to obtain square root of the inverse of Fmat
        svd_Finv <- svd(Finv)
        
        # Square root of Inverse of Fmat
        Finv_root <- svd_Finv$u %*% sqrt(diag(svd_Finv$d, p, p)) %*% t(svd_Finv$u)
        
        # Normalised prediction error
        v_norm[t,] <- Finv_root %*% matrix(v[t,])
      }
    }
    
    # Store loglikelihood
    loglik[i] <- filter_output$loglik
  }
  
  # Should only the loglikelihood be returned?
  if (loglik_only) {
    return(sum(loglik, na.rm = TRUE))
  }
  
  # Calculating test statistics based on the normalised prediction errors
  m1 <- apply(v_norm, MARGIN = 2, FUN = mean, na.rm = TRUE)
  v_norm_centered <- t(t(v_norm) - m1)
  m2 <- apply(v_norm_centered^2, MARGIN = 2, FUN = mean, na.rm = TRUE)
  m3 <- apply(v_norm_centered^3, MARGIN = 2, FUN = mean, na.rm = TRUE)
  m4 <- apply(v_norm_centered^4, MARGIN = 2, FUN = mean, na.rm = TRUE)
  Sstat <- m3 / sqrt(m2^3)
  Kstat <- m4 / m2^2
  Nstat <- N * (Sstat^2 / 6 + (Kstat - 3)^2 / 24)
  
  # Correlogram and Box-Ljung statistic
  obs_init <- ceiling(Initialisation_Steps / p) # Number of observations used for initialisation
  obs <- N - obs_init # Number of observations after initialisation
  if (obs > 0) {
    correlogram <- matrix(0, floor(obs/2), p)
    Box_Ljung <- matrix(0, floor(obs/2), p)
    BL_running <- 0
    # Number of non NA values per column
    N_p <- apply(v_norm_centered, MARGIN = 2, FUN = function(x) {sum(!is.na(x))})
    for (i in 1:floor(obs/2)) {
      correlogram[i,] <- apply(
        as.matrix(v_norm_centered[(i+1):N,]) * as.matrix(v_norm_centered[1:(N-i),]), 
        MARGIN = 2, FUN = sum, na.rm = TRUE
      ) / (N_p * m2)
      BL_running <- BL_running + N_p * (N_p + 2) * correlogram[i,]^2 / (N_p - i)
      Box_Ljung[i,] <- BL_running
    }
    
    # Heteroscedasticity test
    Heteroscedasticity <- matrix(0, floor(obs/3), p)
    group1 <- 0
    group2 <- 0
    for (i in 1:floor(obs/3)) {
      group1 <- group1 + v_norm[obs_init + i,]^2
      group2 <- group2 + v_norm[N + 1 - i,]^2
      Heteroscedasticity[i,] <- group2 / group1
    }
  }
  ################################ Kalman Smoother ################################
  
  # Initialise smoothed state and corresponding variance
  a_smooth <- a_fil # N x m
  V <- P_fil # m x m x N
  
  # Initialise smoothed state disturbance and corresponding variance
  eta <- matrix(0, N, r) # N x r
  eta_var <- array(0, dim = c(r, r, N)) # r x r x N
  
  # Initialise r and N
  r_UT <- matrix(0, N*p + 1, m)
  N_UT <- array(0, dim = c(m, m , N*p + 1))
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
  for (i in (N*p):1) {
    
    # Updating indices
    row <- row - 1
    if (row == 0) {
      row <- p
      t <- t - 1
    }
    
    # When should a transition to the next timepoint be made?
    if (row == 1) {
      timestep <- TRUE
      
      # T, R, and Q matrices only needed when a transition to the next timepoint is made
      if (Tdim < 3) {
        T_input <- T_kal
      } else {
        if (t > 1) {
          T_input <- as.matrix(T_kal[,,t - 1])
        }
      }
      R_input <- R_kal
      Q_input <- Q_kal
    } else {
      timestep <- FALSE
    }
    
    if (Zdim < 3) {
      Z_input <- matrix(Z_kal[row,], 1, dim(Z_kal)[2])
    } else {
      Z_input <- matrix(Z_kal[row,,t], 1, dim(Z_kal)[2])
    }
    
    # For the initialisation steps, other computations are required
    if (i > Initialisation_Steps) {
      
      # Check for missing observation
      if (is.na(y[t, row])) {
        r_UT[i,] <- r_UT[i + 1,]
        N_UT[,,i] <- N_UT[,,i + 1]
      } else {
        
        # Auxiliary computations
        Finv <- 1 / c(Z_input %*% as.matrix(P_star[,,i]) %*% t(Z_input))
        v_UT <- y[t, row] - c(Z_input %*% a[i,])
        K_UT <- as.matrix(P_star[,,i]) %*% t(Z_input) * Finv
        L_UT <- diag(length(Z_input)) - K_UT %*% Z_input
        r_UT[i,] <- t(Z_input) * Finv * v_UT + t(L_UT) %*% r_UT[i + 1,]
        N_UT[,,i] <- t(Z_input) %*% Z_input * Finv + t(L_UT) %*% N_UT[,,i + 1] %*% L_UT
      }
      
      # If a transition to the previous timepoint is made, do some additional computations
      if (timestep) {
        
        # Save r and N for each timepoint and compute smoothed state and variance
        r_vec[t,] <- r_UT[i,]
        Nmat[,,t] <- N_UT[,,i]
        a_smooth[t,] <- a_pred[t,] + P_pred[,,t] %*% r_vec[t,]
        V[,,t] <- P_pred[,,t] - P_pred[,,t] %*% Nmat[,,t] %*% P_pred[,,t]
        
        # Compute smoothed state disturbance and corresponding variance
        eta[t,] <- Q_input %*% t(R_input) %*% r_vec[t + 1,]
        eta_var[,,t] <- Q_input - Q_input %*% t(R_input) %*% Nmat[,,t + 1] %*% R_input %*% Q_input
        
        # r and N for the next step, not valid/needed for t = 1
        if (t > 1) {
          r_UT[i,] <- t(T_input) %*% r_UT[i,]
          N_UT[,,i] <- t(T_input) %*% N_UT[,,i] %*% T_input
        }
      }
      
    } else {
      
      # Check for missing observation
      if (is.na(y[t, row])) {
        r_UT[i,] <- r_UT[i + 1,]
        N_UT[,,i] <- N_UT[,,i + 1]
      } else {
        
        # Auxiliary computations
        v_UT <- y[t, row] - c(Z_input %*% a[i,])
        M_inf <- as.matrix(P_inf[,,i]) %*% t(Z_input)
        M_star <- as.matrix(P_star[,,i]) %*% t(Z_input)
        F_inf <- c(Z_input %*% as.matrix(P_inf[,,i]) %*% t(Z_input))
        F_star <- c(Z_input %*% as.matrix(P_star[,,i]) %*% t(Z_input))
        
        # Check if F_inf is nearly 0
        if (F_inf < 1e-7) {
          F_1 <- 1 / F_star
          K_0 <- M_star * F_1
          L_0 <- diag(length(Z_input)) - K_0 %*% Z_input
          r_UT[i,] <- t(Z_input) * F_1 * v_UT + t(L_0) %*% r_UT[i + 1,]
          N_UT[,,i] <- t(Z_input) %*% Z_input * F_1 + t(L_0) %*% N_UT[,,i + 1] %*% L_0
          N_1 <- N_1 %*% L_0
        } else {
          F_1 <- 1 / F_inf
          F_2 <- -F_1 * F_star * F_1
          K_0 <- M_inf * F_1
          K_1 <- M_star * F_1 + M_inf * F_2
          L_0 <- diag(length(Z_input)) - K_0 %*% Z_input
          L_1 <- -K_1 %*% Z_input
          r_UT[i,] <- t(L_0) %*% r_UT[i + 1,]
          r_1 <- t(Z_input) * F_1 * v_UT + t(L_0) %*% r_1 + t(L_1) %*% r_UT[i + 1,]
          N_UT[,,i] <- t(L_0) %*% N_UT[,,i + 1] %*% L_0
          N_1 <- t(Z_input) %*% Z_input * F_1 + t(L_0) %*% N_1 %*% L_0 + t(L_1) %*% N_UT[,,i + 1] %*% L_0 + t(L_0) %*% N_UT[,,i + 1] %*% L_1
          N_2 <- t(Z_input) %*% Z_input * F_2 + t(L_0) %*% N_2 %*% L_0 + t(L_0) %*% N_1 %*% L_1 + t(L_1) %*% t(N_1) %*% L_0 + t(L_1) %*% N_UT[,,i + 1] %*% L_1
        }
      }
      
      # If a transition to the previous timepoint is made, do some additional computations
      if (timestep) {
        
        # Save r and N for each timepoint and compute smoothed state and variance
        r_vec[t, ] <- r_UT[i,]
        Nmat[,,t] <- N_UT[,,i]
        a_smooth[t,] <- a_pred[t,] + P_star[,,i] %*% r_vec[t, ] + P_inf[,,i] %*% r_1
        V[,,t] <- P_star[,,i] - P_star[,,i] %*% Nmat[,,t] %*% P_star[,,i] - t(P_inf[,,i] %*% N_1 %*% P_star[,,i]) - P_inf[,,i] %*% N_1 %*% P_star[,,i] - P_inf[,,i] %*% N_2 %*% P_inf[,,i]
        
        # Compute smoothed state disturbance and corresponding variance
        eta[t,] <- Q_input %*% t(R_input) %*% r_vec[t + 1,]
        eta_var[,,t] <- Q_input - Q_input %*% t(R_input) %*% Nmat[,,t + 1] %*% R_input %*% Q_input
        
        # r and N for the next step, not valid/needed for t = 1
        if (t > 1) {
          r_UT[i,] <- t(T_input) %*% r_UT[i,]
          r_1 <- t(T_input) %*% r_1
          N_UT[,,i] <- t(T_input) %*% N_UT[,,i] %*% T_input
          N_1 <- t(T_input) %*% N_1 %*% T_input
          N_2 <- t(T_input) %*% N_2 %*% T_input
        }
      }
    }
  }
  
  # Smoothed observation disturbance and corresponding variance
  epsilon <- matrix(a_smooth[,sys_mat$residuals_state], N, p)
  epsilon_var <- array(V[sys_mat$residuals_state, sys_mat$residuals_state,], dim = c(p, p, N))
  
  # Calculating smoothing error e and corresponding variance D
  # Plus T-statistics for both the observation and state equation
  e <- matrix(0, N, p)
  D <- array(0, dim = c(p, p, N))
  Tstat_observation <- matrix(0, N, p)
  Tstat_state <- matrix(0, N + 1, m)
  Hinv <- solve(H)
  for (i in 1:N) {
    e[i,] <- Hinv %*% matrix(epsilon[i,])
    D[,,i] <- Hinv %*% (-as.matrix(epsilon_var[,,i]) + H) %*% Hinv
    Tstat_observation[i,] <- e[i,] / sqrt(diag(as.matrix(D[,,i])))
    Tstat_state[i + 1,] <- r_vec[i + 1,] / sqrt(diag(as.matrix(Nmat[,,i + 1])))
  }
  ######################################################################################
  
  ############## Removing residuals from components and storing components #############
  
  # Removing residuals
  a_pred <- matrix(a_pred[, -sys_mat$residuals_state], N, m - p)
  P_pred <- array(P_pred[-sys_mat$residuals_state, -sys_mat$residuals_state,], dim = c(m - p, m - p, N))
  a_fil <- matrix(a_fil[, -sys_mat$residuals_state], N, m - p)
  P_fil <- array(P_fil[-sys_mat$residuals_state, -sys_mat$residuals_state,], dim = c(m - p, m - p, N))
  P_inf <- array(P_inf[-sys_mat$residuals_state, -sys_mat$residuals_state,], dim = c(m - p, m - p, N))
  P_star <- array(P_star[-sys_mat$residuals_state, -sys_mat$residuals_state,], dim = c(m - p, m - p, N))
  a_fc <- matrix(a_fc[-sys_mat$residuals_state])
  P_fc <- as.matrix(P_fc[-sys_mat$residuals_state, -sys_mat$residuals_state])
  r_vec <- matrix(r_vec[-1,-sys_mat$residuals_state], N, m - p)
  Nmat <- array(Nmat[-sys_mat$residuals_state, -sys_mat$residuals_state, -1], dim = c(m - p, m - p, N))
  a_smooth <- matrix(a_smooth[, -sys_mat$residuals_state], N, m - p)
  V <- array(V[-sys_mat$residuals_state, -sys_mat$residuals_state,], dim = c(m - p, m - p, N))
  eta <- matrix(eta[, -sys_mat$r_residuals_state], N, r - p)
  eta_var <- array(eta_var[-sys_mat$r_residuals_state, -sys_mat$r_residuals_state,], dim = c(r - p, r - p, N))
  Tstat_state <- matrix(Tstat_state[-1, -sys_mat$residuals_state], N, m - p)
  
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
  diagnostics$Initialisation_Steps <- Initialisation_Steps
  diagnostics$loglik <- sum(loglik, na.rm = TRUE)
  diagnostics$AIC <- 1/N * (-2 * diagnostics$loglik + 2 * (diffuse_num + length(param)))
  diagnostics$BIC <- 1/N * (-2 * diagnostics$loglik + (diffuse_num + length(param)) * log(N))
  diagnostics$r <- r_vec
  diagnostics$N <- Nmat
  diagnostics$e <- e
  diagnostics$D <- D
  diagnostics$Tstat_observation <- Tstat_observation
  diagnostics$Tstat_state <- Tstat_state
  diagnostics$v_normalised <- v_norm
  diagnostics$Skewness <- Sstat
  diagnostics$Kurtosis <- Kstat
  diagnostics$Normality <- Nstat
  diagnostics$Normality_criticalvalue <- qchisq(0.95, df = 2)
  if (obs > 0) {
    diagnostics$correlogram <- correlogram
    diagnostics$Box_Ljung <- Box_Ljung
    diagnostics$Box_Ljung_criticalvalues <- matrix(qchisq(0.95, df = 1:floor(obs/2)))
    diagnostics$Heteroscedasticity <- Heteroscedasticity
    diagnostics$Heteroscedasticity_criticalvalues <- cbind(
      qf(0.025, df1 = 1:floor(obs/3), df2 = 1:floor(obs/3)), 
      qf(0.975, df1 = 1:floor(obs/3), df2 = 1:floor(obs/3))
    )
  }
  ######################################################################################
  
  #### Adjusting dimensions of Z matrices of components and adding fitted components of the model ####
  
  # Local Level
  if (local_level_ind & !slope_ind & is.null(level_addvar_list) & is.null(slope_addvar_list)) {
    tempZ <- matrix(0, p, m - p)
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      tempZ[1:length(Z_padded$level)] <- Z_padded$level
      predicted$level[i,] <- tempZ %*% as.matrix(a_pred[i,])
      filtered$level[i,] <- tempZ %*% as.matrix(a_fil[i,])
      smoothed$level[i,] <- tempZ %*% as.matrix(a_smooth[i,])
    }
    Z_padded$level <- tempZ
  }
  
  # Local Level + Slope
  if (slope_ind & is.null(level_addvar_list) & is.null(slope_addvar_list)) {
    tempZ <- matrix(0, p, m - p)
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      tempZ[1:length(Z_padded$level)] <- Z_padded$level
      predicted$level[i,] <- tempZ %*% as.matrix(a_pred[i,])
      filtered$level[i,] <- tempZ %*% as.matrix(a_fil[i,])
      smoothed$level[i,] <- tempZ %*% as.matrix(a_smooth[i,])
    }
    Z_padded$level <- tempZ
  }
  
  # BSM
  if (length(BSM_vec) > 0) {
    for (s in BSM_vec) {
      tempZ <- matrix(0, p, m - p)
      predicted[[paste0('BSM', s)]] <- matrix(0, N, p)
      filtered[[paste0('BSM', s)]] <- matrix(0, N, p)
      smoothed[[paste0('BSM', s)]] <- matrix(0, N, p)
      for (i in 1:N) {
        tempZ[1:length(Z_padded[[paste0('BSM', s)]])] <- Z_padded[[paste0('BSM', s)]]
        predicted[[paste0('BSM', s)]][i,] <- tempZ %*% as.matrix(a_pred[i,])
        filtered[[paste0('BSM', s)]][i,] <- tempZ %*% as.matrix(a_fil[i,])
        smoothed[[paste0('BSM', s)]][i,] <- tempZ %*% as.matrix(a_smooth[i,])
      }
      Z_padded[[paste0('BSM', s)]] <- tempZ
    }
  }
  
  # Explanatory variables
  if (!is.null(addvar_list)) {
    tempZ <- array(0, dim = c(p, m - p, N))
    predicted$addvar <- matrix(0, N, p)
    filtered$addvar <- matrix(0, N, p)
    smoothed$addvar <- matrix(0, N, p)
    for (i in 1:N) {
      tempZ[,,i][1:length(Z_padded$addvar[,,i])] <- Z_padded$addvar[,,i]
      predicted$addvar[i,] <- matrix(tempZ[,,i], p, m - p) %*% as.matrix(a_pred[i,])
      filtered$addvar[i,] <- matrix(tempZ[,,i], p, m - p) %*% as.matrix(a_fil[i,])
      smoothed$addvar[i,] <- matrix(tempZ[,,i], p, m - p) %*% as.matrix(a_smooth[i,])
    }
    Z_padded$addvar <- tempZ
    filtered$addvar_coeff <- a_fil[,sys_mat$addvar_state]
    filtered$addvar_coeff_se <- t(apply(P_fil, 3, function(x) {sqrt(diag(as.matrix(x))[sys_mat$addvar_state])}))
    smoothed$addvar_coeff <- a_smooth[,sys_mat$addvar_state]
    smoothed$addvar_coeff_se <- t(apply(V, 3, function(x) {sqrt(diag(as.matrix(x))[sys_mat$addvar_state])}))
  }
  
  # level_addvar
  if (!is.null(level_addvar_list) & is.null(slope_addvar_list) & !slope_ind) {
    tempZ <- matrix(0, p, m - p)
    predicted$level <- matrix(0, N, p)
    filtered$level <- matrix(0, N, p)
    smoothed$level <- matrix(0, N, p)
    for (i in 1:N) {
      tempZ[1:length(Z_padded$level)] <- Z_padded$level
      predicted$level[i,] <- tempZ %*% as.matrix(a_pred[i,])
      filtered$level[i,] <- tempZ %*% as.matrix(a_fil[i,])
      smoothed$level[i,] <- tempZ %*% as.matrix(a_smooth[i,])
    }
    Z_padded$level <- tempZ
    filtered$level_addvar_coeff <- a_fil[,sys_mat$level_addvar_state]
    filtered$level_addvar_coeff_se <- t(apply(P_fil, 3, function(x) {sqrt(diag(as.matrix(x))[sys_mat$level_addvar_state])}))
    smoothed$level_addvar_coeff <- a_smooth[,sys_mat$level_addvar_state]
    smoothed$level_addvar_coeff_se <- t(apply(V, 3, function(x) {sqrt(diag(as.matrix(x))[sys_mat$level_addvar_state])}))
  }
  
  # slope_addvar
  if (!is.null(slope_addvar_list) | (!is.null(level_addvar_list) & slope_ind)) {
    tempZ <- matrix(0, p, m - p)
    predicted$slope <- matrix(0, N, p)
    filtered$slope <- matrix(0, N, p)
    smoothed$slope <- matrix(0, N, p)
    for (i in 1:N) {
      tempZ[1:length(Z_padded$level)] <- Z_padded$level
      predicted$slope[i,] <- tempZ %*% as.matrix(a_pred[i,])
      filtered$slope[i,] <- tempZ %*% as.matrix(a_fil[i,])
      smoothed$slope[i,] <- tempZ %*% as.matrix(a_smooth[i,])
    }
    Z_padded$level <- tempZ
    filtered$level_addvar_coeff <- a_fil[,sys_mat$slope_addvar_state]
    filtered$level_addvar_coeff_se <- t(apply(P_fil, 3, function(x) {sqrt(diag(as.matrix(x))[sys_mat$slope_addvar_state])}))
    smoothed$level_addvar_coeff <- a_smooth[,sys_mat$slope_addvar_state]
    smoothed$level_addvar_coeff_se <- t(apply(V, 3, function(x) {sqrt(diag(as.matrix(x))[sys_mat$slope_addvar_state])}))
  }
  
  # Cycle
  if (cycle_ind) {
    for (j in seq_along(format_cycle_list)) {
      tempZ <- matrix(0, p, m - p)
      predicted[[paste0('Cycle', j)]] <- matrix(0, N, p)
      filtered[[paste0('Cycle', j)]] <- matrix(0, N, p)
      smoothed[[paste0('Cycle', j)]] <- matrix(0, N, p)
      for (i in 1:N) {
        tempZ[1:length(Z_padded[[paste0('Cycle', j)]])] <- Z_padded[[paste0('Cycle', j)]]
        predicted[[paste0('Cycle', j)]][i,] <- tempZ %*% as.matrix(a_pred[i,])
        filtered[[paste0('Cycle', j)]][i,] <- tempZ %*% as.matrix(a_fil[i,])
        smoothed[[paste0('Cycle', j)]][i,] <- tempZ %*% as.matrix(a_smooth[i,])
      }
      Z_padded[[paste0('Cycle', j)]] <- tempZ
    }
  }
  ####################################################################################################
  
  # Filling system_matrices
  system_matrices$H <- sys_mat$H
  system_matrices$Z <- sys_mat$Z
  system_matrices$T <- sys_mat$Tmat
  system_matrices$R <- sys_mat$R
  system_matrices$Q <- sys_mat$Q
  system_matrices$Q_loading_matrix <- sys_mat$Q_loading_matrix
  system_matrices$Q_diagonal_matrix <- sys_mat$Q_diagonal_matrix
  system_matrices$Q_correlation_matrix <- sys_mat$Q_correlation_matrix
  system_matrices$Q_stdev_matrix <- sys_mat$Q_stdev_matrix
  system_matrices$lambda <- sys_mat$lambda
  system_matrices$rho <- sys_mat$rho
  system_matrices$a1 <- sys_mat$a1
  system_matrices$P_inf <- sys_mat$P_inf
  system_matrices$P_star <- sys_mat$P_star
  system_matrices$Z_padded <- Z_padded
  system_matrices$state_label <- sys_mat$state_label
  
  # Returning the result
  result <- list()
  result$function_call <- c(list(param = param, y = y), sys_mat$function_call)
  result$system_matrices <- system_matrices
  result$predicted <- predicted
  result$filtered <- filtered
  result$smoothed <- smoothed
  result$diagnostics <- diagnostics
  return(result)
}