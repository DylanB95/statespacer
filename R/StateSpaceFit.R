#' State Space Model Fitting
#'
#' Fits a State Space model as specified by the user.
#'
#' @param y N x p matrix containing the N observations of the p
#'   dependent variables.
#' @param H_format Format of the H system matrix,
#'   the variance - covariance matrix of the observation equation.
#' @param local_level_ind Boolean indicating whether a local level should
#'   be added to the state space model.
#' @param slope_ind Boolean indicating whether a local level + slope should
#'   be added to the state space model.
#' @param BSM_vec Vector containing the BSM seasonalities that have to be added
#'   to the state space model.
#' @param cycle_ind Boolean indicating whether a cycle has to be added
#'   to the state space model.
#' @param addvar_list A list containing the explanatory variables for each of
#'   the dependent variables. The list should contain p (number of dependent
#'   variables) elements. Each element of the list should be a N x k_p matrix,
#'   with k_p being the number of explanatory variables for the pth
#'   dependent variable. If no explanatory variables should be added for one
#'   of the dependent variables, then set the corresponding element to `NULL`.
#' @param level_addvar_list A list containing the explanatory variables for
#'   each of the dependent variables. The list should contain p (number of
#'   dependent variables) elements. Each element of the list should be a
#'   N x k_p matrix, with k_p being the number of explanatory variables
#'   for the pth dependent variable. If no explanatory variables should be
#'   added for one of the dependent variables, then set the corresponding
#'   element to `NULL`.
#' @param arima_list Specifications of the ARIMA components, should be a list
#'   containing vectors of length 3 with the following format: c(AR, I, MA).
#'   Should be a list to allow different ARIMA models for different sets of
#'   dependent variables. Note: The AR and MA coefficients are
#'   constrained such that the AR component is stationary, and the MA
#'   component is invertible.
#'   See \insertCite{ansley1986note;textual}{statespacer} for details about
#'   the transformation used.
#' @param sarima_list Specifications of the SARIMA components, should be a list
#'   containing lists that contain 4 named vectors. Vectors should be named:
#'   "s", "ar", "i", "ma". Should be a list of lists to allow different SARIMA
#'   models for different sets of dependent variables. Note: The AR and MA
#'   coefficients are constrained such that the AR components are stationary,
#'   and the MA components are invertible.
#'   See \insertCite{ansley1986note;textual}{statespacer} for details about
#'   the transformation used. Note: For multivariate models, the order of "s"
#'   matters, as matrix multiplication is not commutative!
#' @param exclude_level Vector containing the dependent variables that should
#'   not get a local level.
#' @param exclude_slope Vector containing the dependent variables that should
#'   not get a slope.
#' @param exclude_BSM_list List of vectors, each vector containing the
#'   dependent variables that should not get the corresponding BSM component.
#' @param exclude_cycle_list The dependent variables that should not get the
#'   corresponding cycle component. Should be a list of vectors to allow
#'   different dependent variables to be excluded for different cycles.
#' @param exclude_arima_list The dependent variables that should not be
#'   involved in the corresponding ARIMA component. Should be a list of
#'   vectors to allow different dependent variables to be excluded for
#'   different ARIMA components.
#' @param exclude_sarima_list The dependent variables that should not be
#'   involved in the corresponding SARIMA component. Should be a list of
#'   vectors to allow different dependent variables to be excluded for
#'   different SARIMA components.
#' @param damping_factor_ind Boolean indicating whether a damping factor should
#'   be included. Must be a vector if multiple cycles are included,
#'   to indicate which cycles should include a damping factor.
#' @param format_level Format of the Q_level system matrix
#'   the variance - covariance matrix of the level state equation.
#' @param format_slope Format of the Q_slope system matrix,
#'   the variance - covariance matrix of the slope state equation.
#' @param format_BSM_list Format of the Q_BSM system matrix,
#'   the variance - covariance matrix of the BSM state equation. Should be a
#'   list to allow different formats for different seasonality periods.
#' @param format_cycle_list Format of the Q_cycle system matrix,
#'   the variance - covariance matrix of the cycle state equation. Should be a
#'   list to allow different formats for different cycles.
#' @param format_addvar Format of the Q_addvar system matrix, the
#'   variance - covariance matrix of the explanatory variables state equation.
#' @param format_level_addvar Format of the Q_level_addvar system matrix, the
#'   variance - covariance matrix of the explanatory variables of the level
#'   state equation.
#' @param method Method that should be used by the \code{\link[stats]{optim}}
#'   or \code{\link[optimx]{optimr}} function to estimate the parameters.
#' @param initial Initial values for the parameter search, allowed to be a
#'   vector or just one number.
#' @param control A list of control parameters for the
#'   \code{\link[stats]{optim}} or \code{\link[optimx]{optimr}} function.
#'
#' @details
#' To fit the specified State Space model, it might be beneficial to scale
#' the dependent variables or pay careful attention to the initial values.
#' If an error occurs, try to scale the dependents by a bigger number, or
#' try different initial values. Initial values should not be too big, as
#' some parameters use the transformation exp(2x) to ensure non-negative
#' values, they should also not be too small as some variances might be
#' relatively close to 0, relative to the magnitude of y.
#' Note: after fitting the model, remember to scale the estimates back!
#' Variances should be multiplied by the square of the scaling number!
#'
#' @return
#' A list containing:
#' * function_call: A list containing the input to the function.
#' * system_matrices: A list containing the system matrices of
#'   the State Space model.
#' * predicted: A list containing the predicted components of
#'   the State Space model.
#' * filtered: A list containing the filtered components of
#'   the State Space model.
#' * smoothed: A list containing the smoothed components of
#'   the State Space model.
#' * diagnostics: A list containing items useful for diagnostical tests.
#' * optim: A list containing the variables that are returned by the
#'   \code{\link[stats]{optim}} or \code{\link[optimx]{optimr}} function.
#' * loglik_fun: Function that returns the loglikelihood of the
#'   specified State Space model, as a function of its parameters.
#' * standard_errors: A list containing the standard errors of
#'   the parameters of the State Space model.
#'
#' For extensive details about the object returned,
#' see \code{vignette("dictionary", package = "statespacer")}.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#' @references
#' \insertRef{durbin2012time}{statespacer}
#'
#' \insertRef{ansley1986note}{statespacer}
#'
#' @examples
#' # Fits a local level model for the Nile data
#' library(datasets)
#' y <- matrix(Nile)
#' fit <- StateSpaceFit(initial = 1, y = y / 100, local_level_ind = TRUE)
#'
#' # Plots the filtered estimates
#' plot(1871:1970, 100 * fit$function_call$y, type = 'p',
#'      ylim = c(500, 1400), xlab = NA, ylab = NA)
#' lines(1871:1970, 100 * fit$filtered$level, type = 'l')
#' lines(1871:1970, 100 * fit$filtered$level +
#'                  1.644854 * 100 * sqrt(fit$filtered$P[1,1,]),
#'       type = 'l', col = 'gray')
#' lines(1871:1970, 100 * fit$filtered$level -
#'                  1.644854 * 100 * sqrt(fit$filtered$P[1,1,]),
#'       type = 'l', col = 'gray')
#'
#' # Plots the smoothed estimates
#' plot(1871:1970, 100 * fit$function_call$y, type = 'p',
#'      ylim = c(500, 1400), xlab = NA, ylab = NA)
#' lines(1871:1970, 100 * fit$smoothed$level,type='l')
#' lines(1871:1970, 100 * fit$smoothed$level +
#'                  1.644854 * 100 * sqrt(fit$smoothed$V[1,1,]),
#'       type = 'l', col = 'gray')
#' lines(1871:1970, 100 * fit$smoothed$level -
#'                  1.644854 * 100 * sqrt(fit$smoothed$V[1,1,]),
#'       type = 'l', col = 'gray')
#'
#' @export
#' @importFrom Rdpack reprompt
StateSpaceFit <- function(y,
                          H_format = NULL,
                          local_level_ind = FALSE,
                          slope_ind = FALSE,
                          BSM_vec = NULL,
                          cycle_ind = FALSE,
                          addvar_list = NULL,
                          level_addvar_list = NULL,
                          arima_list = NULL,
                          sarima_list = NULL,
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
                          format_level_addvar = NULL,
                          method = "BFGS",
                          initial = 0,
                          control = list()) {

  # Check whether optimr is available
  # Note: Negative loglikelihood will be minimised,
  #       equivalent to maximising loglikelihood
  if (requireNamespace("optimx", quietly = TRUE)) {
    optim_fun <- optimx::optimr
    control$maximize <- FALSE
    control$fnscale <- NULL
  } else {
    optim_fun <- stats::optim
    control$fnscale <- 1
    control$maximize <- NULL
  }

  # Default trace = TRUE
  if (is.null(control$trace)) {
    control$trace <- TRUE
  }

  # N = Number of observations
  N <- dim(y)[1]

  # p = Number of dependent variables
  p <- dim(y)[2]

  # Construct the system matrices that do not require any parameters
  sys_mat <- GetSysMat(p = p,
                       param = NULL,
                       update_part = FALSE,
                       H_format = H_format,
                       local_level_ind = local_level_ind,
                       slope_ind = slope_ind,
                       BSM_vec = BSM_vec,
                       cycle_ind = cycle_ind,
                       addvar_list = addvar_list,
                       level_addvar_list = level_addvar_list,
                       arima_list = arima_list,
                       sarima_list = sarima_list,
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

  # Parameter indices and numbers
  param_indices <- sys_mat$param_indices
  param_num_list <- sys_mat$param_num_list

  # Number of state parameters
  m <- dim(sys_mat$a_kal)[1]

  # Initialising state vector
  a <- t(matrix(sys_mat$a_kal, m, N * p)) # N*p x m

  # Check if P_inf is already 0
  initialisation <- !all(abs(sys_mat$P_inf_kal) < 1e-7)

  # Initialise P_inf
  P_inf <- array(sys_mat$P_inf_kal, dim = c(m, m, N * p)) # m x m x N*p

  # Initialising loglikelihood vector
  loglik <- rep(0, N * p)

  # Initialising Q Matrix that depends on the parameters
  Q_kal <- NULL

  # System matrices of components
  T_list <- sys_mat$Tmat
  R_list <- sys_mat$R
  Q_list <- sys_mat$Q
  Q_list2 <- sys_mat$Q2
  temp_list <- sys_mat$temp_list
  H <- sys_mat$H$H
  Z_kal <- sys_mat$Z_kal
  T_kal <- sys_mat$T_kal
  R_kal <- sys_mat$R_kal
  P_star <- sys_mat$P_star_kal

  # Dimensions of the system matrices
  Zdim <- length(dim(Z_kal))
  Tdim <- length(dim(T_kal))

  ### Function that returns the LogLikelihood ###
  LogLikelihood <- function(param) {

    # H Matrix
    if (param_num_list$H > 0) {
      H <- Cholesky(
        param = param[param_indices$H],
        format = H_format,
        decompositions = FALSE
      )

      # Add H matrix to Q matrix
      Q_kal <- BlockMatrix(Q_kal, H)

      # Add H matrix to P_star matrix
      P_star <- BlockMatrix(H, P_star)

    } else {
      Q_kal <- BlockMatrix(Q_kal, H)
    }

    ## Constructing Q Matrix ##

    # Local Level
    if (local_level_ind & !slope_ind & is.null(level_addvar_list)) {
      if (param_num_list$level > 0) {
        update <- LocalLevel(p = p,
                             exclude_level =  exclude_level,
                             fixed_part = FALSE,
                             update_part = TRUE,
                             param = param[param_indices$level],
                             format_level = format_level,
                             decompositions = FALSE
        )
        Q_kal <- BlockMatrix(Q_kal, update$Q)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level)
      }
    }

    # Local Level + Slope
    if (slope_ind & is.null(level_addvar_list)) {
      if ((param_num_list$level + param_num_list$slope) > 0) {
        update <- Slope(p = p,
                        exclude_level = exclude_level,
                        exclude_slope = exclude_slope,
                        fixed_part = FALSE,
                        update_part = TRUE,
                        param = param[param_indices$level],
                        format_level = format_level,
                        format_slope = format_slope,
                        decompositions = FALSE
        )
      }
      if (param_num_list$level > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_level)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level)
      }
      if (param_num_list$slope > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_slope)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$slope)
      }
    }

    # BSM
    if (length(BSM_vec) > 0) {
      for (i in seq_along(BSM_vec)) {
        s <- BSM_vec[i]
        if (param_num_list[[paste0('BSM', s)]] > 0) {
          update <- BSM(p = p,
                        s = s,
                        exclude_BSM = exclude_BSM_list[[i]],
                        fixed_part = FALSE,
                        update_part = TRUE,
                        param = param[param_indices[[paste0('BSM', s)]]],
                        format_BSM = format_BSM_list[[i]],
                        decompositions = FALSE
          )
          Q_kal <- BlockMatrix(Q_kal, update$Q)
        } else {
          Q_kal <- BlockMatrix(Q_kal, Q_list2[[paste0('BSM', s)]])
        }
      }
    }

    # Explanatory Variables
    if (!is.null(addvar_list)) {
      if (param_num_list$addvar > 0) {
        update <- AddVar(p = p,
                         addvar_list = addvar_list,
                         fixed_part = FALSE,
                         update_part = TRUE,
                         param = param[param_indices$addvar],
                         format_addvar = format_addvar,
                         decompositions = FALSE
        )
        Q_kal <- BlockMatrix(Q_kal, update$Q)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$addvar)
      }
    }

    # Local Level + Explanatory Variables
    if (!is.null(level_addvar_list) & !slope_ind) {
      if ((param_num_list$level + param_num_list$level_addvar) > 0) {
        update <- LevelAddVar(p = p,
                              exclude_level = exclude_level,
                              level_addvar_list = level_addvar_list,
                              fixed_part = FALSE,
                              update_part = TRUE,
                              param = param[param_indices$level],
                              format_level = format_level,
                              format_level_addvar = format_level_addvar,
                              decompositions = FALSE
        )
      }
      if (param_num_list$level > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_level)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level)
      }
      if (param_num_list$level_addvar > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_level_addvar)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level_addvar)
      }
    }

    # Local Level + Explanatory Variables + Slope
    if (!is.null(level_addvar_list) & slope_ind) {
      if ((param_num_list$level +
           param_num_list$slope +
           param_num_list$level_addvar) > 0) {
        update <- SlopeAddVar(p = p,
                              exclude_level = exclude_level,
                              exclude_slope = exclude_slope,
                              level_addvar_list = level_addvar_list,
                              fixed_part = FALSE,
                              update_part = TRUE,
                              param = param[param_indices$level],
                              format_level = format_level,
                              format_slope = format_slope,
                              format_level_addvar = format_level_addvar,
                              decompositions = FALSE
        )
      }
      if (param_num_list$level > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_level)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level)
      }
      if (param_num_list$slope > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_slope)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$slope)
      }
      if (param_num_list$level_addvar > 0) {
        Q_kal <- BlockMatrix(Q_kal, update$Q_level_addvar)
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level_addvar)
      }
    }

    # Cycle
    if (cycle_ind) {
      for (i in seq_along(format_cycle_list)) {
        update <- Cycle(p = p,
                        exclude_cycle = exclude_cycle_list[[i]],
                        damping_factor_ind = damping_factor_ind[i],
                        fixed_part = FALSE,
                        update_part = TRUE,
                        param = param[param_indices[[paste0('Cycle', i)]]],
                        format_cycle = format_cycle_list[[i]],
                        decompositions = FALSE
        )
        if (param_num_list[[paste0('Cycle', i)]] > (1 + damping_factor_ind[i])) {
          Q_kal <- BlockMatrix(Q_kal, update$Q)
          if (damping_factor_ind[i]) {
            P_star <- BlockMatrix(P_star, update$P_star)
          }
        } else {
          Q_kal <- BlockMatrix(Q_kal, Q_list2[[paste0('Cycle', i)]])
        }
        if (Tdim < 3) {
          T_kal <- BlockMatrix(T_kal, update$Tmat)
        } else {
          T_kal <- array(
            apply(
              T_kal, 3,
              function(x) BlockMatrix(as.matrix(x), update$Tmat)
            ),
            dim = c(sum(dim(T_kal)[1], dim(update$Tmat)[1]),
                    sum(dim(T_kal)[2], dim(update$Tmat)[2]),
                    N
            )
          )
        }
      }
    }

    # ARIMA
    if (!is.null(arima_list)) {
      for (i in seq_along(arima_list)) {
        update <- ARIMA(p = p,
                        arima_spec = arima_list[[i]],
                        exclude_arima = exclude_arima_list[[i]],
                        fixed_part = FALSE,
                        update_part = TRUE,
                        param = param[param_indices[[paste0('ARIMA', i)]]],
                        decompositions = FALSE,
                        T1 = temp_list[[paste0('ARIMA', i)]]$T1,
                        T2 = temp_list[[paste0('ARIMA', i)]]$T2,
                        T3 = temp_list[[paste0('ARIMA', i)]]$T3,
                        R1 = temp_list[[paste0('ARIMA', i)]]$R1,
                        R2 = temp_list[[paste0('ARIMA', i)]]$R2
        )
        Q_kal <- BlockMatrix(Q_kal, update$Q)
        P_star <- BlockMatrix(P_star, update$P_star)
        if (arima_list[[i]][1] == 0 & arima_list[[i]][3] == 0) {
          if (Tdim < 3) {
            T_kal <- BlockMatrix(T_kal, T_list[[paste0('ARIMA', i)]])
          } else {
            T_kal <- array(
              apply(
                T_kal, 3,
                function(x) BlockMatrix(as.matrix(x), T_list[[paste0('ARIMA', i)]])
              ),
              dim = c(sum(dim(T_kal)[1], dim(T_list[[paste0('ARIMA', i)]])[1]),
                      sum(dim(T_kal)[2], dim(T_list[[paste0('ARIMA', i)]])[2]),
                      N
              )
            )
          }
          R_kal <- BlockMatrix(R_kal, R_list[[paste0('ARIMA', i)]])
        } else {
          if (Tdim < 3) {
            T_kal <- BlockMatrix(T_kal, update$Tmat)
          } else {
            T_kal <- array(
              apply(
                T_kal, 3,
                function(x) BlockMatrix(as.matrix(x), update$Tmat)
              ), dim = c(sum(dim(T_kal)[1], dim(update$Tmat)[1]),
                         sum(dim(T_kal)[2], dim(update$Tmat)[2]),
                         N
                       )
            )
          }
          R_kal <- BlockMatrix(R_kal, update$R)
        }
      }
    }

    # SARIMA
    if (!is.null(sarima_list)) {
      for (i in seq_along(sarima_list)) {
        update <- SARIMA(p = p,
                         sarima_spec = sarima_list[[i]],
                         exclude_sarima = exclude_sarima_list[[i]],
                         fixed_part = FALSE,
                         update_part = TRUE,
                         param = param[param_indices[[paste0('SARIMA', i)]]],
                         decompositions = FALSE,
                         T1 = temp_list[[paste0('SARIMA', i)]]$T1,
                         T2 = temp_list[[paste0('SARIMA', i)]]$T2,
                         T3 = temp_list[[paste0('SARIMA', i)]]$T3,
                         R1 = temp_list[[paste0('SARIMA', i)]]$R1,
                         R2 = temp_list[[paste0('SARIMA', i)]]$R2
        )
        Q_kal <- BlockMatrix(Q_kal, update$Q)
        P_star <- BlockMatrix(P_star, update$P_star)
        if (sum(sarima_list[[i]]$ar) == 0 & sum(sarima_list[[i]]$ma) == 0) {
          if (Tdim < 3) {
            T_kal <- BlockMatrix(T_kal, T_list[[paste0('SARIMA', i)]])
          } else {
            T_kal <- array(
              apply(
                T_kal, 3,
                function(x) BlockMatrix(as.matrix(x), T_list[[paste0('SARIMA', i)]])
              ),
              dim = c(sum(dim(T_kal)[1], dim(T_list[[paste0('SARIMA', i)]])[1]),
                      sum(dim(T_kal)[2], dim(T_list[[paste0('SARIMA', i)]])[2]),
                      N
              )
            )
          }
          R_kal <- BlockMatrix(R_kal, R_list[[paste0('SARIMA', i)]])
        } else {
          if (Tdim < 3) {
            T_kal <- BlockMatrix(T_kal, update$Tmat)
          } else {
            T_kal <- array(
              apply(
                T_kal, 3,
                function(x) BlockMatrix(as.matrix(x), update$Tmat)
              ),
              dim = c(sum(dim(T_kal)[1], dim(update$Tmat)[1]),
                      sum(dim(T_kal)[2], dim(update$Tmat)[2]),
                      N
              )
            )
          }
          R_kal <- BlockMatrix(R_kal, update$R)
        }
      }
    }

    # Initialise P_star
    P_star <- array(P_star, dim = c(m, m, N * p)) # m x m x N*p

    ###### Applying Kalman Filter with exact initialisation ######

    # Keep track of which time index and row index to use
    # For selection of system matrices
    t <- 1
    row <- 0

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

        # T, R, and Q matrices only needed when a transition to
        # the next timepoint is made
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
        Z_input <- Z_kal[row,, drop = FALSE]
      } else {
        Z_input <- matrix(Z_kal[row,,t], nrow = 1)
      }

      # Apply KalmanEI in initialisation stage, else KalmanUT
      if (initialisation) {

        # Calling the Kalman Filter with exact initialisation
        filter_output <- KalmanEI(y = y[t, row],
                                  a = matrix(a[i,]),
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
          initialisation <- !all(abs(filter_output$P_inf) < 1e-7)
        }

      } else {

        # Calling the Kalman Filter
        filter_output <- KalmanUT(y = y[t, row],
                                  a = matrix(a[i,]),
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
      }

      # Store loglikelihood
      loglik[i] <- filter_output$loglik
    }

    # Return the negative average loglikelihood
    # Note: The average is computed as sum / N, using mean would divide by N*p
    return(-sum(loglik, na.rm = TRUE) / N)
  }

  # Checking the number of initial parameters specified
  if (length(initial) < sys_mat$param_num) {
    warning(
      paste0(
        "Number of initial parameters is less than the required ",
        "amount of parameters (", sys_mat$param_num, "), ",
        "recycling the initial parameters the required amount of times."
      )
    )
    initial <- rep(
      initial,
      ceiling(sys_mat$param_num / length(initial))
    )[1:sys_mat$param_num]
  } else if (length(initial) > sys_mat$param_num) {
    warning(
      paste0(
        "Number of initial parameters is greater than the required ",
        "amount of parameters (", sys_mat$param_num, "), ",
        "only using the first ", sys_mat$param_num, " initial parameters."
      )
    )
    initial <- initial[1:sys_mat$param_num]
  }

  # Keeping track of the elapsed time of the optim function
  t1 <- Sys.time()
  message(paste0("Starting the optimisation procedure at: ", t1))

  # Optimising parameters
  fit <- optim_fun(par = initial,
                   fn = LogLikelihood,
                   method = method,
                   control = control
  )

  # Elapsed time
  t2 <- Sys.time()
  message(paste("Finished the optimisation procedure at:", t2))
  message(paste0("Time difference of ", t2 - t1, " ", units(t2 - t1), "\n"))

  # (Adjusted) Input parameters that were passed on to StateSpaceFit
  function_call <- c(list(y = y),
                     sys_mat$function_call,
                     list(method = method, initial = initial, control = control)
  )

  # List that will be returned by the function
  result <- do.call(
    StateSpaceEval,
    c(list(param = fit$par, y = y), sys_mat$function_call)
  )
  result$function_call <- function_call
  result$optim <- fit
  result$loglik_fun <- function(param) -N * LogLikelihood(param)

  # Indices of the parameters for each of the components
  result$diagnostics$param_indices <- param_indices

  # Check if numDeriv is available, and return standard_errors if this is the case
  if (requireNamespace("numDeriv", quietly = TRUE)) {

    # Hessian of the loglikelihood evaluated at the ML estimates of the parameters
    result$diagnostics$hessian <- numDeriv::hessian(
      func = result$loglik_fun,
      x = fit$par
    )

    # Jacobian of the transformed parameters
    jacobian <- do.call(numDeriv::jacobian,
                        c(list(func = TransformParam, x = fit$par, p = p),
                          sys_mat$function_call
                        )
    )

    # Standard errors of the transformed parameters
    std_errors <- sqrt(
      diag(jacobian %*% -solve(result$diagnostics$hessian) %*% t(jacobian))
    )

    # Structured standard errors of the transformed parameters
    result$standard_errors <- StructParam(
      param = std_errors,
      sys_mat = result$system_matrices
    )
  } else{
    warning(
      paste(
        "Install \"numDeriv\" if standard errors of the estimated",
        "parameters are required."
      )
    )
  }

  return(result)
}
