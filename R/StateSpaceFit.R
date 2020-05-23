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
#'   containing vectors of length 3 with the following format: `c(AR, I, MA)`.
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
#' @param self_spec_list A list containing the specification of the self
#'   specified component. See the Details section of `StateSpaceFit()` for
#'   extensive details about the format that must be followed for this
#'   argument.
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
#' @param collapse Boolean indicating whether the observation vector should be
#'   collapsed. Should only be set to `TRUE` if the dimensionality of the
#'   observation vector exceeds the dimensionality of the state vector.
#'   If this is the case, computational gains can be achieved by collapsing
#'   the observation vector.
#' @param standard_errors Boolean indicating whether standard errors should be
#'   computed. \pkg{numDeriv} must be installed in order to compute the
#'   standard errors! Defaults to TRUE if \pkg{numDeriv} is available.
#'
#' @details
#' To fit the specified State Space model, one occasionally has to pay careful
#' attention to the initial values supplied. See
#' \code{vignette("dictionary", "statespacer")} for details.
#' Initial values should not be too large, as some parameters use the
#' transformation exp(2x) to ensure non-negative values, they should also not
#' be too small as some variances might become relatively too close to 0,
#' relative to the magnitude of y.
#'
#' If a component is specified without a `format`, then the format defaults to a
#' diagonal `format`.
#'
#' `self_spec_list` provides a means to incorporate a self-specified component
#' into the State Space model. This argument can only contain any of the
#' following items, of which some are mandatory:
#' * `H_spec`: Boolean indicating whether the H matrix is self-specified.
#'   Should be `TRUE`, if you want to specify the H matrix yourself.
#' * `state_num` (mandatory): The number of state parameters introduced by the
#'   self-specified component. Must be 0 if only `H` is self-specified.
#' * `param_num`: The number of parameters needed by the self-specified
#'   component. Must be specified and greater than 0 if parameters are needed.
#' * `sys_mat_fun`: A function returning a list of system matrices that are
#'   constructed using the parameters. Must have `param` as an argument. The
#'   items in the list returned should have any of the following names: Z,
#'   Tmat, R, Q, a1, P_star, H. Note: Only the system matrices that depend on
#'   the parameters should be returned by the function!
#' * `sys_mat_input`: A list containing additional arguments to `sys_mat_fun`.
#' * `Z`: The Z system matrix if it does not depend on the parameters.
#' * `Tmat`: The T system matrix if it does not depend on the parameters.
#' * `R`: The R system matrix if it does not depend on the parameters.
#' * `Q`: The Q system matrix if it does not depend on the parameters.
#' * `a1`: The initial guess of the state vector. Must be a matrix
#'   with one column.
#' * `P_inf`: The initial diffuse part of the variance - covariance
#'   matrix of the initial state vector. Must be a matrix.
#' * `P_star`: The initial non-diffuse part of the variance - covariance
#'   matrix of the initial state vector if it does not depend on the
#'   parameters. Must be a matrix.
#' * `H`: The H system matrix if it does not depend on the parameters.
#' * `transform_fun`: Function that returns transformed parameters for which
#'   standard errors have to be computed.
#' * `transform_input`: A list containing additional arguments to
#'   `transform_fun.`
#' * `state_only`: The indices of the self specified state that do not play a
#'   role in the observation equations, but only in the state equations. Should
#'   only be used if you want to use `collapse = TRUE` and have some state
#'   parameters that do not play a role in the observation equations. Does not
#'   have to be specified for `collapse = FALSE`.
#'
#' Note: System matrices should only be specified once and need to be
#' specified once! That is, system matrices that are returned by `sys_mat_fun`
#' should not be specified directly, and vice versa. So, system matrices need
#' to be either specified directly, or be returned by `sys_mat_fun`. An
#' exception holds for the case where you \strong{only} want to specify `H`
#' yourself. This will not be checked, so be aware of erroneous output if you
#' do not follow the guidelines of specifying `self_spec_list`. If time-varying
#' system matrices are required, return an array for the time-varying system
#' matrix instead of a matrix.
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
#' * `optim`: A list containing the variables that are returned by the
#'   \code{\link[stats]{optim}} or \code{\link[optimx]{optimr}} function.
#' * `loglik_fun`: Function that returns the loglikelihood of the
#'   specified State Space model, as a function of its parameters.
#' * `standard_errors` (optional): A list containing the standard errors of
#'   the parameters of the State Space model.
#'
#' For extensive details about the object returned,
#' see \code{vignette("dictionary", "statespacer")}.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#'
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
                          format_level_addvar = NULL,
                          method = "BFGS",
                          initial = 0,
                          control = list(),
                          collapse = FALSE,
                          standard_errors = NULL) {

  # Check whether standard_errors should be computed
  if (is.null(standard_errors)) {
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      standard_errors <- TRUE
    } else {
      standard_errors <- FALSE
    }
  } else if (standard_errors) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop(
        "\"numDeriv\" must be installed if standard errors are required.",
        call. = FALSE
      )
    }
  }

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
  p2 <- p

  # Construct the system matrices that do not require any parameters
  sys_mat <- GetSysMat(p = p,
                       param = NULL,
                       update_part = FALSE,
                       add_residuals = !collapse,
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

  # Parameter indices and numbers
  param_indices <- sys_mat$param_indices
  param_num_list <- sys_mat$param_num_list

  # Number of state parameters
  m <- length(sys_mat$state_label)

  if (collapse) {

    # Parameters in the state that will not be used for collapsing
    m2 <- m - length(self_spec_list$state_only)

    # Check if dimensionality of the observation vector is larger than
    # the dimensionality of the state vector
    if (p <= m2) {
      stop(
        paste(
          "`collapse = TRUE` can only be set if the dimensionality of the",
          "observation vector is larger than the dimensionality of the",
          "state vector. Please set `collapse = FALSE`."
        ),
        call. = FALSE
      )
    }

    sys_mat$a_kal <- rbind(matrix(0, m2, 1), sys_mat$a_kal)
    sys_mat$P_inf_kal <- BlockMatrix(matrix(0, m2, m2), sys_mat$P_inf_kal)
    Z_kal_tf <- cbind(diag(1, m2, m2), diag(1, m2, m2),
                      matrix(0, m2, length(self_spec_list$state_only)))
    sys_mat$T_kal <- CombineTRQ(matrix(0, m2, m2), sys_mat$T_kal)
    sys_mat$R_kal <- CombineTRQ(diag(1, m2, m2), sys_mat$R_kal)
    p <- m2

    # Dealing with NA for collapsing transformation
    y_temp <- y
    y_temp[is.na(y)] <- 0
  }
  m <- m + p # Residuals in state

  # Check if P_inf is already 0
  initialisation <- !all(abs(sys_mat$P_inf_kal) < 1e-7)

  # Initialising state vector
  a <- sys_mat$a_kal
  if (length(a) == m) {
    a <- t(matrix(a, m, N * p)) # N*p x m
  }

  # Initialise P_inf
  P_inf <- array(sys_mat$P_inf_kal, dim = c(m, m, N * p)) # m x m x N*p

  # Initialising loglikelihood vector
  loglik <- rep(0, N * p)

  # Initialising Q Matrix that depends on the parameters
  Q_kal <- NULL

  # System matrices of components
  T_list <- sys_mat$Tmat
  R_list <- sys_mat$R
  Q_list <- sys_mat[["Q"]]
  Q_list2 <- sys_mat$Q2
  temp_list <- sys_mat$temp_list
  H <- sys_mat[["H"]][["H"]]
  Z_kal <- sys_mat$Z_kal
  T_kal <- sys_mat$T_kal
  R_kal <- sys_mat$R_kal
  P_star <- sys_mat$P_star_kal

  # Initialise for collapsing
  y_kal <- y
  loglik_add <- 0

  ### Function that returns the LogLikelihood ###
  LogLikelihood <- function(param) {

    ## Constructing Q Matrix ##

    # Local Level
    if (local_level_ind & !slope_ind & is.null(level_addvar_list)) {
      if (param_num_list$level > 0) {
        update <- LocalLevel(p = p2,
                             exclude_level =  exclude_level,
                             fixed_part = FALSE,
                             update_part = TRUE,
                             param = param[param_indices$level],
                             format_level = format_level,
                             decompositions = FALSE
        )
        Q_kal <- BlockMatrix(Q_kal, update[["Q"]])
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$level)
      }
    }

    # Local Level + Slope
    if (slope_ind & is.null(level_addvar_list)) {
      if ((param_num_list$level + param_num_list$slope) > 0) {
        update <- Slope(p = p2,
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
          update <- BSM(p = p2,
                        s = s,
                        exclude_BSM = exclude_BSM_list[[i]],
                        fixed_part = FALSE,
                        update_part = TRUE,
                        param = param[param_indices[[paste0('BSM', s)]]],
                        format_BSM = format_BSM_list[[i]],
                        decompositions = FALSE
          )
          Q_kal <- BlockMatrix(Q_kal, update[["Q"]])
        } else {
          Q_kal <- BlockMatrix(Q_kal, Q_list2[[paste0('BSM', s)]])
        }
      }
    }

    # Explanatory Variables
    if (!is.null(addvar_list)) {
      if (param_num_list$addvar > 0) {
        update <- AddVar(p = p2,
                         addvar_list = addvar_list,
                         fixed_part = FALSE,
                         update_part = TRUE,
                         param = param[param_indices$addvar],
                         format_addvar = format_addvar,
                         decompositions = FALSE
        )
        Q_kal <- BlockMatrix(Q_kal, update[["Q"]])
      } else {
        Q_kal <- BlockMatrix(Q_kal, Q_list$addvar)
      }
    }

    # Local Level + Explanatory Variables
    if (!is.null(level_addvar_list) & !slope_ind) {
      if ((param_num_list$level + param_num_list$level_addvar) > 0) {
        update <- LevelAddVar(p = p2,
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
        update <- SlopeAddVar(p = p2,
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
        update <- Cycle(p = p2,
                        exclude_cycle = exclude_cycle_list[[i]],
                        damping_factor_ind = damping_factor_ind[i],
                        fixed_part = FALSE,
                        update_part = TRUE,
                        param = param[param_indices[[paste0('Cycle', i)]]],
                        format_cycle = format_cycle_list[[i]],
                        decompositions = FALSE
        )
        if (param_num_list[[paste0('Cycle', i)]] > (1 + damping_factor_ind[i])) {
          Q_kal <- BlockMatrix(Q_kal, update[["Q"]])
          if (damping_factor_ind[i]) {
            P_star <- BlockMatrix(P_star, update$P_star)
          }
        } else {
          Q_kal <- BlockMatrix(Q_kal, Q_list2[[paste0('Cycle', i)]])
        }
        T_kal <- CombineTRQ(T_kal, update$Tmat)
      }
    }

    # ARIMA
    if (!is.null(arima_list)) {
      for (i in seq_along(arima_list)) {
        update <- ARIMA(p = p2,
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
        Q_kal <- BlockMatrix(Q_kal, update[["Q"]])
        P_star <- BlockMatrix(P_star, update$P_star)
        if (arima_list[[i]][1] == 0 & arima_list[[i]][3] == 0) {
          T_kal <- CombineTRQ(T_kal, T_list[[paste0('ARIMA', i)]])
          R_kal <- BlockMatrix(R_kal, R_list[[paste0('ARIMA', i)]])
        } else {
          T_kal <- CombineTRQ(T_kal, update$Tmat)
          R_kal <- BlockMatrix(R_kal, update$R)
        }
      }
    }

    # SARIMA
    if (!is.null(sarima_list)) {
      for (i in seq_along(sarima_list)) {
        update <- SARIMA(p = p2,
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
        Q_kal <- BlockMatrix(Q_kal, update[["Q"]])
        P_star <- BlockMatrix(P_star, update$P_star)
        if (sum(sarima_list[[i]]$ar) == 0 & sum(sarima_list[[i]]$ma) == 0) {
          T_kal <- CombineTRQ(T_kal, T_list[[paste0('SARIMA', i)]])
          R_kal <- BlockMatrix(R_kal, R_list[[paste0('SARIMA', i)]])
        } else {
          T_kal <- CombineTRQ(T_kal, update$Tmat)
          R_kal <- BlockMatrix(R_kal, update$R)
        }
      }
    }

    # Self Specified
    if (!is.null(self_spec_list)) {

      # System matrices that depend on parameters
      if (!is.null(self_spec_list$sys_mat_fun)) {
        input <- list(param = param[param_indices$self_spec])
        if (!is.null(self_spec_list$sys_mat_input)) {
          input <- c(input, self_spec_list$sys_mat_input)
        }
        update <- do.call(self_spec_list$sys_mat_fun, input)

        # Adding to full system matrices
        if (!is.null(update$Z)) {
          Z_kal <- CombineZ(Z_kal, update$Z)
        }
        if (!is.null(update$Tmat)) {
          T_kal <- CombineTRQ(T_kal, update$Tmat)
        }
        if (!is.null(update$R)) {
          R_kal <- CombineTRQ(R_kal, update$R)
        }
        if (!is.null(update[["Q"]])) {
          Q_kal <- CombineTRQ(Q_kal, update[["Q"]])
        }
        if (!is.null(update$a1)) {
          a <- rbind(a, update$a1)
          a <- t(matrix(a, m, N * p))
        }
        if (!is.null(update$P_star)) {
          P_star <- BlockMatrix(P_star, update$P_star)
        }
      }

      # System matrices that do not depend on parameters
      if (!is.null(T_list$self_spec)) {
        T_kal <- CombineTRQ(T_kal, T_list$self_spec)
      }
      if (!is.null(R_list$self_spec)) {
        R_kal <- CombineTRQ(R_kal, R_list$self_spec)
      }
      if (!is.null(Q_list$self_spec)) {
        Q_kal <- CombineTRQ(Q_kal, Q_list$self_spec)
      }
      if (!is.null(self_spec_list$P_star)) {
        P_star <- BlockMatrix(P_star, self_spec_list$P_star)
      }

      # Check for state only parameters
      if (!is.null(self_spec_list$state_only) & collapse) {
        state_only_indices <- dim(a)[2] - p - self_spec_list$state_num +
          self_spec_list$state_only
        if (is.matrix(Z_kal)) {
          Z_kal <- Z_kal[, -state_only_indices, drop = FALSE]
        } else {
          Z_kal <- Z_kal[, -state_only_indices,, drop = FALSE]
        }
      }
    }

    # H Matrix
    if (!sys_mat$H_spec) {
      if (param_num_list$H > 0) {
        H <- Cholesky(
          param = param[param_indices$H],
          format = H_format,
          decompositions = FALSE
        )
        if (!collapse) {
          P_star <- BlockMatrix(H, P_star)
        }
      }
    } else if (is.null(self_spec_list[["H"]])) {
      H <- update[["H"]]
      if (!collapse) {
        if (is.matrix(H)) {
          P_star <- BlockMatrix(H, P_star)
        } else {
          P_star <- BlockMatrix(H[,,1], P_star)
        }
      }
    }
    if (!collapse) {
      Q_kal <- CombineTRQ(H, Q_kal)
    }

    # Collapse observation vector
    if (collapse) {
      if (is.matrix(H) & is.matrix(Z_kal)) {
        Hinv <- solve(H)
        ZtHinv <- t(Z_kal) %*% Hinv
        A_star <- solve(ZtHinv %*% Z_kal) %*% ZtHinv
        y_kal <- y_temp %*% t(A_star)
        H_star <- A_star %*% H %*% t(A_star)
        P_star <- BlockMatrix(H_star, P_star)
        Q_kal <- CombineTRQ(H_star, Q_kal)

        # loglikelihood contribution
        for (i in 1:N) {
          e <- y_temp[i,] - Z_kal %*% y_kal[i,]
          loglik_add <- loglik_add - 0.5 * t(e) %*% Hinv %*% e
        }
        y_kal[y_kal == 0] <- NA
        loglik_add <- loglik_add + N * 0.5 * log(det(H_star) / det(H)) -
          0.5 * (sum(!is.na(y)) - sum(!is.na(y_kal))) * log(2 * pi)
      } else if (is.matrix(H) & !is.matrix(Z_kal)) {
        y_kal <- matrix(0, N, m2)
        H_star <- array(0, dim = c(m2, m2, N))
        Hinv <- solve(H)
        detH <- det(H)
        for (i in 1:N) {
          Z_t <- matrix(Z_kal[,,i], nrow = p2)
          ZtHinv <- t(Z_t) %*% Hinv
          A_star <- solve(ZtHinv %*% Z_t) %*% ZtHinv
          y_kal[i,] <- A_star %*% y_temp[i,]
          H_star[,,i] <- A_star %*% H %*% t(A_star)
          e <- y_temp[i,] - Z_t %*% y_kal[i,]
          loglik_add <- loglik_add +
            0.5 * log(det(H_star[,,i]) / detH) -
            0.5 * t(e) %*% Hinv %*% e
        }
        y_kal[y_kal == 0] <- NA
        loglik_add <- loglik_add -
          0.5 * (sum(!is.na(y)) - sum(!is.na(y_kal))) * log(2 * pi)
        P_star <- BlockMatrix(H_star[,,1], P_star)
        Q_kal <- CombineTRQ(H_star, Q_kal)
      } else if (!is.matrix(H) & is.matrix(Z_kal)) {
        y_kal <- matrix(0, N, m2)
        H_star <- array(0, dim = c(m2, m2, N))
        for (i in 1:N) {
          H_t <- as.matrix(H[,,i])
          Hinv <- solve(H_t)
          detH <- det(H_t)
          ZtHinv <- t(Z_kal) %*% Hinv
          A_star <- solve(ZtHinv %*% Z_kal) %*% ZtHinv
          y_kal[i,] <- A_star %*% y_temp[i,]
          H_star[,,i] <- A_star %*% H_t %*% t(A_star)
          e <- y_temp[i,] - Z_kal %*% y_kal[i,]
          loglik_add <- loglik_add +
            0.5 * log(det(H_star[,,i]) / detH) -
            0.5 * t(e) %*% Hinv %*% e
        }
        y_kal[y_kal == 0] <- NA
        loglik_add <- loglik_add -
          0.5 * (sum(!is.na(y)) - sum(!is.na(y_kal))) * log(2 * pi)
        P_star <- BlockMatrix(H_star[,,1], P_star)
        Q_kal <- CombineTRQ(H_star, Q_kal)
      } else {
        y_kal <- matrix(0, N, m2)
        H_star <- array(0, dim = c(m2, m2, N))
        for (i in 1:N) {
          H_t <- as.matrix(H[,,i])
          Hinv <- solve(H_t)
          detH <- det(H_t)
          Z_t <- matrix(Z_kal[,,i], nrow = p2)
          ZtHinv <- t(Z_t) %*% Hinv
          A_star <- solve(ZtHinv %*% Z_t) %*% ZtHinv
          y_kal[i,] <- A_star %*% y_temp[i,]
          H_star[,,i] <- A_star %*% H_t %*% t(A_star)
          e <- y_temp[i,] - Z_t %*% y_kal[i,]
          loglik_add <- loglik_add +
            0.5 * log(det(H_star[,,i]) / detH) -
            0.5 * t(e) %*% Hinv %*% e
        }
        y_kal[y_kal == 0] <- NA
        loglik_add <- loglik_add -
          0.5 * (sum(!is.na(y)) - sum(!is.na(y_kal))) * log(2 * pi)
        P_star <- BlockMatrix(H_star[,,1], P_star)
        Q_kal <- CombineTRQ(H_star, Q_kal)
      }
      Z_kal <- Z_kal_tf
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
        if (is.matrix(T_kal)) {
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

      if (is.matrix(Z_kal)) {
        Z_input <- Z_kal[row,, drop = FALSE]
      } else {
        Z_input <- matrix(Z_kal[row,,t], nrow = 1)
      }

      # Apply KalmanEI in initialisation stage, else KalmanUT
      if (initialisation) {

        # Calling the Kalman Filter with exact initialisation
        filter_output <- KalmanEI(y = y_kal[t, row],
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
        filter_output <- KalmanUT(y = y_kal[t, row],
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
    return(-sum(loglik, na.rm = TRUE) / N - loglik_add / N)
  }

  # Checking the number of initial parameters specified
  if (length(initial) < sys_mat$param_num) {
    warning(
      paste0(
        "Number of initial parameters is less than the required ",
        "amount of parameters (", sys_mat$param_num, "), ",
        "recycling the initial parameters the required amount of times."
      ),
      call. = FALSE
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
      ),
      call. = FALSE
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
                     list(method = method,
                          initial = initial,
                          control = control,
                          collapse = collapse,
                          standard_errors = standard_errors)
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

  # Compute standard errors if required
  if (standard_errors) {

    # Initialise standard_errors list
    standard_errors <- list()

    # Hessian of the loglikelihood evaluated at the ML estimates
    hessian <- numDeriv::hessian(
      func = result$loglik_fun,
      x = fit$par
    )
    result$diagnostics$hessian <- hessian
    min_hessian_inv <- -solve(hessian)

    # Local Level
    if (local_level_ind & !slope_ind & is.null(level_addvar_list)) {
      if (param_num_list$level > 0) {
        TransformFun <- function(param) {
          update <- LocalLevel(p = p2,
                               exclude_level =  exclude_level,
                               fixed_part = FALSE,
                               update_part = TRUE,
                               param = param,
                               format_level = format_level,
                               decompositions = TRUE
          )
          result <- c(update[["Q"]],
                      update$loading_matrix, update$diagonal_matrix,
                      update$correlation_matrix, update$stdev_matrix
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices$level]
        )
        hess_subset <- min_hessian_inv[param_indices$level,
                                       param_indices$level,
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))
        dimQ <- dim(result$system_matrices$Q$level)[1]
        se_index <- 1:(dimQ * dimQ)

        # Q
        standard_errors$Q$level <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_loading_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_loading_matrix$level <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_diagonal_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_diagonal_matrix$level <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_correlation_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_correlation_matrix$level <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_stdev_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_stdev_matrix$level <- matrix(
          std_errors[se_index], dimQ, dimQ
        )
      }
    }

    # Local Level + Slope
    if (slope_ind & is.null(level_addvar_list)) {
      if ((param_num_list$level + param_num_list$slope) > 0) {
        TransformFun <- function(param) {
          update <- Slope(p = p2,
                          exclude_level = exclude_level,
                          exclude_slope = exclude_slope,
                          fixed_part = FALSE,
                          update_part = TRUE,
                          param = param,
                          format_level = format_level,
                          format_slope = format_slope,
                          decompositions = TRUE
          )
          result <- c(update$Q_level,
                      update$loading_matrix_level, update$diagonal_matrix_level,
                      update$correlation_matrix_level, update$stdev_matrix_level,
                      update$Q_slope,
                      update$loading_matrix_slope, update$diagonal_matrix_slope,
                      update$correlation_matrix_slope, update$stdev_matrix_slope
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices$level]
        )
        hess_subset <- min_hessian_inv[param_indices$level,
                                       param_indices$level,
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))

        # level
        if (param_num_list$level > 0) {
          dimQ <- dim(result$system_matrices$Q$level)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_level
          standard_errors$Q$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          std_errors <- std_errors[-se_index]
        }

        # slope
        if (param_num_list$slope > 0) {
          dimQ <- dim(result$system_matrices$Q$slope)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_slope
          standard_errors$Q$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )
        }
      }
    }

    # BSM
    if (length(BSM_vec) > 0) {
      for (i in seq_along(BSM_vec)) {
        s <- BSM_vec[i]
        if (param_num_list[[paste0('BSM', s)]] > 0) {
          TransformFun <- function(param) {
            update <- BSM(p = p2,
                          s = s,
                          exclude_BSM = exclude_BSM_list[[i]],
                          fixed_part = FALSE,
                          update_part = TRUE,
                          transform_only = TRUE,
                          param = param,
                          format_BSM = format_BSM_list[[i]],
                          decompositions = TRUE
            )
            result <- c(update$Q_BSM,
                        update$loading_matrix, update$diagonal_matrix,
                        update$correlation_matrix, update$stdev_matrix
            )
            return(result)
          }
          jacobian <- numDeriv::jacobian(
            func = TransformFun,
            x = fit$par[param_indices[[paste0('BSM', s)]]]
          )
          hess_subset <- min_hessian_inv[param_indices[[paste0('BSM', s)]],
                                         param_indices[[paste0('BSM', s)]],
                                         drop = FALSE]
          std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))
          dimQ <- dim(result$system_matrices$Q[[paste0('BSM', s)]])[1]
          se_index <- 1:(dimQ * dimQ)

          # Q
          standard_errors$Q[[paste0('BSM', s)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix[[paste0('BSM', s)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix[[paste0('BSM', s)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix[[paste0('BSM', s)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix[[paste0('BSM', s)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )
        }
      }
    }

    # Explanatory Variables
    if (!is.null(addvar_list)) {
      if (param_num_list$addvar > 0) {
        TransformFun <- function(param) {
          update <- AddVar(p = p2,
                           addvar_list = addvar_list,
                           fixed_part = FALSE,
                           update_part = TRUE,
                           param = param,
                           format_addvar = format_addvar,
                           decompositions = TRUE
          )
          result <- c(update[["Q"]],
                      update$loading_matrix, update$diagonal_matrix,
                      update$correlation_matrix, update$stdev_matrix
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices$addvar]
        )
        hess_subset <- min_hessian_inv[param_indices$addvar,
                                       param_indices$addvar,
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))
        dimQ <- dim(result$system_matrices$Q$addvar)[1]
        se_index <- 1:(dimQ * dimQ)

        # Q
        standard_errors$Q$addvar <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_loading_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_loading_matrix$addvar <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_diagonal_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_diagonal_matrix$addvar <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_correlation_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_correlation_matrix$addvar <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_stdev_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_stdev_matrix$addvar <- matrix(
          std_errors[se_index], dimQ, dimQ
        )
      }
    }

    # Local Level + Explanatory Variables
    if (!is.null(level_addvar_list) & !slope_ind) {
      if ((param_num_list$level + param_num_list$level_addvar) > 0) {
        TransformFun <- function(param) {
          update <- LevelAddVar(p = p2,
                                exclude_level = exclude_level,
                                level_addvar_list = level_addvar_list,
                                fixed_part = FALSE,
                                update_part = TRUE,
                                param = param,
                                format_level = format_level,
                                format_level_addvar = format_level_addvar,
                                decompositions = TRUE
          )
          result <- c(update$Q_level,
                      update$loading_matrix_level,
                      update$diagonal_matrix_level,
                      update$correlation_matrix_level,
                      update$stdev_matrix_level,
                      update$Q_level_addvar,
                      update$loading_matrix_level_addvar,
                      update$diagonal_matrix_level_addvar,
                      update$correlation_matrix_level_addvar,
                      update$stdev_matrix_level_addvar
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices$level]
        )
        hess_subset <- min_hessian_inv[param_indices$level,
                                       param_indices$level,
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))

        # level
        if (param_num_list$level > 0) {
          dimQ <- dim(result$system_matrices$Q$level)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_level
          standard_errors$Q$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          std_errors <- std_errors[-se_index]
        }

        # level_addvar
        if (param_num_list$level_addvar > 0) {
          dimQ <- dim(result$system_matrices$Q$level_addvar)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_level_addvar
          standard_errors$Q$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )
        }
      }
    }

    # Local Level + Explanatory Variables + Slope
    if (!is.null(level_addvar_list) & slope_ind) {
      if ((param_num_list$level +
           param_num_list$slope +
           param_num_list$level_addvar) > 0) {
        TransformFun <- function(param) {
          update <- SlopeAddVar(p = p2,
                                exclude_level = exclude_level,
                                exclude_slope = exclude_slope,
                                level_addvar_list = level_addvar_list,
                                fixed_part = FALSE,
                                update_part = TRUE,
                                param = param,
                                format_level = format_level,
                                format_slope = format_slope,
                                format_level_addvar = format_level_addvar,
                                decompositions = TRUE
          )
          result <- c(update$Q_level,
                      update$loading_matrix_level,
                      update$diagonal_matrix_level,
                      update$correlation_matrix_level,
                      update$stdev_matrix_level,
                      update$Q_slope,
                      update$loading_matrix_slope,
                      update$diagonal_matrix_slope,
                      update$correlation_matrix_slope,
                      update$stdev_matrix_slope,
                      update$Q_level_addvar,
                      update$loading_matrix_level_addvar,
                      update$diagonal_matrix_level_addvar,
                      update$correlation_matrix_level_addvar,
                      update$stdev_matrix_level_addvar
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices$level]
        )
        hess_subset <- min_hessian_inv[param_indices$level,
                                       param_indices$level,
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))

        # level
        if (param_num_list$level > 0) {
          dimQ <- dim(result$system_matrices$Q$level)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_level
          standard_errors$Q$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_level
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$level <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          std_errors <- std_errors[-se_index]
        }

        # slope
        if (param_num_list$slope > 0) {
          dimQ <- dim(result$system_matrices$Q$slope)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_slope
          standard_errors$Q$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_slope
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$slope <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          std_errors <- std_errors[-se_index]
        }

        # level_addvar
        if (param_num_list$level_addvar > 0) {
          dimQ <- dim(result$system_matrices$Q$level_addvar)[1]
          se_index <- 1:(dimQ * dimQ)

          # Q_level_addvar
          standard_errors$Q$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix_level_addvar
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix$level_addvar <- matrix(
            std_errors[se_index], dimQ, dimQ
          )
        }
      }
    }

    # Cycle
    if (cycle_ind) {
      for (i in seq_along(format_cycle_list)) {
        TransformFun <- function(param) {
          update <- Cycle(p = p2,
                          exclude_cycle = exclude_cycle_list[[i]],
                          damping_factor_ind = damping_factor_ind[i],
                          fixed_part = FALSE,
                          update_part = TRUE,
                          transform_only = TRUE,
                          param = param,
                          format_cycle = format_cycle_list[[i]],
                          decompositions = TRUE
          )
          result <- c(update$lambda, update$rho, update$Q_cycle,
                      update$loading_matrix, update$diagonal_matrix,
                      update$correlation_matrix, update$stdev_matrix
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices[[paste0('Cycle', i)]]]
        )
        hess_subset <- min_hessian_inv[param_indices[[paste0('Cycle', i)]],
                                       param_indices[[paste0('Cycle', i)]],
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))

        # lambda
        standard_errors$lambda[[paste0('Cycle', i)]] <- std_errors[1]
        std_errors <- std_errors[-1]

        # rho
        if (damping_factor_ind[i]) {
          standard_errors$rho[[paste0('Cycle', i)]] <- std_errors[1]
          std_errors <- std_errors[-1]
        }

        if (param_num_list[[paste0('Cycle', i)]] > (1 + damping_factor_ind[i])) {
          dimQ <- dim(result$system_matrices$Q[[paste0('Cycle', i)]])[1]
          se_index <- 1:(dimQ * dimQ)

          # Q
          standard_errors$Q[[paste0('Cycle', i)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_loading_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_loading_matrix[[paste0('Cycle', i)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_diagonal_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_diagonal_matrix[[paste0('Cycle', i)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_correlation_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_correlation_matrix[[paste0('Cycle', i)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )

          # Q_stdev_matrix
          std_errors <- std_errors[-se_index]
          standard_errors$Q_stdev_matrix[[paste0('Cycle', i)]] <- matrix(
            std_errors[se_index], dimQ, dimQ
          )
        }
      }
    }

    # ARIMA
    if (!is.null(arima_list)) {
      for (i in seq_along(arima_list)) {
        TransformFun <- function(param) {
          update <- ARIMA(p = p2,
                          arima_spec = arima_list[[i]],
                          exclude_arima = exclude_arima_list[[i]],
                          fixed_part = FALSE,
                          update_part = TRUE,
                          transform_only = TRUE,
                          param = param,
                          decompositions = TRUE
          )
          result <- c(update$ar, update$ma, update$Q,
                      update$loading_matrix, update$diagonal_matrix,
                      update$correlation_matrix, update$stdev_matrix
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices[[paste0('ARIMA', i)]]]
        )
        hess_subset <- min_hessian_inv[param_indices[[paste0('ARIMA', i)]],
                                       param_indices[[paste0('ARIMA', i)]],
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))

        # AR
        if (!is.null(result$system_matrices$AR[[paste0('ARIMA', i)]])) {
          se_index <- 1:length(result$system_matrices$AR[[paste0('ARIMA', i)]])
          if (is.vector(result$system_matrices$AR[[paste0('ARIMA', i)]])) {
            standard_errors$AR[[paste0('ARIMA', i)]] <- std_errors[se_index]
          } else {
            standard_errors$AR[[paste0('ARIMA', i)]] <- array(
              std_errors[se_index],
              dim = dim(result$system_matrices$AR[[paste0('ARIMA', i)]])
            )
          }
          std_errors <- std_errors[-se_index]
        }

        # MA
        if (!is.null(result$system_matrices$MA[[paste0('ARIMA', i)]])) {
          se_index <- 1:length(result$system_matrices$MA[[paste0('ARIMA', i)]])
          if (is.vector(result$system_matrices$MA[[paste0('ARIMA', i)]])) {
            standard_errors$MA[[paste0('ARIMA', i)]] <- std_errors[se_index]
          } else {
            standard_errors$MA[[paste0('ARIMA', i)]] <- array(
              std_errors[se_index],
              dim = dim(result$system_matrices$MA[[paste0('ARIMA', i)]])
            )
          }
          std_errors <- std_errors[-se_index]
        }

        dimQ <- dim(result$system_matrices$Q[[paste0('ARIMA', i)]])[1]
        se_index <- 1:(dimQ * dimQ)

        # Q
        standard_errors$Q[[paste0('ARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_loading_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_loading_matrix[[paste0('ARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_diagonal_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_diagonal_matrix[[paste0('ARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_correlation_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_correlation_matrix[[paste0('ARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_stdev_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_stdev_matrix[[paste0('ARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )
      }
    }

    # SARIMA
    if (!is.null(sarima_list)) {
      for (i in seq_along(sarima_list)) {
        TransformFun <- function(param) {
          update <- SARIMA(p = p2,
                           sarima_spec = sarima_list[[i]],
                           exclude_sarima = exclude_sarima_list[[i]],
                           fixed_part = FALSE,
                           update_part = TRUE,
                           transform_only = TRUE,
                           param = param,
                           decompositions = TRUE
          )
          result <- c()
          if (length(update$sar) > 0) {
            for (sar in update$sar) {
              result <- c(result, sar)
            }
          }
          if (length(update$sma) > 0) {
            for (sma in update$sma) {
              result <- c(result, sma)
            }
          }
          result <- c(result, update$Q,
                      update$loading_matrix, update$diagonal_matrix,
                      update$correlation_matrix, update$stdev_matrix
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices[[paste0('SARIMA', i)]]]
        )
        hess_subset <- min_hessian_inv[param_indices[[paste0('SARIMA', i)]],
                                       param_indices[[paste0('SARIMA', i)]],
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))

        # SAR
        if (!is.null(result$system_matrices$SAR[[paste0('SARIMA', i)]])) {
          for (j in seq_along(result$system_matrices$SAR[[paste0('SARIMA', i)]])) {
            sar <- result$system_matrices$SAR[[paste0('SARIMA', i)]][[j]]
            se_index <- 1:length(sar)
            if (is.vector(sar)) {
              standard_errors$SAR[[paste0('SARIMA', i)]][[j]] <- std_errors[se_index]
            } else {
              standard_errors$SAR[[paste0('SARIMA', i)]][[j]] <- array(
                std_errors[se_index],
                dim = dim(sar)
              )
            }
            std_errors <- std_errors[-se_index]
          }
          names(standard_errors$SAR[[paste0('SARIMA', i)]]) <- names(
            result$system_matrices$SAR[[paste0('SARIMA', i)]]
          )
        }

        # SMA
        if (!is.null(result$system_matrices$SMA[[paste0('SARIMA', i)]])) {
          for (j in seq_along(result$system_matrices$SMA[[paste0('SARIMA', i)]])) {
            sma <- result$system_matrices$SMA[[paste0('SARIMA', i)]][[j]]
            se_index <- 1:length(sma)
            if (is.vector(sma)) {
              standard_errors$SMA[[paste0('SARIMA', i)]][[j]] <- std_errors[se_index]
            } else {
              standard_errors$SMA[[paste0('SARIMA', i)]][[j]] <- array(
                std_errors[se_index],
                dim = dim(sma)
              )
            }
            std_errors <- std_errors[-se_index]
          }
          names(standard_errors$SMA[[paste0('SARIMA', i)]]) <- names(
            result$system_matrices$SMA[[paste0('SARIMA', i)]]
          )
        }

        dimQ <- dim(result$system_matrices$Q[[paste0('SARIMA', i)]])[1]
        se_index <- 1:(dimQ * dimQ)

        # Q
        standard_errors$Q[[paste0('SARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_loading_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_loading_matrix[[paste0('SARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_diagonal_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_diagonal_matrix[[paste0('SARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_correlation_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_correlation_matrix[[paste0('SARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )

        # Q_stdev_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$Q_stdev_matrix[[paste0('SARIMA', i)]] <- matrix(
          std_errors[se_index], dimQ, dimQ
        )
      }
    }

    # Self Specified
    if (!is.null(self_spec_list$transform_fun)) {
      input <- list(func = self_spec_list$transform_fun,
                    x = fit$par[param_indices$self_spec]
      )
      if (!is.null(self_spec_list$transform_input)) {
        input <- c(input, self_spec_list$transform_input)
      }
      jacobian <- do.call(numDeriv::jacobian, input)
      hess_subset <- min_hessian_inv[param_indices$self_spec,
                                     param_indices$self_spec,
                                     drop = FALSE]
      standard_errors$self_spec <- sqrt(
        diag(jacobian %*% hess_subset %*% t(jacobian))
      )
    }

    # H Matrix
    if (!sys_mat$H_spec) {
      if (param_num_list$H > 0) {
        TransformFun <- function(param) {
          update <- Cholesky(
            param = param,
            format = H_format,
            decompositions = TRUE
          )
          result <- c(update$cov_mat,
                      update$loading_matrix, update$diagonal_matrix,
                      update$correlation_matrix, update$stdev_matrix
          )
          return(result)
        }
        jacobian <- numDeriv::jacobian(
          func = TransformFun,
          x = fit$par[param_indices$H]
        )
        hess_subset <- min_hessian_inv[param_indices$H,
                                       param_indices$H,
                                       drop = FALSE]
        std_errors <- sqrt(diag(jacobian %*% hess_subset %*% t(jacobian)))
        dimH <- dim(result$system_matrices$H$H)[1]
        se_index <- 1:(dimH * dimH)

        # H
        standard_errors$H$H <- matrix(
          std_errors[se_index], dimH, dimH
        )

        # loading_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$H$loading_matrix <- matrix(
          std_errors[se_index], dimH, dimH
        )

        # diagonal_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$H$diagonal_matrix <- matrix(
          std_errors[se_index], dimH, dimH
        )

        # correlation_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$H$correlation_matrix <- matrix(
          std_errors[se_index], dimH, dimH
        )

        # stdev_matrix
        std_errors <- std_errors[-se_index]
        standard_errors$H$stdev_matrix <- matrix(
          std_errors[se_index], dimH, dimH
        )
      }
    }

    # Add standard_errors to result
    result$standard_errors <- standard_errors
  }

  return(result)
}
