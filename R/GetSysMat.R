#' Construct the System Matrices of a State Space Model
#' 
#' Constructs the system matrices of a State Space model, 
#' as specified by the user.
#' 
#' @param p Number of dependent variables in the model.
#' @param update_part Boolean indicating whether the system matrices should be
#'   constructed that depend on the parameters.
#' @inheritParams StateSpaceFit
#' @inheritParams StateSpaceEval
#' 
#' @return 
#' A list containing the system matrices, among other useful variables.
#'  
#' @noRd
GetSysMat <- function(p,
                      param = NULL,
                      update_part = TRUE,
                      H_format = NULL,
                      local_level_ind = FALSE,
                      slope_ind = FALSE,
                      BSM_vec = NULL,
                      cycle_ind = FALSE,
                      addvar_list = NULL,
                      level_addvar_list = NULL,
                      slope_addvar_list = NULL,
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
                      format_level_addvar = NULL) {

  # Keeping track of the number of state parameters m
  m <- 0
  
  # Keeping track of the number of parameters supplied to the 
  # loglikelihood function
  param_num <- 0
  
  # Keeping track of which parameters to use for the components
  index_par <- 1
  
  #### Adding components to the state space model that are specified ####
  
  # Initialising lists to store matrices of components
  Z_list <- list()
  T_list <- list()
  R_list <- list()
  Q_list <- list()
  Q_list2 <- list() # For full Q matrix, in case it differs from Q_list
  a_list <- list()
  Pinf_list <- list()
  Pstar_list <- list()
  temp_list <- list() # For fixed components of the system matrices
  
  # Initialising lists to store Cholesky decomposition of Q system matrices
  L_list <- list()
  D_list <- list()
  
  # Initialising lists to store correlation and standard deviations of 
  # Q system matrices
  corr_list <- list()
  stdev_list <- list()
  
  # List to store the number of parameters for each of the components
  param_num_list <- list()
  
  # List to store the indices of the parameters for each of the components
  param_indices <- list()
  
  # Initialising label of the state vector
  state_label <- c()
  
  # Initialising list to store the padded versions of the Z system matrix
  Z_padded_list <- list()
  
  # Initialising lists to store lambda and rho parameters of the cycles
  lambda_list <- NULL
  rho_list <- NULL
  
  # Initialising lists to store AR and MA coefficients of the ARIMA components
  ar_list <- NULL
  ma_list <- NULL
  
  # Initialising lists to store AR and MA coefficients of the SARIMA components
  sar_list <- NULL
  sma_list <- NULL
  
  # Initialising list to store state indices of coefficients
  coeff_list <- list()
  
  # Initialising matrices for the Kalman filter
  Z_kal <- NULL
  T_kal <- NULL
  R_kal <- NULL
  Q_kal <- NULL
  a <- NULL
  P_inf <- NULL
  P_star <- NULL
  
  # Dimensions of the matrices start at 2, is overwritten if a component 
  # has dimension 3
  Zdim <- 2 
  Tdim <- 2
  Rdim <- 2 
  Qdim <- 2 
  
  #### Residuals ####
  Z_kal <- cbind(Z_kal, diag(1, p, p))
  T_kal <- BlockMatrix(T_kal, matrix(0, p, p))
  R_kal <- BlockMatrix(R_kal, diag(1, p, p))
  
  # How many parameters are needed for H matrix?
  if (is.null(H_format)) {
    H_format <- diag(1, p, p)
  }
  # Only entries in lower triangle of matrix count
  H_format[upper.tri(H_format)] <- 0
  param_num_list$H <- sum(H_format != 0 & lower.tri(H_format, diag = TRUE))
  param_indices$H <- 1:param_num_list$H
  param_num <- param_num + param_num_list$H
  
  H_list <- list()
  if (update_part) {
    if (param_num_list$H > 0) {
      update <- Cholesky(
        param = param[param_indices$H], 
        format = H_format,
        decompositions = TRUE
      )
      H_list <- list(
        H = update$cov_mat, 
        loading_matrix = update$loading_matrix, 
        diagonal_matrix = update$diagonal_matrix,
        correlation_matrix = update$correlation_matrix,
        stdev_matrix = update$stdev_matrix
      )
      H <- H_list$H
      Q_kal <- BlockMatrix(Q_kal, H)
      P_star <- BlockMatrix(P_star, H)
    } 
  }
  if (param_num_list$H == 0) {
    H <- matrix(0, p, p)
    H_list <- list(H = H)
    Q_kal <- BlockMatrix(Q_kal, H)
    P_star <- BlockMatrix(P_star, H)
  }
  
  a <- rbind(a, matrix(0, p, 1))
  P_inf <- BlockMatrix(P_inf, matrix(0, p, p))
  
  # Indices of the residuals in the state vector
  residuals_state <- 1:p
  
  # Updating indices for the first parameter to use for the next component
  index_par <- index_par + param_num_list$H
  
  #### Local Level ####
  if (local_level_ind & !slope_ind & 
      is.null(level_addvar_list) & is.null(slope_addvar_list)) {
    
    # Constructing matrix with 0s using former m dimension
    zero_mat <- matrix(0, p, m)
    
    # Check which dependent variables will be excluded, removing doubles as well
    exclude_level <- exclude_level[which(exclude_level >= 1 & exclude_level <= p)]
    exclude_level <- unique(exclude_level)
    
    # The number of dependent variables that should get a local level
    n_level <- p - length(exclude_level)
    
    # Saving former m
    m_old <- m
    
    # Updating m
    m <- m + n_level
    
    # Label of the state parameters
    state_label <- c(
      state_label, 
      paste0("Local Level y", (1:p)[which(!1:p %in% exclude_level)])
    )
    
    # How many parameters in param vector are meant for the level component?
    if (is.null(format_level)) {
      format_level <- diag(1, n_level, n_level)
    }
    # Only entries in lower triangle of matrix count
    format_level[upper.tri(format_level)] <- 0 
    param_num_list$level <- sum(format_level != 0 & 
                                lower.tri(format_level, diag = TRUE)
    )
    
    # Last index of parameter that should be used for the component
    index_par2 <- index_par + param_num_list$level - 1
    
    # Indices of parameters that should be used for the component
    param_indices$level <- index_par:index_par2
    
    # Keeping track of how many parameters the State Space model uses
    param_num <- param_num + param_num_list$level
    
    # Calling the proper function to obtain system matrices
    update <- LocalLevel(p = p,
                         exclude_level = exclude_level,
                         fixed_part = TRUE,
                         update_part = update_part,
                         param = param[param_indices$level],
                         format_level = format_level,
                         decompositions = TRUE
    )
    
    # Updating indices for the first parameter to use for the next component
    index_par <- index_par2 + 1
    
    # Storing system matrices
    Z_list$level <- update$Z
    Z_kal <- cbind(Z_kal, update$Z)
    T_list$level <- update$Tmat
    T_kal <- BlockMatrix(T_kal, update$Tmat)
    R_list$level <- update$R
    R_kal <- BlockMatrix(R_kal, update$R)
    a_list$level <- update$a1
    a <- rbind(a, update$a1)
    Pinf_list$level <- update$P_inf
    P_inf <- BlockMatrix(P_inf, update$P_inf)
    Pstar_list$level <- update$P_star
    P_star <- BlockMatrix(P_star, update$P_star)
    Z_padded_list$level <- cbind(zero_mat, update$Z)
    if (param_num_list$level == 0 | update_part) {
      Q_list$level <- update$Q
      Q_kal <- BlockMatrix(Q_kal, update$Q)
      L_list$level <- update$loading_matrix
      D_list$level <- update$diagonal_matrix
      corr_list$level <- update$correlation_matrix
      stdev_list$level <- update$stdev_matrix
    }
  }
  
  #### Local Level + Slope ####
  if (slope_ind & is.null(level_addvar_list) & is.null(slope_addvar_list)) {
    
    # Constructing matrix with 0s using former m dimension
    zero_mat <- matrix(0, p, m)
    
    # Check which dependent variables will be excluded, removing doubles as well
    exclude_level <- exclude_level[which(exclude_level >= 1 & exclude_level <= p)]
    exclude_level <- unique(exclude_level)
    
    # The number of dependent variables that should get a local level
    n_level <- p - length(exclude_level)
    
    # Check which dependent variables will not get a slope, removing doubles as well
    exclude_slope <- exclude_slope[which(exclude_slope >= 1 & exclude_slope <= p)]
    exclude_slope <- c(exclude_level, exclude_slope)
    exclude_slope <- unique(exclude_slope)
    
    # The number of local levels that should get a slope
    n_slope <- p - length(exclude_slope)
    
    # Saving former m
    m_old <- m
    
    # Updating m
    m <- m + n_level + n_slope
    
    # Label of the state parameters
    state_label <- c(
      state_label, 
      paste0("Local Level y", (1:p)[which(!1:p %in% exclude_level)])
    )
    state_label <- c(
      state_label, 
      paste0("Slope y", (1:p)[which(!1:p %in% exclude_slope)])
    )
    
    # How many parameters in param vector are meant for the level component?
    if (is.null(format_level)) {
      format_level <- diag(1, n_level, n_level)
    }
    # Only entries in lower triangle of matrix count
    format_level[upper.tri(format_level)] <- 0 
    param_num_list$level <- sum(format_level != 0 & 
                                lower.tri(format_level, diag = TRUE)
    )
    
    # How many parameters in param vector are meant for the slope component?
    if (is.null(format_slope)) {
      format_slope <- diag(1, n_slope, n_slope)
    }
    # Only entries in lower triangle of matrix count
    format_slope[upper.tri(format_slope)] <- 0 
    param_num_list$slope <- sum(format_slope != 0 & 
                                lower.tri(format_slope, diag = TRUE)
    )
    
    # Last index of parameter that should be used for the component
    index_par2 <- index_par + param_num_list$level + param_num_list$slope - 1
    
    # Indices of parameters that should be used for the component
    param_indices$slope <- index_par:index_par2
    
    # Keeping track of how many parameters the State Space model needs
    param_num <- param_num + param_num_list$level + param_num_list$slope
    
    # Calling the proper function to obtain system matrices
    update <- Slope(p = p,
                    exclude_level = exclude_level,
                    exclude_slope = exclude_slope,
                    fixed_part = TRUE,
                    update_part = update_part,
                    param = param[param_indices$slope],
                    format_level = format_level,
                    format_slope = format_slope,
                    decompositions = TRUE
    )
    
    # Updating indices for the first parameter to use for the next component
    index_par <- index_par2 + 1
    
    # Storing system matrices
    Z_list$slope <- update$Z
    Z_kal <- cbind(Z_kal, update$Z)
    T_list$slope <- update$Tmat
    T_kal <- BlockMatrix(T_kal, update$Tmat)
    R_list$slope <- update$R
    R_kal <- BlockMatrix(R_kal, update$R)
    a_list$slope <- update$a1
    a <- rbind(a, update$a1)
    Pinf_list$slope <- update$P_inf
    P_inf <- BlockMatrix(P_inf, update$P_inf)
    Pstar_list$slope <- update$P_star
    P_star <- BlockMatrix(P_star, update$P_star)
    Z_padded_list$level <- cbind(zero_mat, update$Z)
    if (param_num_list$level == 0 | update_part) {
      Q_list$level <- update$Q_level
      Q_kal <- BlockMatrix(Q_kal, update$Q_level)
      L_list$level <- update$loading_matrix_level
      D_list$level <- update$diagonal_matrix_level
      corr_list$level <- update$correlation_matrix_level
      stdev_list$level <- update$stdev_matrix_level
    }
    if (param_num_list$slope == 0 | update_part) {
      Q_list$slope <- update$Q_slope
      Q_kal <- BlockMatrix(Q_kal, update$Q_slope)
      L_list$slope <- update$loading_matrix_slope
      D_list$slope <- update$diagonal_matrix_slope
      corr_list$slope <- update$correlation_matrix_slope
      stdev_list$slope <- update$stdev_matrix_slope
    }
  }
  
  #### BSM ####
  if (length(BSM_vec) > 0) {
    for (i in seq_along(BSM_vec)) {
      
      # Constructing matrix with 0s using former m dimension
      zero_mat <- matrix(0, p, m)
      
      # Period of the season
      s <- BSM_vec[i]
      
      # Check which dependent variables will be excluded, removing doubles as well
      exclude_BSM_list[[i]] <- exclude_BSM_list[[i]][
        which(exclude_BSM_list[[i]] >= 1 & exclude_BSM_list[[i]] <= p)
      ]
      exclude_BSM_list[[i]] <- unique(exclude_BSM_list[[i]])
      
      # The number of dependent variables that should get a BSM`s` component
      n_BSM <- p - length(exclude_BSM_list[[i]])
      
      # Saving former m
      m_old <- m
      
      # Updating m
      m <- m + (s - 1) * n_BSM
      
      # Label of the state parameters
      state_label <- c(
        state_label, 
        rep(
          paste0("BSM", s, " y", (1:p)[which(!1:p %in% exclude_BSM_list[[i]])]),
          s - 1
        )
      )
      
      # How many parameters in param vector are meant for the BSM`s` component?
      if (is.null(format_BSM_list[[i]])) {
        format_BSM_list[[i]] <- diag(1, n_BSM, n_BSM)
      }
      # Only entries in lower triangle of matrix count
      format_BSM_list[[i]][upper.tri(format_BSM_list[[i]])] <- 0 
      param_num_list[[paste0('BSM', s)]] <- sum(
        format_BSM_list[[i]] != 0 & 
        lower.tri(format_BSM_list[[i]], diag = TRUE)
      )
      
      # Last index of parameter that should be used for the component
      index_par2 <- index_par + param_num_list[[paste0('BSM', s)]] - 1
      
      # Indices of parameters that should be used for the component
      param_indices[[paste0('BSM', s)]] <- index_par:index_par2
      
      # Keeping track of how many parameters the State Space model needs
      param_num <- param_num + param_num_list[[paste0('BSM', s)]]
      
      # Calling the proper function to obtain system matrices
      update <- BSM(p = p, 
                    s = s, 
                    exclude_BSM = exclude_BSM_list[[i]],
                    fixed_part = TRUE,
                    update_part = update_part,
                    param = param[param_indices[[paste0('BSM', s)]]],
                    format_BSM = format_BSM_list[[i]],
                    decompositions = TRUE
      )
      
      # Updating indices for the first parameter to use for the next component
      index_par <- index_par2 + 1
      
      # Storing system matrices
      Z_list[[paste0('BSM', s)]] <- update$Z
      Z_kal <- cbind(Z_kal, update$Z)
      T_list[[paste0('BSM', s)]] <- update$Tmat
      T_kal <- BlockMatrix(T_kal, update$Tmat)
      R_list[[paste0('BSM', s)]] <- update$R
      R_kal <- BlockMatrix(R_kal, update$R)
      a_list[[paste0('BSM', s)]] <- update$a1
      a <- rbind(a, update$a1)
      Pinf_list[[paste0('BSM', s)]] <- update$P_inf
      P_inf <- BlockMatrix(P_inf, update$P_inf)
      Pstar_list[[paste0('BSM', s)]] <- update$P_star
      P_star <- BlockMatrix(P_star, update$P_star)
      Z_padded_list[[paste0('BSM', s)]] <- cbind(zero_mat, update$Z)
      if (param_num_list[[paste0('BSM', s)]] == 0 | update_part) {
        Q_list[[paste0('BSM', s)]] <- update$Q_BSM
        Q_list2[[paste0('BSM', s)]] <- update$Q
        Q_kal <- BlockMatrix(Q_kal, update$Q)
        L_list[[paste0('BSM', s)]] <- update$loading_matrix
        D_list[[paste0('BSM', s)]] <- update$diagonal_matrix
        corr_list[[paste0('BSM', s)]] <- update$correlation_matrix
        stdev_list[[paste0('BSM', s)]] <- update$stdev_matrix
      }
    }
  }
  
  #### Explanatory Variables ####
  if (!is.null(addvar_list)) {
    
    # Constructing matrix with 0s using former m dimension
    zero_mat <- matrix(0, p, m)
    
    # Number of coefficients
    k <- sum(
      sapply(addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[2]}})
    )
    
    # Number of observations
    N <- max(
      sapply(addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[1]}})
    )
    
    # Saving former m
    m_old <- m
    
    # Updating m
    m <- m + k
    
    # Indices of the coefficients in the state vector
    coeff_list$addvar_state <- (m_old + 1):m
    
    # Function for naming the explanatory variables
    names_fun <- function(i) {
      if (!is.null(addvar_list[[i]]) & is.null(colnames(addvar_list[[i]]))) {
        return(
          rep(paste0("Explanatory Variable of y", i), dim(addvar_list[[i]])[2])
        )
      } else {
        return(colnames(addvar_list[[i]]))
      }
    }
    
    # Label of the state parameters
    state_label <- c(state_label, do.call(c, lapply(1:p, names_fun)))
    
    # How many parameters in param vector are meant for the explanatory variables?
    if (is.null(format_addvar)) {
      format_addvar <- matrix(0, k, k)
    }
    # Only entries in lower triangle of matrix count
    format_addvar[upper.tri(format_addvar)] <- 0 
    param_num_list$addvar <- sum(
      format_addvar != 0 & 
      lower.tri(format_addvar, diag = TRUE)
    )
    
    # Last index of parameter that should be used for the component
    index_par2 <- index_par + param_num_list$addvar - 1
    
    # Indices of parameters that should be used for the component
    param_indices$addvar <- index_par:index_par2
    
    # Keeping track of how many parameters the State Space model needs
    param_num <- param_num + param_num_list$addvar
    
    # Calling the proper function to obtain system matrices
    update <- AddVar(p = p, 
                     addvar_list = addvar_list,
                     fixed_part = TRUE,
                     update_part = update_part,
                     param = param[param_indices$addvar],
                     format_addvar = format_addvar,
                     decompositions = TRUE
    )
    
    # Updating indices for the first parameter to use for the next component
    index_par <- index_par2 + 1
    
    # Storing system matrices
    Z_list$addvar <- update$Z
    Zdim <- 3
    Z_kal <- array(
      apply(
        update$Z, 3, 
        function(x) {cbind(
          Z_kal, 
          matrix(x, p, k)
        )}
      ), 
      dim = c(p, sum(dim(Z_kal)[2], k), N)
    )
    T_list$addvar <- update$Tmat
    T_kal <- BlockMatrix(T_kal, update$Tmat)
    R_list$addvar <- update$R
    R_kal <- BlockMatrix(R_kal, update$R)
    a_list$addvar <- update$a1
    a <- rbind(a, update$a1)
    Pinf_list$addvar <- update$P_inf
    P_inf <- BlockMatrix(P_inf, update$P_inf)
    Pstar_list$addvar <- update$P_star
    P_star <- BlockMatrix(P_star, update$P_star)
    Z_padded_list$addvar <- array(
      apply(
        update$Z, 3, 
        function(x) {cbind(
          zero_mat, 
          matrix(x, p, k)
        )}
      ), 
      dim = c(p, sum(dim(zero_mat)[2], k), N)
    )
    if (param_num_list$addvar == 0 | update_part) {
      Q_list$addvar <- update$Q
      Q_kal <- BlockMatrix(Q_kal, update$Q)
      L_list$addvar <- update$loading_matrix
      D_list$addvar <- update$diagonal_matrix
      corr_list$addvar <- update$correlation_matrix
      stdev_list$addvar <- update$stdev_matrix
    }
  }
  
  #### Local Level + Explanatory Variables ####
  if (!is.null(level_addvar_list) & is.null(slope_addvar_list) & !slope_ind) {
    
    # Constructing matrix with 0s using former m dimension
    zero_mat <- matrix(0, p, m)
    
    # Check which dependent variables will be excluded, removing doubles as well
    exclude_level <- exclude_level[which(exclude_level >= 1 & exclude_level <= p)]
    exclude_level <- unique(exclude_level)
    
    # The number of dependent variables that should get a local level
    n_level <- p - length(exclude_level)
    
    # Number of coefficients
    k_level <- sum(
      sapply(level_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[2]}})
    )
    
    # Number of observations
    N <- max(
      sapply(level_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[1]}})
    )
    
    # Saving former m
    m_old <- m
    
    # Updating m
    m <- m + n_level + k_level
    
    # Indices of the coefficients in the state vector
    coeff_list$level_addvar_state <- (m_old + n_level + 1):m
    
    # Function for naming the explanatory variables
    names_fun <- function(i) {
      if (!is.null(level_addvar_list[[i]]) & 
          is.null(colnames(level_addvar_list[[i]]))) {
        return(
          rep(
            paste0("Explanatory Variable in level y", i), 
            dim(level_addvar_list[[i]])[2]
          )
        )
      } else {
        return(colnames(level_addvar_list[[i]]))
      }
    }
    
    # Label of the state parameters
    state_label <- c(
      state_label, 
      paste0("Local Level y", (1:p)[which(!1:p %in% exclude_level)])
    )
    state_label <- c(state_label, do.call(c, lapply(1:p, names_fun)))
    
    # How many parameters in param vector are meant for the level component?
    if (is.null(format_level)) {
      format_level <- diag(1, n_level, n_level)
    }
    # Only entries in lower triangle of matrix count
    format_level[upper.tri(format_level)] <- 0 
    param_num_list$level <- sum(
      format_level != 0 & 
      lower.tri(format_level, diag = TRUE)
    )
    
    # How many parameters in param vector are meant for the explanatory variables?
    if (is.null(format_level_addvar)) {
      format_level_addvar <- matrix(0, k_level, k_level)
    }
    # Only entries in lower triangle of matrix count
    format_level_addvar[upper.tri(format_level_addvar)] <- 0 
    param_num_list$level_addvar <- sum(
      format_level_addvar != 0 & 
      lower.tri(format_level_addvar, diag = TRUE)
    )
    
    # Last index of parameter that should be used for the component
    index_par2 <- index_par + param_num_list$level + 
      param_num_list$level_addvar - 1
    
    # Indices of parameters that should be used for the component
    param_indices$level_addvar <- index_par:index_par2
    
    # Keeping track of how many parameters the State Space model needs
    param_num <- param_num + param_num_list$level + param_num_list$level_addvar
    
    # Calling the proper function to obtain system matrices
    update <- LevelAddVar(p = p,
                          exclude_level = exclude_level,
                          level_addvar_list = level_addvar_list,
                          fixed_part = TRUE,
                          update_part = update_part,
                          param = param[param_indices$level_addvar],
                          format_level = format_level,
                          format_level_addvar = format_level_addvar,
                          decompositions = TRUE
    )
    
    # Updating indices for the first parameter to use for the next component
    index_par <- index_par2 + 1
    
    # Storing system matrices
    Z_list$level_addvar <- update$Z
    if (Zdim < 3) {
      Z_kal <- cbind(Z_kal, update$Z)
    } else {
      Z_kal <- array(
        apply(
          Z_kal, 3, 
          function(x) {cbind(
            matrix(x, p, dim(Z_kal)[2]), 
            update$Z
          )}
        ), 
        dim = c(p, sum(dim(Z_kal)[2], k_level), N)
      )
    }
    T_list$level_addvar <- update$Tmat
    Tdim <- 3
    T_kal <- array(
      apply(
        update$Tmat, 3, 
        function(x) BlockMatrix(T_kal, as.matrix(x))
      ), 
      dim = c(sum(dim(T_kal)[1], k_level), sum(dim(T_kal)[2], k_level), N)
    )
    R_list$level_addvar <- update$R
    R_kal <- BlockMatrix(R_kal, update$R)
    a_list$level_addvar <- update$a1
    a <- rbind(a, update$a1)
    Pinf_list$level_addvar <- update$P_inf
    P_inf <- BlockMatrix(P_inf, update$P_inf)
    Pstar_list$level_addvar <- update$P_star
    P_star <- BlockMatrix(P_star, update$P_star)
    Z_padded_list$level <- cbind(zero_mat, update$Z)
    if (param_num_list$level == 0 | update_part) {
      Q_list$level <- update$Q_level
      Q_kal <- BlockMatrix(Q_kal, update$Q_level)
      L_list$level <- update$loading_matrix_level
      D_list$level <- update$diagonal_matrix_level
      corr_list$level <- update$correlation_matrix_level
      stdev_list$level <- update$stdev_matrix_level
    }
    if (param_num_list$level_addvar == 0 | update_part) {
      Q_list$level_addvar <- update$Q_level_addvar
      Q_kal <- BlockMatrix(Q_kal, update$Q_level_addvar)
      L_list$level_addvar <- update$loading_matrix_level_addvar
      D_list$level_addvar <- update$diagonal_matrix_level_addvar
      corr_list$level_addvar <- update$correlation_matrix_level_addvar
      stdev_list$level_addvar <- update$stdev_matrix_level_addvar
    }
  }
  
  #### Local Level + Explanatory Variables + Slope ####
  if (!is.null(slope_addvar_list) | (!is.null(level_addvar_list) & slope_ind)) {
    
    # Constructing matrix with 0s using former m dimension
    zero_mat <- matrix(0, p, m)
    
    # Check which dependent variables will be excluded, removing doubles as well
    exclude_level <- exclude_level[which(exclude_level >= 1 & exclude_level <= p)]
    exclude_level <- unique(exclude_level)
    
    # The number of dependent variables that should get a local level
    n_level <- p - length(exclude_level)
    
    # Check which dependent variables will not get a slope, removing doubles as well
    exclude_slope <- exclude_slope[which(exclude_slope >= 1 & exclude_slope <= p)]
    exclude_slope <- c(exclude_level, exclude_slope)
    exclude_slope <- unique(exclude_slope)
    
    # The number of local levels that should get a slope
    n_slope <- p - length(exclude_slope)
    
    # Number of coefficients
    k_level <- sum(
      sapply(slope_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[2]}})
    )
    
    # Number of observations
    N <- max(
      sapply(slope_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[1]}})
    )
    
    # Saving former m
    m_old <- m
    
    # Updating m
    m <- m + n_level + n_slope + k_level
    
    # Indices of the coefficients in the state vector
    coeff_list$slope_addvar_state <- (m_old + n_level + n_slope + 1):m
    
    # Function for naming the explanatory variables
    names_fun <- function(i) {
      if (!is.null(slope_addvar_list[[i]]) & 
          is.null(colnames(slope_addvar_list[[i]]))) {
        return(
          rep(
            paste0("Explanatory Variable in level y", i), 
            dim(slope_addvar_list[[i]])[2]
          )
        )
      } else {
        return(colnames(slope_addvar_list[[i]]))
      }
    }
    
    # Label of the state parameters
    state_label <- c(
      state_label, 
      paste0("Local Level y", (1:p)[which(!1:p %in% exclude_level)])
    )
    state_label <- c(
      state_label, 
      paste0("Slope y", (1:p)[which(!1:p %in% exclude_slope)])
    )
    state_label <- c(state_label, do.call(c, lapply(1:p, names_fun)))
    
    # How many parameters in param vector are meant for the level component?
    if (is.null(format_level)) {
      format_level <- diag(1, n_level, n_level)
    }
    # Only entries in lower triangle of matrix count
    format_level[upper.tri(format_level)] <- 0 
    param_num_list$level <- sum(
      format_level != 0 & 
      lower.tri(format_level, diag = TRUE)
    )
    
    # How many parameters in param vector are meant for the slope component?
    if (is.null(format_slope)) {
      format_slope <- diag(1, n_slope, n_slope)
    }
    # Only entries in lower triangle of matrix count
    format_slope[upper.tri(format_slope)] <- 0 
    param_num_list$slope <- sum(
      format_slope != 0 & 
      lower.tri(format_slope, diag = TRUE)
    )
    
    # How many parameters in param vector are meant for the explanatory variables?
    if (is.null(format_level_addvar)) {
      format_level_addvar <- matrix(0, k_level, k_level)
    }
    # Only entries in lower triangle of matrix count
    format_level_addvar[upper.tri(format_level_addvar)] <- 0 
    param_num_list$level_addvar <- sum(
      format_level_addvar != 0 & 
      lower.tri(format_level_addvar, diag = TRUE)
    )
    
    # Last index of parameter that should be used for the component
    index_par2 <- index_par + param_num_list$level + param_num_list$slope + 
      param_num_list$level_addvar - 1
    
    # Indices of parameters that should be used for the component
    param_indices$slope_addvar <- index_par:index_par2
    
    # Keeping track of how many parameters the State Space model needs
    param_num <- param_num + param_num_list$level + param_num_list$slope + 
      param_num_list$level_addvar
    
    # Calling the proper function to obtain system matrices
    update <- SlopeAddVar(p = p,
                          exclude_level = exclude_level,
                          exclude_slope = exclude_slope,
                          slope_addvar_list = slope_addvar_list,
                          fixed_part = TRUE,
                          update_part = update_part,
                          param = param[param_indices$slope_addvar],
                          format_level = format_level,
                          format_slope = format_slope,
                          format_level_addvar = format_level_addvar,
                          decompositions = TRUE
    )
    
    # Updating indices for the first parameter to use for the next component
    index_par <- index_par2 + 1
    
    # Storing system matrices
    Z_list$slope_addvar <- update$Z
    if (Zdim < 3) {
      Z_kal <- cbind(Z_kal, update$Z)
    } else {
      Z_kal <- array(
        apply(
          Z_kal, 3, 
          function(x) {cbind(
            matrix(x, p, dim(Z_kal)[2]), 
            update$Z
          )}
        ), 
        dim = c(p, sum(dim(Z_kal)[2], k_level), N)
      )
    }
    T_list$slope_addvar <- update$Tmat
    Tdim <- 3
    T_kal <- array(
      apply(
        update$Tmat, 3, 
        function(x) BlockMatrix(T_kal, as.matrix(x))
      ), 
      dim = c(sum(dim(T_kal)[1], k_level), sum(dim(T_kal)[2], k_level), N)
    )
    R_list$slope_addvar <- update$R
    R_kal <- BlockMatrix(R_kal, update$R)
    a_list$slope_addvar <- update$a1
    a <- rbind(a, update$a1)
    Pinf_list$slope_addvar <- update$P_inf
    P_inf <- BlockMatrix(P_inf, update$P_inf)
    Pstar_list$slope_addvar <- update$P_star
    P_star <- BlockMatrix(P_star, update$P_star)
    Z_padded_list$level <- cbind(zero_mat, update$Z)
    if (param_num_list$level == 0 | update_part) {
      Q_list$level <- update$Q_level
      Q_kal <- BlockMatrix(Q_kal, update$Q_level)
      L_list$level <- update$loading_matrix_level
      D_list$level <- update$diagonal_matrix_level
      corr_list$level <- update$correlation_matrix_level
      stdev_list$level <- update$stdev_matrix_level
    }
    if (param_num_list$slope == 0 | update_part) {
      Q_list$slope <- update$Q_slope
      Q_kal <- BlockMatrix(Q_kal, update$Q_slope)
      L_list$slope <- update$loading_matrix_slope
      D_list$slope <- update$diagonal_matrix_slope
      corr_list$slope <- update$correlation_matrix_slope
      stdev_list$slope <- update$stdev_matrix_slope
    }
    if (param_num_list$level_addvar == 0 | update_part) {
      Q_list$level_addvar <- update$Q_level_addvar
      Q_kal <- BlockMatrix(Q_kal, update$Q_level_addvar)
      L_list$level_addvar <- update$loading_matrix_level_addvar
      D_list$level_addvar <- update$diagonal_matrix_level_addvar
      corr_list$level_addvar <- update$correlation_matrix_level_addvar
      stdev_list$level_addvar <- update$stdev_matrix_level_addvar
    }
  }
  
  #### Cycle ####
  if (cycle_ind) {
    
    # Initialising lists to store lambda and rho parameters of the cycles
    lambda_list <- list()
    if (sum(damping_factor_ind) > 0) {
      rho_list <- list()
    }
    
    for (i in seq_along(format_cycle_list)) {
      
      # Constructing matrix with 0s using former m dimension
      zero_mat <- matrix(0, p, m)
      
      # Check which dependent variables will be excluded, removing doubles as well
      exclude_cycle_list[[i]] <- exclude_cycle_list[[i]][
        which(exclude_cycle_list[[i]] >= 1 & exclude_cycle_list[[i]] <= p)
      ]
      exclude_cycle_list[[i]] <- unique(exclude_cycle_list[[i]])
      
      # The number of dependent variables that should get a cycle
      n_cycle <- p - length(exclude_cycle_list[[i]])
      
      # Saving former m
      m_old <- m
      
      # Updating m
      m <- m + 2 * n_cycle
      
      # Label of the state parameters
      state_label <- c(
        state_label, 
        rep(
          paste0("Cycle", i, " y", (1:p)[which(!1:p %in% exclude_cycle_list[[i]])]),
          2
        )
      )
      
      # How many parameters in param vector are meant for the cycle component?
      if (is.null(format_cycle_list[[i]])) {
        format_cycle_list[[i]] <- diag(1, n_cycle, n_cycle)
      }
      # Only entries in lower triangle of matrix count
      format_cycle_list[[i]][upper.tri(format_cycle_list[[i]])] <- 0 
      param_num_list[[paste0('Cycle', i)]] <- 1 + damping_factor_ind[i] + 
        sum(format_cycle_list[[i]] != 0 & 
            lower.tri(format_cycle_list[[i]], diag = TRUE)
        )
      
      # Last index of parameter that should be used for the component
      index_par2 <- index_par + param_num_list[[paste0('Cycle', i)]] - 1
      
      # Indices of parameters that should be used for the component
      param_indices[[paste0('Cycle', i)]] <- index_par:index_par2
      
      # Keeping track of how many parameters the State Space model needs
      param_num <- param_num + param_num_list[[paste0('Cycle', i)]]
      
      # Calling the proper function to obtain system matrices
      update <- Cycle(p = p,
                      exclude_cycle = exclude_cycle_list[[i]],
                      damping_factor_ind = damping_factor_ind[i],
                      fixed_part = TRUE,
                      update_part = update_part,
                      param = param[param_indices[[paste0('Cycle', i)]]],
                      format_cycle = format_cycle_list[[i]],
                      decompositions = TRUE
      )
      
      # Updating indices for the first parameter to use for the next component
      index_par <- index_par2 + 1
      
      # Storing system matrices
      Z_list[[paste0('Cycle', i)]] <- update$Z
      if (Zdim < 3) {
        Z_kal <- cbind(Z_kal, update$Z)
      } else {
        Z_kal <- array(
          apply(
            Z_kal, 3, 
            function(x) {cbind(
              matrix(x, p, dim(Z_kal)[2]), 
              update$Z
            )}
          ), 
          dim = c(p, sum(dim(Z_kal)[2], 2 * n_cycle), N)
        )
      }
      R_list[[paste0('Cycle', i)]] <- update$R
      R_kal <- BlockMatrix(R_kal, update$R)
      a_list[[paste0('Cycle', i)]] <- update$a1
      a <- rbind(a, update$a1)
      Pinf_list[[paste0('Cycle', i)]] <- update$P_inf
      P_inf <- BlockMatrix(P_inf, update$P_inf)
      Z_padded_list[[paste0('Cycle', i)]] <- cbind(zero_mat, update$Z)
      if (param_num_list[[paste0('Cycle', i)]] == (1 + damping_factor_ind[i]) | 
          update_part) {
        Q_list[[paste0('Cycle', i)]] <- update$Q_cycle
        Q_list2[[paste0('Cycle', i)]] <- update$Q
        Q_kal <- BlockMatrix(Q_kal, update$Q)
        L_list[[paste0('Cycle', i)]] <- update$loading_matrix
        D_list[[paste0('Cycle', i)]] <- update$diagonal_matrix
        corr_list[[paste0('Cycle', i)]] <- update$correlation_matrix
        stdev_list[[paste0('Cycle', i)]] <- update$stdev_matrix
      }
      if (param_num_list[[paste0('Cycle', i)]] == (1 + damping_factor_ind[i]) | 
          !damping_factor_ind[i] | update_part) {
        Pstar_list[[paste0('Cycle', i)]] <- update$P_star
        P_star <- BlockMatrix(P_star, update$P_star)
      }
      if (update_part) {
        T_list[[paste0('Cycle', i)]] <- update$Tmat
        if (Tdim < 3) {
          T_kal <- BlockMatrix(T_kal, update$Tmat)
        } else {
          T_kal <- array(
            apply(T_kal, 3, 
                  function(x) BlockMatrix(as.matrix(x), update$Tmat)
            ), 
            dim = c(sum(dim(T_kal)[1], 2 * n_cycle), 
                    sum(dim(T_kal)[2], 2 * n_cycle), 
                    N
            )
          )
        }
        lambda_list[[paste0('Cycle', i)]] <- update$lambda
        if (damping_factor_ind[i]) {
          rho_list[[paste0('Cycle', i)]] <- update$rho
        }
      }
    }
  }
  
  #### ARIMA ####
  if (!is.null(arima_list)) {
    
    # Initialising lists to store AR and MA coefficients of the ARIMA components
    ar_list <- list()
    ma_list <- list()
    
    for (i in seq_along(arima_list)) {
      
      # Constructing matrix with 0s using former m dimension
      zero_mat <- matrix(0, p, m)
      
      # Check which dependent variables will be excluded, removing doubles as well
      exclude_arima_list[[i]] <- exclude_arima_list[[i]][
        which(exclude_arima_list[[i]] >= 1 & exclude_arima_list[[i]] <= p)
      ]
      exclude_arima_list[[i]] <- unique(exclude_arima_list[[i]])
      
      # The number of dependent variables that are involved in the ARIMA component
      n_arima <- p - length(exclude_arima_list[[i]])
      
      # The number of state parameters
      state_num <- (max(arima_list[[i]][1], arima_list[[i]][3] + 1) + 
                      arima_list[[i]][2]) * n_arima
      
      # Saving former m
      m_old <- m
      
      # Updating m
      m <- m + state_num
      
      # Label of the state parameters
      state_label <- c(state_label, rep(paste0("ARIMA", i), state_num))
      
      # How many parameters in param vector are meant for the ARIMA component?
      param_num_list[[paste0('ARIMA', i)]] <- 0.5 * n_arima * (n_arima + 1) +
        n_arima^2 * (arima_list[[i]][1] + arima_list[[i]][3])
      
      # Last index of parameter that should be used for the component
      index_par2 <- index_par + param_num_list[[paste0('ARIMA', i)]] - 1
      
      # Indices of parameters that should be used for the component
      param_indices[[paste0('ARIMA', i)]] <- index_par:index_par2
      
      # Keeping track of how many parameters the State Space model needs
      param_num <- param_num + param_num_list[[paste0('ARIMA', i)]]
      
      # Calling the proper function to obtain system matrices
      update <- ARIMA(p = p,
                      arima_spec = arima_list[[i]],
                      exclude_arima = exclude_arima_list[[i]],
                      fixed_part = TRUE,
                      update_part = update_part,
                      param = param[param_indices[[paste0('ARIMA', i)]]],
                      decompositions = TRUE
      )
      
      # Updating indices for the first parameter to use for the next component
      index_par <- index_par2 + 1
      
      # Storing system matrices
      Z_list[[paste0('ARIMA', i)]] <- update$Z
      if (Zdim < 3) {
        Z_kal <- cbind(Z_kal, update$Z)
      } else {
        Z_kal <- array(
          apply(
            Z_kal, 3, 
            function(x) {cbind(
              matrix(x, p, dim(Z_kal)[2]), 
              update$Z
            )}
          ), 
          dim = c(p, sum(dim(Z_kal)[2], state_num), N)
        )
      }
      a_list[[paste0('ARIMA', i)]] <- update$a1
      a <- rbind(a, update$a1)
      Pinf_list[[paste0('ARIMA', i)]] <- update$P_inf
      P_inf <- BlockMatrix(P_inf, update$P_inf)
      Z_padded_list[[paste0('ARIMA', i)]] <- cbind(zero_mat, update$Z)
      if (update_part) {
        if (arima_list[[i]][1] > 0) {
          ar_list[[paste0('ARIMA', i)]] <- update$ar
        }
        if (arima_list[[i]][3] > 0) {
          ma_list[[paste0('ARIMA', i)]] <- update$ma
        }
        R_list[[paste0('ARIMA', i)]] <- update$R
        R_kal <- BlockMatrix(R_kal, update$R)
        T_list[[paste0('ARIMA', i)]] <- update$Tmat
        if (Tdim < 3) {
          T_kal <- BlockMatrix(T_kal, update$Tmat)
        } else {
          T_kal <- array(
            apply(T_kal, 3, 
                  function(x) BlockMatrix(as.matrix(x), update$Tmat)
            ), 
            dim = c(
              sum(dim(T_kal)[1], state_num), 
              sum(dim(T_kal)[2], state_num), 
              N
            )
          )
        }
        Q_list[[paste0('ARIMA', i)]] <- update$Q
        Q_kal <- BlockMatrix(Q_kal, update$Q)
        L_list[[paste0('ARIMA', i)]] <- update$loading_matrix
        D_list[[paste0('ARIMA', i)]] <- update$diagonal_matrix
        corr_list[[paste0('ARIMA', i)]] <- update$correlation_matrix
        stdev_list[[paste0('ARIMA', i)]] <- update$stdev_matrix
        Pstar_list[[paste0('ARIMA', i)]] <- update$P_star
        P_star <- BlockMatrix(P_star, update$P_star)
      }
      if (!update_part & arima_list[[i]][1] == 0 & arima_list[[i]][3] == 0) {
        R_list[[paste0('ARIMA', i)]] <- update$R
        T_list[[paste0('ARIMA', i)]] <- update$Tmat
      }
      if (!update_part & (arima_list[[i]][1] > 0 | arima_list[[i]][3] > 0)) {
        temp_list[[paste0('ARIMA', i)]] <- list(
          T1 = update$T1,
          T2 = update$T2,
          T3 = update$T3,
          R1 = update$R1,
          R2 = update$R2
        )
      }
    }
  }
  
  #### SARIMA ####
  if (!is.null(sarima_list)) {
    
    # Initialising lists to store AR and MA coefficients of the SARIMA components
    sar_list <- list()
    sma_list <- list()
    
    for (i in seq_along(sarima_list)) {
      
      # Constructing matrix with 0s using former m dimension
      zero_mat <- matrix(0, p, m)
      
      # Check which dependent variables will be excluded, removing doubles as well
      exclude_sarima_list[[i]] <- exclude_sarima_list[[i]][
        which(exclude_sarima_list[[i]] >= 1 & exclude_sarima_list[[i]] <= p)
      ]
      exclude_sarima_list[[i]] <- unique(exclude_sarima_list[[i]])
      
      # The number of dependent variables that are involved in the SARIMA component
      n_sarima <- p - length(exclude_sarima_list[[i]])
      
      # The number of state parameters
      state_num <- (max(
        sum(sarima_list[[i]]$s * sarima_list[[i]]$ar), 
        sum(sarima_list[[i]]$s * sarima_list[[i]]$ma) + 1
      ) + sum(sarima_list[[i]]$s * sarima_list[[i]]$i)) * n_sarima
      
      # Saving former m
      m_old <- m
      
      # Updating m
      m <- m + state_num
      
      # Label of the state parameters
      state_label <- c(state_label, rep(paste0("SARIMA", i), state_num))
      
      # How many parameters in param vector are meant for the SARIMA component?
      param_num_list[[paste0('SARIMA', i)]] <- 0.5 * n_sarima * (n_sarima + 1) +
        n_sarima^2 * (sum(sarima_list[[i]]$ar) + sum(sarima_list[[i]]$ma))
      
      # Last index of parameter that should be used for the component
      index_par2 <- index_par + param_num_list[[paste0('SARIMA', i)]] - 1
      
      # Indices of parameters that should be used for the component
      param_indices[[paste0('SARIMA', i)]] <- index_par:index_par2
      
      # Keeping track of how many parameters the State Space model needs
      param_num <- param_num + param_num_list[[paste0('SARIMA', i)]]
      
      # Calling the proper function to obtain system matrices
      update <- SARIMA(p = p,
                       sarima_spec = sarima_list[[i]],
                       exclude_sarima = exclude_sarima_list[[i]],
                       fixed_part = TRUE,
                       update_part = update_part,
                       param = param[param_indices[[paste0('SARIMA', i)]]],
                       decompositions = TRUE
      )
      
      # Updating indices for the first parameter to use for the next component
      index_par <- index_par2 + 1
      
      # Storing system matrices
      Z_list[[paste0('SARIMA', i)]] <- update$Z
      if (Zdim < 3) {
        Z_kal <- cbind(Z_kal, update$Z)
      } else {
        Z_kal <- array(
          apply(
            Z_kal, 3, 
            function(x) {cbind(
              matrix(x, p, dim(Z_kal)[2]), 
              update$Z
            )}
          ), 
          dim = c(p, sum(dim(Z_kal)[2], state_num), N)
        )
      }
      a_list[[paste0('SARIMA', i)]] <- update$a1
      a <- rbind(a, update$a1)
      Pinf_list[[paste0('SARIMA', i)]] <- update$P_inf
      P_inf <- BlockMatrix(P_inf, update$P_inf)
      Z_padded_list[[paste0('SARIMA', i)]] <- cbind(zero_mat, update$Z)
      if (update_part) {
        if (sum(sarima_list[[i]]$ar) > 0) {
          sar_list[[paste0('SARIMA', i)]] <- update$sar
        }
        if (sum(sarima_list[[i]]$ma) > 0) {
          sma_list[[paste0('SARIMA', i)]] <- update$sma
        }
        R_list[[paste0('SARIMA', i)]] <- update$R
        R_kal <- BlockMatrix(R_kal, update$R)
        T_list[[paste0('SARIMA', i)]] <- update$Tmat
        if (Tdim < 3) {
          T_kal <- BlockMatrix(T_kal, update$Tmat)
        } else {
          T_kal <- array(
            apply(T_kal, 3, 
                  function(x) BlockMatrix(as.matrix(x), update$Tmat)
            ), 
            dim = c(
              sum(dim(T_kal)[1], state_num), 
              sum(dim(T_kal)[2], state_num), 
              N
            )
          )
        }
        Q_list[[paste0('SARIMA', i)]] <- update$Q
        Q_kal <- BlockMatrix(Q_kal, update$Q)
        L_list[[paste0('SARIMA', i)]] <- update$loading_matrix
        D_list[[paste0('SARIMA', i)]] <- update$diagonal_matrix
        corr_list[[paste0('SARIMA', i)]] <- update$correlation_matrix
        stdev_list[[paste0('SARIMA', i)]] <- update$stdev_matrix
        Pstar_list[[paste0('SARIMA', i)]] <- update$P_star
        P_star <- BlockMatrix(P_star, update$P_star)
      }
      if (!update_part & 
          sum(sarima_list[[i]]$ar) == 0 & 
          sum(sarima_list[[i]]$ma) == 0) {
        R_list[[paste0('SARIMA', i)]] <- update$R
        T_list[[paste0('SARIMA', i)]] <- update$Tmat
      }
      if (!update_part & 
          (sum(sarima_list[[i]]$ar) > 0 | sum(sarima_list[[i]]$ma) > 0)) {
        temp_list[[paste0('SARIMA', i)]] <- list(
          T1 = update$T1,
          T2 = update$T2,
          T3 = update$T3,
          R1 = update$R1,
          R2 = update$R2
        )
      }
    }
  }
  
  # Input arguments that specify the state space model
  function_call <- list(H_format = H_format,
                        local_level_ind = local_level_ind,
                        slope_ind = slope_ind,
                        BSM_vec = BSM_vec,
                        cycle_ind = cycle_ind, 
                        addvar_list = addvar_list,
                        level_addvar_list = level_addvar_list,
                        slope_addvar_list = slope_addvar_list,
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
                        format_level_addvar = format_level_addvar)
  
  # Constructing the list to return
  result <- list(
    function_call = function_call,
    H = H_list,
    Z = Z_list,
    Tmat = T_list,
    R = R_list,
    Q = Q_list,
    Q2 = Q_list2,
    temp_list = temp_list,
    Q_loading_matrix = L_list,
    Q_diagonal_matrix = D_list,
    Q_correlation_matrix = corr_list,
    Q_stdev_matrix = stdev_list,
    a1 = a_list,
    P_inf = Pinf_list,
    P_star = Pstar_list,
    Z_padded = Z_padded_list,
    param_indices = param_indices,
    state_label = state_label,
    Z_kal = Z_kal,
    T_kal = T_kal,
    R_kal = R_kal,
    Q_kal = Q_kal,
    a_kal = a,
    P_inf_kal = P_inf,
    P_star_kal = P_star,
    residuals_state = residuals_state,
    param_num = param_num,
    param_num_list = param_num_list
  )
  result$lambda <- lambda_list
  result$rho <- rho_list
  if (length(ar_list) > 0) {
    result$AR <- ar_list
  }
  if (length(ma_list) > 0) {
    result$MA <- ma_list
  }
  if (length(sar_list) > 0) {
    result$SAR <- sar_list
  }
  if (length(sma_list) > 0) {
    result$SMA <- sma_list
  }
  result$addvar_state <- coeff_list$addvar_state
  result$level_addvar_state <- coeff_list$level_addvar_state
  result$slope_addvar_state <- coeff_list$slope_addvar_state
  
  # Return the result
  return(result)
}