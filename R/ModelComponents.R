#' Construct the System Matrices of a Local Level Component
#' 
#' Constructs the system matrices of a Local Level component.
#' 
#' @param fixed_part Boolean indicating whether the system matrices should be 
#'   constructed that do not depend on the parameters.
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceFit
#' @inheritParams StateSpaceEval
#' @inheritParams Cholesky
#' 
#' @return 
#' A list containing the system matrices.
#'  
#' @noRd
LocalLevel <- function(p = 1, 
                       exclude_level = NULL, 
                       fixed_part = TRUE, 
                       update_part = TRUE, 
                       param = rep(1, p), 
                       format_level = diag(1, p, p), 
                       decompositions = TRUE) {

  # The number of dependent variables that should get a local level
  n_level <- p - length(exclude_level)
  
  # Variable that is used to check if Q should be a matrix of 0s
  diag_level <- sum(diag(format_level) != 0)
  
  # Initialising the list to return
  result <- list() 
  
  if (fixed_part) {
    
    # Z matrix = a matrix with rows and columns having at most one 1, 
    #            dimension = p x n_level
    Z <- diag(1, p, p)
    if (n_level < p) {
      Z <- matrix(Z[, -exclude_level], p, n_level)
    }
    result$Z <- Z
    
    # T matrix = a diagonal matrix with ones on the diagonal, 
    #            dimension = n_level
    result$Tmat <- diag(1, n_level, n_level)
    
    # R matrix = a diagonal matrix with ones on the diagonal, 
    #            dimension = n_level
    result$R <- diag(1, n_level, n_level)
    
    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, n_level, 1)
    result$P_inf <- diag(1, n_level, n_level)
    result$P_star <- matrix(0, n_level, n_level)
    
    # Check if Q depends on parameters
    if (diag_level == 0) {
      result$Q <- matrix(0, n_level, n_level)
    }
  }
  
  if (update_part & diag_level > 0) {
    
    # Check whether the number of rows of format_level is valid
    if (dim(format_level)[1] != n_level) {
      stop(
        paste0(
          "The number of rows of `format_level` for the local level component ",
          "must be equal to the number of dependent variables minus ",
          "the number of excluded dependent variables."
        )
      )
    }
    
    # Using Cholesky function to get a valid variance - covariance matrix 
    # for the Q matrix
    Q <- Cholesky(
      param = param, 
      format = format_level, 
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

#' Construct the System Matrices of a Local Level + Slope Component
#' 
#' Constructs the system matrices of a Local Level + Slope component.
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
Slope <- function(p = 1,
                  exclude_level = NULL,
                  exclude_slope = NULL,
                  fixed_part = TRUE, 
                  update_part = TRUE,
                  param = rep(1, 2 * p),
                  format_level = diag(1, p, p),
                  format_slope = diag(1, p, p),
                  decompositions = TRUE) {

  # The number of dependent variables that should get a local level
  n_level <- p - length(exclude_level)
  
  # The number of local levels that should get a slope
  n_slope <- p - length(exclude_slope)
  
  # Variable that is used to check if Q_level should be a matrix of 0s
  diag_level <- sum(diag(format_level) != 0)
  
  # Variable that is used to check if Q_slope should be a matrix of 0s
  diag_slope <- sum(diag(format_slope) != 0)
  
  # Initialising the list to return
  result <- list() 
  
  if (fixed_part) {
    
    # Z matrix = matrix of p rows and (n_level + n_slope) columns
    #   First n_level columns containing only one element equal to 1, 0s elsewhere 
    #   Followed by n_slope columns containing 0s
    Z <- cbind(diag(1, p, p), matrix(0, p, n_slope))
    if (n_level < p) {
      Z <- matrix(Z[, -exclude_level], p, n_level + n_slope)
    }
    result$Z <- Z
    
    # T matrix = a (n_level + n_slope) x (n_level + n_slope) matrix
    Tmat <- rbind(
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
      Tmat <- matrix(
        Tmat[-exclude, -exclude], 
        n_level + n_slope, 
        n_level + n_slope
      )
    }
    result$Tmat <- Tmat
    
    # R matrix = a diagonal (n_level + n_slope) x (n_level + n_slope) matrix, 
    #            containing ones on the diagonal
    result$R <- diag(1, (n_level + n_slope), (n_level + n_slope))
    
    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, n_level + n_slope, 1)
    result$P_inf <- diag(1, n_level + n_slope, n_level + n_slope)
    result$P_star <- matrix(0, n_level + n_slope, n_level + n_slope)
    
    # Check if Q_level depends on parameters
    if (diag_level == 0) {
      result$Q_level <- matrix(0, n_level, n_level)
    }
    
    # Check if Q_slope depends on parameters
    if (diag_slope == 0) {
      result$Q_slope <- matrix(0, n_slope, n_slope)
    }
  }
  
  if (update_part & (diag_level + diag_slope) > 0) {
    
    if (diag_level > 0) {
      
      # Check whether the number of rows of format_level is valid
      if (dim(format_level)[1] != n_level) {
        stop(
          paste(
            "The number of rows of `format_level` for the level + slope component",
            "must be equal to the number of dependent variables minus",
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
            "The number of rows of `format_slope` for the level + slope component",
            "must be equal to the number of local levels minus",
            "the number of excluded local levels."
          )
        )
      }
      
      # Using Cholesky function to get a valid variance - covariance matrix 
      # for the Q matrix for the slope
      Q_slope <- Cholesky(
        param = param[(split + 1) : length(param)],
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
    }
  }
  
  return(result)
}

#' Construct the System Matrices of a Cycle Component
#' 
#' Constructs the system matrices of a Cycle component.
#' 
#' @param exclude_cycle The dependent variables that should not get a Cycle component.
#' @param damping_factor_ind Boolean indicating whether a damping factor should be included.
#' @param format_cycle `format` argument for the \code{\link{Cholesky}} function.
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceEval
#' @inheritParams LocalLevel
#' @inheritParams Cholesky
#' 
#' @return 
#' A list containing the system matrices.
#'  
#' @noRd
Cycle <- function(p = 1,
                  exclude_cycle = NULL,
                  damping_factor_ind = TRUE,
                  fixed_part = TRUE, 
                  update_part = TRUE,
                  param = rep(1, p + 1 + damping_factor_ind),
                  format_cycle = diag(1, p, p),
                  decompositions = TRUE) {

  # The number of dependent variables that should get a cycle
  n_cycle <- p - length(exclude_cycle)
  
  # Variable that is used to check if Q should be a matrix of 0s
  diag_cycle <- sum(diag(format_cycle) != 0)
  
  # Initialising the list to return
  result <- list() 
  
  if (fixed_part) {
    
    # Z matrix = matrix of p rows and 2 * n_cycle columns
    #   First n_cycle columns containing only one element equal to 1, 0s elsewhere
    #   Second n_cycle columns containing zeroes
    Z <- cbind(diag(1, p, p), matrix(0, p, n_cycle))
    if (n_cycle < p) {
      Z <- matrix(Z[, -exclude_cycle], p, 2 * n_cycle)
    }
    result$Z <- Z
    
    # R matrix = a diagonal 2 * n_cycle x 2 * n_cycle matrix, 
    #            containing ones on the diagonal
    result$R <- diag(1, 2 * n_cycle, 2 * n_cycle)
    
    # Initial guess for the State vector
    result$a1 <- matrix(0, 2 * n_cycle, 1)
    
    # Check if Q depends on parameters
    if (diag_cycle == 0) {
      result$Q_cycle <- matrix(0, n_cycle, n_cycle)
      result$Q <- matrix(0, 2 * n_cycle, 2 * n_cycle)
    }
    
    # Check when initial state of the cycle is diffuse
    if (diag_cycle == 0 | !damping_factor_ind) {
      result$P_star <- matrix(0, 2 * n_cycle, 2 * n_cycle)
      result$P_inf <- matrix(1, 2 * n_cycle, 2 * n_cycle)
    } else {
      result$P_inf <- matrix(0, 2 * n_cycle, 2 * n_cycle)
    }
  }
  
  if (update_part) {
    
    # Frequency lambda (> 0)
    result$lambda <- exp(param[1])
    if (is.na(result$lambda)) {
      stop("Not enough parameters supplied.")
    }
    
    
    # T matrix = a 2 * n_cycle x 2 * n_cycle matrix
    result$Tmat <- rbind(
      cbind(
        diag(cos(result$lambda), n_cycle, n_cycle), 
        diag(sin(result$lambda), n_cycle, n_cycle)
      ), 
      cbind(
        diag(-sin(result$lambda), n_cycle, n_cycle), 
        diag(cos(result$lambda), n_cycle, n_cycle)
      )
    )
    
    # Damping factor rho (between 0 and 1)
    if (damping_factor_ind) {
      result$rho <- 1 / (1 + exp(param[2]))
      if (is.na(result$rho)) {
        stop("Not enough parameters supplied.")
      }
      result$Tmat <- result$rho * result$Tmat
      paramQ <- param[-(1:2)]
    } else {
      paramQ <- param[-1]
    }
    
    if (diag_cycle > 0) {
      
      # Check whether the number of rows of format_cycle is valid
      if (dim(format_cycle)[1] != n_cycle) {
        stop(
          paste(
            "The number of rows of `format_cycle` for the cycle component",
            "must be equal to the number of dependent variables minus",
            "the number of excluded dependent variables."
          )
        )
      }
      
      # Using Cholesky function to get a valid variance - covariance matrix for the Q matrix
      Q_cycle <- Cholesky(
        param = paramQ, 
        format = format_cycle, 
        decompositions = decompositions
      )
      
      # Check what to return
      if (decompositions) {
        result$Q_cycle <- Q_cycle$cov_mat
        result$loading_matrix <- Q_cycle$loading_matrix
        result$diagonal_matrix <- Q_cycle$diagonal_matrix
        result$correlation_matrix <- Q_cycle$correlation_matrix
        result$stdev_matrix <- Q_cycle$stdev_matrix
        result$Q <- BlockMatrix(Q_cycle$cov_mat, Q_cycle$cov_mat)
      } else {
        result$Q <- BlockMatrix(Q_cycle, Q_cycle)
      }
      
      # Initial uncertainty of the state vector
      if (damping_factor_ind) {
        result$P_star <- result$Q / (1 - result$rho^2)
      }
    }
  }
  
  return(result)
}

#' Construct the System Matrices of a BSM Seasonality Component
#' 
#' Constructs the system matrices of a BSM seasonality component.
#' 
#' @param s The number of periods for the BSM seasonality.
#' @param exclude_BSM The dependent variables that should not get a BSM component.
#' @param format_BSM `format` argument for the \code{\link{Cholesky}} function.
#' @inheritParams GetSysMat
#' @inheritParams LocalLevel
#' @inheritParams StateSpaceEval
#' @inheritParams Cholesky
#'
#' @return 
#' A list containing the system matrices.
#'  
#' @noRd
BSM <- function(p = 1,
                s = 7,
                exclude_BSM = NULL,
                fixed_part = TRUE, 
                update_part = TRUE,
                param = rep(1, p),
                format_BSM = diag(1, p, p),
                decompositions = TRUE) {

  # The number of dependent variables that should get a bsm component
  n_BSM <- p - length(exclude_BSM)
  
  # Variable that is used to check if Q should be a matrix of 0s
  diag_BSM <- sum(diag(format_BSM) != 0)
  
  # Initialising the list to return
  result <- list()
  
  # Number of lambdas
  if ((s %% 2) == 0) {
    
    # If s is even then s/2 - 1 lambdas
    lambda_num <- s / 2 - 1
    
  } else {
    
    # If s is odd then (s-1)/2 lambdas
    lambda_num <- floor(s / 2)
  }
  
  if (fixed_part) {  
    
    # Z matrix = matrix of p rows and (s - 1) * n_BSM columns
    #   First n_BSM * lambda_num columns: lambda_num p x n_BSM matrices 
    #                                     containing only one 1 in each 
    #                                     column, 0s elsewhere
    #   Second n_BSM * lambda_num columns: lambda_num p x n_BSM matrices 
    #                                      containing 0s
    #   (if s is even, then last n_BSM columns: p x n_BSM matrix containing 
    #                                           only one 1 in each column, 
    #                                           0s elsewhere)
    Ztemp <- diag(1, p, p)
    if (n_BSM < p) {
      Ztemp <- matrix(Z[, -exclude_BSM], p, n_BSM)
    }
    Z <- cbind(
      matrix(
        rep(Ztemp, lambda_num), 
        p, lambda_num * n_BSM
      ),
      matrix(0, p, lambda_num * n_BSM)
    )
    if ((s %% 2) == 0) {
      Z <- cbind(Z, Ztemp)
    }
    result$Z <- Z
    
    # lambdas used for the T matrix
    lambda <- 2 * pi * (1:lambda_num) / s
    
    # T matrix = a n_BSM * 2 * lambda_num (+ n_BSM if s is even) 
    #            x n_BSM * 2 * lambda_num (+ n_BSM if s is even) matrix
    T1 <- do.call(BlockMatrix, lapply(lambda, function(x) {diag(cos(x), n_BSM, n_BSM)}))
    T2 <- do.call(BlockMatrix, lapply(lambda, function(x) {diag(sin(x), n_BSM, n_BSM)}))
    Tmat <- rbind(cbind(T1, T2), cbind(-T2, T1))
    if ((s %% 2) == 0) {
      Tmat <- BlockMatrix(Tmat, diag(-1, n_BSM, n_BSM))
    }
    result$Tmat <- Tmat
    
    # R matrix = a diagonal (s - 1) * n_BSM x (s - 1) * n_BSM matrix, 
    #            containing ones on the diagonal
    result$R <- diag(1, (s - 1) * n_BSM, (s - 1) * n_BSM)
    
    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, (s - 1) * n_BSM, 1)
    result$P_inf <- diag(1, (s - 1) * n_BSM, (s - 1) * n_BSM)
    result$P_star <- matrix(0, (s - 1) * n_BSM, (s - 1) * n_BSM)
    
    # Check if Q depends on parameters
    if (diag_BSM == 0) {
      result$Q_BSM <- matrix(0, n_BSM, n_BSM)
      result$Q <- matrix(0, (s - 1) * n_BSM, (s - 1) * n_BSM)
    }
  }
  
  if (update_part & diag_BSM > 0) {
    
    # Check whether the number of rows of format_BSM is valid
    if (dim(format_BSM)[1] != n_BSM) {
      stop(
        paste0(
          "The number of rows of `format_BSM` for the BSM", s, " component ",
          "must be equal to the number of dependent variables minus ",
          "the number of excluded dependent variables."
        )
      )
    }
    
    # Using Cholesky function to get a valid variance - covariance matrix 
    # for the Q matrix
    Q_BSM <- Cholesky(
      param = param, 
      format = format_BSM, 
      decompositions = decompositions
    )
    
    # Check what to return
    if (decompositions) {
      result$Q_BSM <- Q_BSM$cov_mat
      result$loading_matrix <- Q_BSM$loading_matrix
      result$diagonal_matrix <- Q_BSM$diagonal_matrix
      result$correlation_matrix <- Q_BSM$correlation_matrix
      result$stdev_matrix <- Q_BSM$stdev_matrix
      result$Q <- do.call(BlockMatrix, replicate(s - 1, Q_BSM$cov_mat, simplify = FALSE))
    } else {
      result$Q <- do.call(BlockMatrix, replicate(s - 1, Q_BSM, simplify = FALSE))
    }
  }
  
  return(result)
}

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
        function(i) {matrix(Ztemp2[seq(i, p * N, N),], p, k)}, 
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

#' Construct the System Matrices of a Explanatory Variables in the Level Component
#' 
#' Constructs the system matrices of a explanatory variables in the level component.
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
      Z <- matrix(Z[, -exclude_level], p, n_level + k)
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
          matrix(Ttemp2[seq(i, n_level * N, N),], n_level, k)
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
            "explanatory variables component must be equal to the",
            "number of dependent variables minus the number of",
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
      
      # Check whether the number of rows of format_slope is valid
      if (dim(format_level_addvar)[1] != k) {
        stop(
          paste(
            "The number of rows of `format_level_addvar` for the level +",
            "explanatory variables component must be equal to the number",
            "of explanatory variables."
          )
        )
      }
      
      # Using Cholesky function to get a valid variance - covariance matrix 
      # for the Q matrix for the slope
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

#' Construct the System Matrices of a Explanatory Variables in the Level + Slope Component
#' 
#' Constructs the system matrices of a explanatory variables in the level + slope component.
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
                        slope_addvar_list,
                        fixed_part = TRUE,
                        update_part = TRUE,
                        param = rep(1, 2 * p),
                        format_level = diag(1, p, p),
                        format_slope = diag(1, p, p),
                        format_level_addvar = diag(0, 1, 1),
                        decompositions = TRUE) {

  # Check if the slope_addvar list has p elements
  if (length(slope_addvar_list) != p) {
    stop(
      paste(
        "The number of elements in `slope_addvar_list` must be equal",
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
    k <- sum(sapply(slope_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[2]}}))
    
    # Number of observations
    N <- max(sapply(slope_addvar_list, function(X) {if (is.null(X)) {0} else {dim(X)[1]}}))
    
    # Z matrix = matrix of p rows and n_level + n_slope + k columns
    #   First n_level columns containing one 1, 0s elsewhere 
    #   Second n_slope columns containing zeroes
    #   Last k columns containing zeroes
    Z <- cbind(diag(1, p, p), matrix(0, p, n_slope), matrix(0, p, k))
    if (n_level < p) {
      Z <- matrix(Z[, -exclude_level], p, n_level + n_slope + k)
    }
    result$Z <- Z
    
    # Tmat = a (n_level + n_slope + k) x (n_level + n_slope + k) x N array
    Ttemp <- do.call(BlockMatrix, slope_addvar_list)
    index <- 1:N
    Ttemp2 <- NULL
    null_mat <- matrix(0, N, k)
    for (i in seq_along(slope_addvar_list)) {
      if (is.null(slope_addvar_list[[i]])) {
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
      Ttemp3 <- matrix(
        Ttemp3[-exclude, -exclude], 
        n_level + n_slope, 
        n_level + n_slope
      )
    }
    result$Tmat <- sapply(
      seq(N), 
      function(i) {rbind(
        cbind(
          Ttemp3, 
          rbind(
            matrix(
              Ttemp2[seq(i, n_level * N, N),], 
              n_level, k
            ),
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
            "+ slope + explanatory variables component must be",
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
            "+ slope + explanatory variables component must be",
            "equal to the number of local levels minus",
            "the number of excluded local levels."
          )
        )
      }
      
      # Which parameters should be used for the Q matrix corresponding to the slope
      split2 <- split + sum(format_slope != 0 & lower.tri(format_slope, diag = TRUE))
      
      # Using Cholesky function to get a valid variance - covariance matrix 
      # for the Q matrix for the level
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
      
      # Check whether the number of rows of format_slope is valid
      if (dim(format_level_addvar)[1] != k) {
        stop(
          paste(
            "The number of rows of `format_level_addvar` for the level",
            "+ slope + explanatory variables component must be",
            "equal to the number of explanatory variables."
          )
        )
      }
      
      # Using Cholesky function to get a valid variance - covariance matrix 
      # for the Q matrix for the slope
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

#' Construct the System Matrices of an ARIMA Component
#' 
#' Constructs the system matrices of an ARIMA component.
#' 
#' @param arima_spec specification of the ARIMA part, should be a vector of 
#'   length 3 with the following format: c(AR, I, MA).
#' @param exclude_arima The dependent variables that should not get an ARIMA component.
#' @param R_stationary The R system matrix that will be used for calculating 
#'   the initialisation matrix P_star
#' @inheritParams LocalLevel
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceEval
#' @inheritParams Cholesky
#' 
#' @return 
#' A list containing the system matrices.
#'  
#' @noRd
ARIMA <- function(p = 1, 
                  arima_spec = c(1, 0, 0), 
                  exclude_arima = NULL, 
                  fixed_part = TRUE, 
                  update_part = TRUE, 
                  param = rep(1, p^2 + 0.5 * p * (p+1)), 
                  decompositions = TRUE,
                  R_stationary = NULL) {
  
  # Check for erroneous input 
  if (length(arima_spec) != 3) {
    stop("`arima_spec` must be a vector of length 3.")
  }
  if (arima_spec[2] < 0) {
    stop("Number of differencing in arima_spec must be greater than or equal to 0.")
  }
  
  # The number of dependent variables that should get an ARIMA component
  n_arima <- p - length(exclude_arima)
  
  # Initialising the list to return
  result <- list() 
  
  if (fixed_part) {
    
    # Z matrix
    Z <- diag(1, p, p)
    if (n_arima < p) {
      Z <- matrix(Z[, -exclude_arima], p, n_arima)
    }
    if (arima_spec[2] > 0) {
      Z <- do.call(cbind, replicate(1 + arima_spec[2], Z, simplify = FALSE))
    }
    if (arima_spec[1] > 1) {
      Z <- cbind(Z, matrix(0, p, n_arima * (arima_spec[1] - 1)))
    }
    Z <- cbind(Z, matrix(0, p, n_arima))
    if (arima_spec[3] > 1) {
      Z <- cbind(Z, matrix(0, p, n_arima * (arima_spec[3] - 1)))
    }
    result$Z <- Z
    
    # R matrix
    R <- diag(1, n_arima, n_arima)
    if (arima_spec[1] > 1) {
      R <- rbind(R, matrix(0, (arima_spec[1] - 1) * n_arima, n_arima))
    }
    R <- rbind(R, diag(1, n_arima, n_arima))
    if (arima_spec[3] > 1) {
      R <- rbind(R, matrix(0, (arima_spec[3] - 1) * n_arima, n_arima))
    }
    R_stationary <- R
    result$R_stationary <- R_stationary
    if (arima_spec[2] > 0) {
      R <- rbind(matrix(0, arima_spec[2] * n_arima, n_arima), R)
    }
    result$R <- R
    
    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, dim(R)[1], 1)
    result$P_inf <- BlockMatrix(
      diag(1, arima_spec[2] * n_arima, arima_spec[2] * n_arima),
      matrix(0, dim(R)[1] - arima_spec[2] * n_arima, dim(R)[1] - arima_spec[2] * n_arima)
    )
    
    # T matrix is fixed if no coefficients are needed
    if (arima_spec[1] == 0 & arima_spec[3] == 0) {
      T1 <- matrix(0, 2 * n_arima, 2 * n_arima)
      if (arima_spec[2] > 0) {
        T1 <- cbind(
          matrix(0, 2 * n_arima, arima_spec[2] * n_arima),
          T1
        )
        T2 <- cbind(
          diag(1, n_arima, n_arima),
          matrix(0, n_arima, n_arima)
        )
        T3 <- matrix(0, 0, 2 * n_arima + arima_spec[2] * n_arima)
        for (i in arima_spec[2]:1) {
          T3 <- rbind(
            T3,
            cbind(
              matrix(0, n_arima, (arima_spec[2] - i) * n_arima),
              do.call(cbind, replicate(i, diag(1, n_arima, n_arima), simplify = FALSE)),
              T2
            )
          )
        }
        T1 <- rbind(T3, T1)
      }
      result$Tmat <- T1
    }
  }
  
  if (update_part) {
    
    # Check for number of parameters
    param <- param[which(!is.na(param))]
    needed <- 0.5 * n_arima * (n_arima + 1) + 
              (arima_spec[1] + arima_spec[3]) * n_arima^2
    if (length(param) < needed) {
      stop("Not enough parameters supplied.")
    }
    if (length(param) > needed) {
      stop("Too many parameters supplied.")
    }
    
    # Parameters for the variance - covariance matrix
    Q_indices <- 1:(0.5 * n_arima * (n_arima + 1))
    Qparam <- param[Q_indices]
    
    # Using Cholesky function to get a valid variance - covariance matrix 
    # for the Q matrix
    Q <- Cholesky(
      param = Qparam,
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
    
    # T matrix
    if (arima_spec[1] > 0 | arima_spec[3] > 0) {
      if (n_arima == 1) {
        A <- param[-Q_indices]
      } else {
        A <- array(
          param[-Q_indices], 
          dim = c(n_arima, n_arima, arima_spec[1] + arima_spec[3])
        )
      }
      coeff <- CoeffARMA(
        A = A, variance = result$Q, 
        ar = arima_spec[1], ma = arima_spec[3]
      )
      coeff_row <- matrix(0, n_arima, 0)
      T1 <- matrix(0, 0, 0)
      if (arima_spec[1] > 0) {
        result$ar <- coeff$ar
        coeff_row <- cbind(coeff_row, matrix(coeff$ar, n_arima, n_arima * arima_spec[1]))
        T1 <- BlockMatrix(
          do.call(
            BlockMatrix, 
            replicate(
              arima_spec[1] - 1, 
              diag(1, n_arima, n_arima), 
              simplify = FALSE
            )
          ),
          matrix(0, n_arima, n_arima)
        )
      } else {
        coeff_row <- cbind(coeff_row, matrix(0, n_arima, n_arima))
        T1 <- matrix(0, n_arima, n_arima)
      }
      if (arima_spec[3] > 0) {
        result$ma <- coeff$ma
        coeff_row <- cbind(coeff_row, matrix(coeff$ma, n_arima, n_arima * arima_spec[3]))
        T1 <- BlockMatrix(
          T1,
          do.call(
            BlockMatrix, 
            replicate(
              arima_spec[3] - 1, 
              diag(1, n_arima, n_arima), 
              simplify = FALSE
            )
          )
        )
        T1 <- cbind(T1, matrix(0, dim(T1)[1], n_arima))
      } else {
        coeff_row <- cbind(coeff_row, matrix(0, n_arima, n_arima))
        T1 <- cbind(T1, matrix(0, dim(T1)[1], n_arima))
      }
      T1 <- rbind(coeff_row, T1)
      T_stationary <- T1
      if (arima_spec[2] > 0) {
        T1 <- cbind(
          matrix(0, dim(T1)[1], arima_spec[2] * n_arima),
          T1
        )
        T2 <- cbind(
          diag(1, n_arima, n_arima),
          matrix(0, n_arima, dim(T_stationary)[2] - n_arima)
        )
        T3 <- matrix(0, 0, dim(T1)[2])
        for (i in arima_spec[2]:1) {
          T3 <- rbind(
            T3,
            cbind(
              matrix(0, n_arima, (arima_spec[2] - i) * n_arima),
              do.call(cbind, replicate(i, diag(1, n_arima, n_arima), simplify = FALSE)),
              T2
            )
          )
        }
        T1 <- rbind(T3, T1)
      }
      result$Tmat <- T1
    }
    
    # Initial uncertainty for the stationary part
    if (arima_spec[1] == 0 & arima_spec[3] == 0) {
      result$P_star <- BlockMatrix(
        matrix(0, arima_spec[2] * n_arima, arima_spec[2] * n_arima),
        result$Q,
        result$Q
      )
    } else {
      T_kronecker <- kronecker(T_stationary, T_stationary)
      Tinv <- solve(diag(1, dim(T_kronecker), dim(T_kronecker)) - T_kronecker)
      vecRQR <- matrix(R_stationary %*% result$Q %*% t(R_stationary))
      vecPstar <- Tinv %*% vecRQR
      result$P_star <- BlockMatrix(
        matrix(0, arima_spec[2] * n_arima, arima_spec[2] * n_arima),
        matrix(vecPstar, dim(T_stationary)[1], dim(T_stationary)[2])
      )
    }
  }
  
  return(result)
}