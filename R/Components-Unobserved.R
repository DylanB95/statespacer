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
      Z <- Z[, -exclude_level, drop = FALSE]
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
      Z <- Z[, -exclude_level, drop = FALSE]
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
      Tmat <- Tmat[-exclude, -exclude, drop = FALSE]
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
#' @param exclude_cycle The dependent variables that should not get a 
#'   Cycle component.
#' @param damping_factor_ind Boolean indicating whether a damping factor 
#'   should be included.
#' @param format_cycle `format` argument for the 
#'   \code{\link{Cholesky}} function.
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
      Z <- Z[, -exclude_cycle, drop = FALSE]
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
      
      # Using Cholesky function to get a valid variance - covariance matrix
      # for the Q matrix
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
#' @param exclude_BSM The dependent variables that should not get a 
#'   BSM component.
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
  
  # Check for erroneous input
  if (s < 3) {
    stop("The period of the BSM component must be greater than or equal to 3.")
  }
  
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
      Ztemp <- Z[, -exclude_BSM, drop = FALSE]
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
    T1 <- do.call(
      BlockMatrix, 
      lapply(lambda, function(x) {diag(cos(x), n_BSM, n_BSM)})
    )
    T2 <- do.call(
      BlockMatrix, 
      lapply(lambda, function(x) {diag(sin(x), n_BSM, n_BSM)})
    )
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
      result$Q <- do.call(
        BlockMatrix, 
        replicate(s - 1, Q_BSM$cov_mat, simplify = FALSE)
      )
    } else {
      result$Q <- do.call(
        BlockMatrix, 
        replicate(s - 1, Q_BSM, simplify = FALSE)
      )
    }
  }
  
  return(result)
}