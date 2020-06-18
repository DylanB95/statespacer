#' Construct the System Matrices of an ARIMA Component
#'
#' Constructs the system matrices of an ARIMA component.
#'
#' @param arima_spec Specification of the ARIMA part, should be a vector of
#'   length 3 with the following format: `c(AR, I, MA)`.
#' @param exclude_arima The dependent variables that should not get an ARIMA component.
#' @param T1 Fixed part of Tmat, only used when `fixed_part = FALSE`.
#' @param T2 Fixed part of Tmat, only used when `fixed_part = FALSE`.
#' @param T3 Fixed part of Tmat, only used when `fixed_part = FALSE`.
#' @param R1 Fixed part of R, only used when `fixed_part = FALSE`.
#' @param R2 Fixed part of R, only used when `fixed_part = FALSE`.
#' @inheritParams LocalLevel
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceEval
#' @inheritParams Cholesky
#' @inheritParams Cycle
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
                  transform_only = FALSE,
                  param = rep(1, p^2 + 0.5 * p * (p + 1)),
                  decompositions = TRUE,
                  T1 = NULL,
                  T2 = NULL,
                  T3 = NULL,
                  R1 = NULL,
                  R2 = NULL) {

  # Check for erroneous input
  if (length(arima_spec) != 3) {
    stop(
      "The ARIMA specification must be a vector of length 3.",
      call. = FALSE
    )
  }
  if (arima_spec[[2]] < 0) {
    stop(
      paste(
        "Number of differencing of the ARIMA component must",
        "be greater than or equal to 0."
      ),
      call. = FALSE
    )
  }

  # The number of dependent variables that should get an ARIMA component
  n_arima <- p - length(exclude_arima)

  # Initialising the list to return
  result <- list()

  # Number of stationary elements in the state vector per included dependent
  r <- max(arima_spec[[1]], arima_spec[[3]] + 1)

  if (fixed_part) {

    # Z matrix
    Z <- diag(1, p, p)
    if (n_arima < p) {
      Z <- Z[, -exclude_arima, drop = FALSE]
    }
    if (arima_spec[[2]] > 0) {
      Z <- do.call(cbind, replicate(1 + arima_spec[[2]], Z, simplify = FALSE))
    }
    Z <- cbind(Z, matrix(0, p, n_arima * (r - 1)))
    result$Z <- Z

    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, (r + arima_spec[[2]]) * n_arima, 1)
    result$P_inf <- BlockMatrix(
      diag(1, arima_spec[[2]] * n_arima, arima_spec[[2]] * n_arima),
      matrix(0, r * n_arima, r * n_arima)
    )

    # T and R matrices are fixed if no coefficients are needed
    if (arima_spec[[1]] == 0 & arima_spec[[3]] == 0) {
      T1 <- matrix(0, n_arima, n_arima)
      if (arima_spec[[2]] > 0) {
        T1 <- cbind(
          matrix(0, n_arima, arima_spec[[2]] * n_arima),
          T1
        )
        T2 <- diag(1, n_arima, n_arima)
        T3 <- matrix(0, 0, n_arima + arima_spec[[2]] * n_arima)
        for (i in arima_spec[[2]]:1) {
          T3 <- rbind(
            T3,
            cbind(
              matrix(0, n_arima, (arima_spec[[2]] - i) * n_arima),
              do.call(
                cbind,
                replicate(i, diag(1, n_arima, n_arima), simplify = FALSE)
              ),
              T2
            )
          )
        }
        T1 <- rbind(T3, T1)
      }
      result$Tmat <- T1
      result$R <- rbind(
        matrix(0, arima_spec[[2]] * n_arima, n_arima),
        diag(1, n_arima, n_arima)
      )
    } else {

      # Fixed part of Tmat, T1
      T1 <- diag(1, (r - 1) * n_arima, (r - 1) * n_arima)
      T1 <- rbind(T1, matrix(0, n_arima, (r - 1) * n_arima))
      T1 <- cbind(matrix(0, r * n_arima, n_arima), T1)
      result$T1 <- T1

      # Fixed part of R, R1
      R1 <- diag(1, n_arima, n_arima)
      result$R1 <- R1

      # Non-stationary part of Tmat and R
      if (arima_spec[[2]] > 0) {
        T2 <- matrix(0, dim(T1)[[1]], arima_spec[[2]] * n_arima)
        Ttemp <- cbind(
          diag(1, n_arima, n_arima),
          matrix(0, n_arima, dim(T1)[[2]] - n_arima)
        )
        T3 <- matrix(0, 0, dim(T1)[[2]] + dim(T2)[[2]])
        for (i in arima_spec[[2]]:1) {
          T3 <- rbind(
            T3,
            cbind(
              matrix(0, n_arima, (arima_spec[[2]] - i) * n_arima),
              do.call(
                cbind,
                replicate(i, diag(1, n_arima, n_arima), simplify = FALSE)
              ),
              Ttemp
            )
          )
        }
        result$T2 <- T2
        result$T3 <- T3
        R2 <- matrix(0, arima_spec[[2]] * n_arima, n_arima)
        result$R2 <- R2
      }
    }
  }

  if (update_part) {

    # Check for number of parameters
    param <- param[!is.na(param)]
    needed <- 0.5 * n_arima * (n_arima + 1) +
      (arima_spec[[1]] + arima_spec[[3]]) * n_arima^2
    if (length(param) < needed) {
      stop("Not enough parameters supplied.", call. = FALSE)
    }
    if (length(param) > needed) {
      stop("Too many parameters supplied.", call. = FALSE)
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

    # T matrix and coefficients in R matrix
    if (arima_spec[[1]] > 0 | arima_spec[[3]] > 0) {

      # Coefficients
      if (n_arima == 1) {
        A <- param[-Q_indices]
      } else {
        A <- array(
          param[-Q_indices],
          dim = c(n_arima, n_arima, arima_spec[[1]] + arima_spec[[3]])
        )
      }
      coeff <- CoeffARMA(
        A = A, variance = result$Q,
        ar = arima_spec[[1]], ma = arima_spec[[3]]
      )

      # T matrix coefficients
      if (arima_spec[[1]] > 0) {
        if (decompositions) {
          result$ar <- coeff$ar
        }
        if (!transform_only) {
          if (n_arima > 1) {
            coeff$ar <- aperm(coeff$ar, c(2, 1, 3))
          }
          T1[1:(n_arima * arima_spec[[1]]), 1:n_arima] <- t(
            matrix(
              coeff$ar,
              n_arima,
              n_arima * arima_spec[[1]]
            )
          )
        }
      }
      if (!transform_only) {
        T_stationary <- T1
      }

      # R matrix coefficients
      if (arima_spec[[3]] > 0) {
        if (decompositions) {
          result$ma <- coeff$ma
        }
        if (!transform_only) {
          if (n_arima > 1) {
            coeff$ma <- aperm(coeff$ma, c(2, 1, 3))
          }
          R1 <- rbind(
            R1,
            t(
              matrix(
                coeff$ma,
                n_arima,
                n_arima * arima_spec[[3]]
              )
            )
          )
        }
      }
      if (!transform_only) {
        R1 <- rbind(
          R1,
          matrix(0, (r - 1 - arima_spec[[3]]) * n_arima, n_arima)
        )
        R_stationary <- R1

        # Non-stationary part
        if (arima_spec[[2]] > 0) {
          T1 <- cbind(T2, T1)
          T1 <- rbind(T3, T1)
          R1 <- rbind(R2, R1)
        }
        result$Tmat <- T1
        result$R <- R1
      }
    }

    if (!transform_only) {

      # Initial uncertainty for the stationary part
      if (arima_spec[[1]] == 0 & arima_spec[[3]] == 0) {
        result$P_star <- BlockMatrix(
          matrix(0, arima_spec[[2]] * n_arima, arima_spec[[2]] * n_arima),
          result$Q
        )
      } else {
        T_kronecker <- kronecker(T_stationary, T_stationary)
        Tinv <- solve(diag(1, dim(T_kronecker)[[1]], dim(T_kronecker)[[2]]) - T_kronecker)
        vecRQR <- matrix(tcrossprod(R_stationary %*% result$Q, R_stationary))
        vecPstar <- Tinv %*% vecRQR
        result$P_star <- BlockMatrix(
          matrix(0, arima_spec[[2]] * n_arima, arima_spec[[2]] * n_arima),
          matrix(vecPstar, dim(T_stationary)[[1]], dim(T_stationary)[[2]])
        )
      }
    }
  }

  return(result)
}

#' Construct the System Matrices of a SARIMA Component
#'
#' Constructs the system matrices of a SARIMA component.
#'
#' @param sarima_spec Specification of the SARIMA part, should be a list of
#'   4 named vectors. Vectors should be named: "s", "ar", "i", "ma".
#' @param exclude_sarima The dependent variables that should not get a SARIMA component.
#' @param T1 Fixed part of Tmat, only used when `fixed_part = FALSE`.
#' @param T2 Fixed part of Tmat, only used when `fixed_part = FALSE`.
#' @param T3 Fixed part of Tmat, only used when `fixed_part = FALSE`.
#' @param R1 Fixed part of R, only used when `fixed_part = FALSE`.
#' @param R2 Fixed part of R, only used when `fixed_part = FALSE`.
#' @inheritParams LocalLevel
#' @inheritParams GetSysMat
#' @inheritParams StateSpaceEval
#' @inheritParams Cholesky
#' @inheritParams Cycle
#'
#' @return
#' A list containing the system matrices.
#'
#' @noRd
SARIMA <- function(p = 1,
                   sarima_spec = list(
                     s = c(1, 12),
                     ar = c(1, 1),
                     i = c(0, 0),
                     ma = c(0, 0)
                   ),
                   exclude_sarima = NULL,
                   fixed_part = TRUE,
                   update_part = TRUE,
                   transform_only = FALSE,
                   param = rep(1, 2 * p^2 + 0.5 * p * (p + 1)),
                   decompositions = TRUE,
                   T1 = NULL,
                   T2 = NULL,
                   T3 = NULL,
                   R1 = NULL,
                   R2 = NULL) {

  # Check for erroneous input
  if (length(sarima_spec) != 4) {
    stop(
      "The SARIMA specification must be a list containing 4 vectors.",
      call. = FALSE
    )
  }
  if (length(sarima_spec$s) != length(sarima_spec$ar) |
    length(sarima_spec$s) != length(sarima_spec$i) |
    length(sarima_spec$s) != length(sarima_spec$ma)) {
    stop(
      "The vectors in the SARIMA specification must be of equal length.",
      call. = FALSE
    )
  }
  if (min(sarima_spec$i) < 0) {
    stop(
      paste(
        "Number of differencing in the SARIMA specification must",
        "be greater than or equal to 0."
      ),
      call. = FALSE
    )
  }

  # The number of dependent variables that should get a SARIMA component
  n_sarima <- p - length(exclude_sarima)

  # Initialising the list to return
  result <- list()

  # Number of stationary elements in the state vector per included dependent
  r <- max(
    sum(sarima_spec$s * sarima_spec$ar),
    sum(sarima_spec$s * sarima_spec$ma) + 1
  )

  # Auxiliary seasonality vector
  s_aux <- rep(sarima_spec$s, sarima_spec$i)

  if (fixed_part) {

    # Z matrix
    Z <- diag(1, p, p)
    if (n_sarima < p) {
      Z <- Z[, -exclude_sarima, drop = FALSE]
    }
    Z_list <- lapply(s_aux, function(x) cbind(matrix(0, p, n_sarima * (x - 1)), Z))
    Z_list <- c(Z_list, list(Z, matrix(0, p, (r - 1) * n_sarima)))
    Z <- do.call(cbind, Z_list)
    result$Z <- Z

    # Initialisation of the State vector and corresponding uncertainty
    result$a1 <- matrix(0, (r + sum(s_aux)) * n_sarima, 1)
    result$P_inf <- BlockMatrix(
      diag(1, sum(s_aux) * n_sarima, sum(s_aux) * n_sarima),
      matrix(0, r * n_sarima, r * n_sarima)
    )

    # T and R matrices are fixed if no coefficients are needed
    if (sum(sarima_spec$ar) == 0 & sum(sarima_spec$ma) == 0) {
      T1 <- matrix(0, n_sarima, n_sarima)
      if (sum(sarima_spec$i) > 0) {
        T1 <- cbind(
          matrix(0, n_sarima, sum(s_aux) * n_sarima),
          T1
        )
        T2 <- matrix(0, 0, n_sarima + sum(s_aux) * n_sarima)
        T_list <- lapply(
          s_aux,
          function(x) cbind(matrix(0, n_sarima, n_sarima * (x - 1)), diag(n_sarima))
        )
        T_list <- c(T_list, list(diag(n_sarima)))
        aux <- 1:length(s_aux)
        for (i in seq_along(s_aux)) {
          T2 <- rbind(
            T2,
            do.call(cbind, T_list),
            cbind(
              matrix(
                0,
                (s_aux[[i]] - 1) * n_sarima,
                sum(s_aux[aux < i]) * n_sarima
              ),
              diag((s_aux[[i]] - 1) * n_sarima),
              matrix(
                0,
                (s_aux[[i]] - 1) * n_sarima,
                (1 + sum(s_aux[aux > i]) + r) * n_sarima
              )
            )
          )
          T_list[[i]] <- matrix(0, dim(T_list[[i]])[[1]], dim(T_list[[i]])[[2]])
        }
        T1 <- rbind(T2, T1)
      }
      result$Tmat <- T1
      result$R <- rbind(
        matrix(0, sum(s_aux) * n_sarima, n_sarima),
        diag(1, n_sarima, n_sarima)
      )
    } else {

      # Fixed part of Tmat, T1
      T1 <- diag(1, (r - 1) * n_sarima, (r - 1) * n_sarima)
      T1 <- rbind(T1, matrix(0, n_sarima, (r - 1) * n_sarima))
      T1 <- cbind(matrix(0, r * n_sarima, n_sarima), T1)
      result$T1 <- T1

      # Fixed part of R, R1
      R1 <- diag(1, n_sarima, n_sarima)
      result$R1 <- R1

      # Non-stationary part of Tmat and R
      if (sum(sarima_spec$i) > 0) {
        T2 <- matrix(0, dim(T1)[[1]], sum(s_aux) * n_sarima)
        T3 <- matrix(0, 0, dim(T1)[[2]] + dim(T2)[[2]])
        T_list <- lapply(
          s_aux,
          function(x) cbind(matrix(0, n_sarima, n_sarima * (x - 1)), diag(n_sarima))
        )
        T_list <- c(
          T_list,
          list(diag(n_sarima), matrix(0, n_sarima, (r - 1) * n_sarima))
        )
        aux <- 1:length(s_aux)
        for (i in seq_along(s_aux)) {
          T3 <- rbind(
            T3,
            do.call(cbind, T_list),
            cbind(
              matrix(
                0,
                (s_aux[[i]] - 1) * n_sarima,
                sum(s_aux[aux < i]) * n_sarima
              ),
              diag((s_aux[[i]] - 1) * n_sarima),
              matrix(
                0,
                (s_aux[[i]] - 1) * n_sarima,
                (1 + sum(s_aux[aux > i]) + r) * n_sarima
              )
            )
          )
          T_list[[i]] <- matrix(0, dim(T_list[[i]])[[1]], dim(T_list[[i]])[[2]])
        }
        result$T2 <- T2
        result$T3 <- T3
        R2 <- matrix(0, sum(s_aux) * n_sarima, n_sarima)
        result$R2 <- R2
      }
    }
  }

  if (update_part) {

    # Check for number of parameters
    param <- param[!is.na(param)]
    needed <- 0.5 * n_sarima * (n_sarima + 1) +
      (sum(sarima_spec$ar) + sum(sarima_spec$ma)) * n_sarima^2
    if (length(param) < needed) {
      stop("Not enough parameters supplied.", call. = FALSE)
    }
    if (length(param) > needed) {
      stop("Too many parameters supplied.", call. = FALSE)
    }

    # Parameters for the variance - covariance matrix
    Q_indices <- 1:(0.5 * n_sarima * (n_sarima + 1))
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

    # Drop parameters that are already used
    param <- param[-Q_indices]

    # T matrix and coefficients in R matrix
    if (sum(sarima_spec$ar) > 0 | sum(sarima_spec$ma) > 0) {
      if (!transform_only) {

        # Initialise coefficient vector/array for the AR part
        if (n_sarima == 1) {
          AR_coeff <- rep(0, sum(sarima_spec$s * sarima_spec$ar))
        } else {
          AR_coeff <- array(
            0,
            dim = c(n_sarima, n_sarima, sum(sarima_spec$s * sarima_spec$ar))
          )
        }
        ar_polynomial <- c()

        # Initialise coefficient vector/array for the MA part
        if (n_sarima == 1) {
          MA_coeff <- rep(0, sum(sarima_spec$s * sarima_spec$ma))
        } else {
          MA_coeff <- array(
            0,
            dim = c(n_sarima, n_sarima, sum(sarima_spec$s * sarima_spec$ma))
          )
        }
        ma_polynomial <- c()
      }

      # Initialise lists for AR and MA coefficients
      if (decompositions) {
        sar <- list()
        sma <- list()
      }

      for (i in seq_along(sarima_spec$s)) {
        if ((sarima_spec$ar[[i]] + sarima_spec$ma[[i]]) == 0) {
          next
        }

        # Coefficients
        ARMA_indices <- 1:((sarima_spec$ar[[i]] + sarima_spec$ma[[i]]) * n_sarima^2)
        if (n_sarima == 1) {
          A <- param[ARMA_indices]
        } else {
          A <- array(
            param[ARMA_indices],
            dim = c(n_sarima, n_sarima, sarima_spec$ar[[i]] + sarima_spec$ma[[i]])
          )
        }
        coeff <- CoeffARMA(
          A = A, variance = result$Q,
          ar = sarima_spec$ar[[i]], ma = sarima_spec$ma[[i]]
        )

        # Update AR coefficient array
        if (sarima_spec$ar[[i]] > 0) {
          if (decompositions) {
            sar[[paste0("S", sarima_spec$s[[i]])]] <- coeff$ar
          }
          if (!transform_only) {
            AR_coeff_old <- AR_coeff
            add_coeff <- sarima_spec$s[[i]] * 1:sarima_spec$ar[[i]]
            if (n_sarima == 1) {
              AR_coeff[add_coeff] <- AR_coeff[add_coeff] - coeff$ar
              for (j in (1:sarima_spec$ar[[i]])) {
                mult_coeff <- sarima_spec$s[[i]] * j + ar_polynomial
                AR_coeff[mult_coeff] <- AR_coeff[mult_coeff] -
                  coeff$ar[[j]] * AR_coeff_old[ar_polynomial]
                add_coeff <- c(add_coeff, mult_coeff)
              }
            } else {
              AR_coeff[, , add_coeff] <- AR_coeff[, , add_coeff, drop = FALSE] - coeff$ar
              for (j in (1:sarima_spec$ar[[i]])) {
                mult_coeff <- sarima_spec$s[[i]] * j + ar_polynomial
                AR_coeff[, , mult_coeff] <- AR_coeff[, , mult_coeff, drop = FALSE] -
                  array(
                    coeff$ar[, , j] %*%
                      matrix(AR_coeff_old[, , ar_polynomial], nrow = n_sarima),
                    dim = c(n_sarima, n_sarima, length(ar_polynomial))
                  )
                add_coeff <- c(add_coeff, mult_coeff)
              }
            }
            ar_polynomial <- unique(c(ar_polynomial, add_coeff))
          }
        }

        # Update MA coefficient array
        if (sarima_spec$ma[[i]] > 0) {
          if (decompositions) {
            sma[[paste0("S", sarima_spec$s[[i]])]] <- coeff$ma
          }
          if (!transform_only) {
            MA_coeff_old <- MA_coeff
            add_coeff <- sarima_spec$s[[i]] * 1:sarima_spec$ma[[i]]
            if (n_sarima == 1) {
              MA_coeff[add_coeff] <- MA_coeff[add_coeff] + coeff$ma
              for (j in (1:sarima_spec$ma[[i]])) {
                mult_coeff <- sarima_spec$s[[i]] * j + ma_polynomial
                MA_coeff[mult_coeff] <- MA_coeff[mult_coeff] +
                  coeff$ma[[j]] * MA_coeff_old[ma_polynomial]
                add_coeff <- c(add_coeff, mult_coeff)
              }
            } else {
              MA_coeff[, , add_coeff] <- MA_coeff[, , add_coeff, drop = FALSE] + coeff$ma
              for (j in (1:sarima_spec$ma[[i]])) {
                mult_coeff <- sarima_spec$s[[i]] * j + ma_polynomial
                MA_coeff[, , mult_coeff] <- MA_coeff[, , mult_coeff, drop = FALSE] +
                  array(
                    coeff$ma[, , j] %*%
                      matrix(MA_coeff_old[, , ma_polynomial], nrow = n_sarima),
                    dim = c(n_sarima, n_sarima, length(ma_polynomial))
                  )
                add_coeff <- c(add_coeff, mult_coeff)
              }
            }
            ma_polynomial <- unique(c(ma_polynomial, add_coeff))
          }
        }

        # Drop parameters that are already used
        param <- param[-ARMA_indices]
      }

      if (!transform_only) {

        # T matrix coefficients
        if (sum(sarima_spec$ar) > 0) {
          if (n_sarima > 1) {
            AR_coeff <- aperm(AR_coeff, c(2, 1, 3))
          }
          T1[1:(n_sarima * sum(sarima_spec$s * sarima_spec$ar)), 1:n_sarima] <- t(
            matrix(
              -AR_coeff,
              n_sarima,
              n_sarima * sum(sarima_spec$s * sarima_spec$ar)
            )
          )
        }
        T_stationary <- T1

        # R matrix coefficients
        if (sum(sarima_spec$ma) > 0) {
          if (n_sarima > 1) {
            MA_coeff <- aperm(MA_coeff, c(2, 1, 3))
          }
          R1 <- rbind(
            R1,
            t(
              matrix(
                MA_coeff,
                n_sarima,
                n_sarima * sum(sarima_spec$s * sarima_spec$ma)
              )
            )
          )
        }
        R1 <- rbind(
          R1,
          matrix(
            0,
            (r - 1 - sum(sarima_spec$s * sarima_spec$ma)) * n_sarima,
            n_sarima
          )
        )
        R_stationary <- R1

        # Non-stationary part
        if (sum(sarima_spec$i) > 0) {
          T1 <- cbind(T2, T1)
          T1 <- rbind(T3, T1)
          R1 <- rbind(R2, R1)
        }
        result$Tmat <- T1
        result$R <- R1
      }

      # Return AR and MA coefficients
      if (decompositions) {
        result$sar <- sar
        result$sma <- sma
      }
    }

    if (!transform_only) {

      # Initial uncertainty for the stationary part
      if (sum(sarima_spec$ar) == 0 & sum(sarima_spec$ma) == 0) {
        result$P_star <- BlockMatrix(
          matrix(
            0,
            sum(sarima_spec$s * sarima_spec$i) * n_sarima,
            sum(sarima_spec$s * sarima_spec$i) * n_sarima
          ),
          result$Q
        )
      } else {
        T_kronecker <- kronecker(T_stationary, T_stationary)
        Tinv <- solve(diag(1, dim(T_kronecker)[[1]], dim(T_kronecker)[[2]]) - T_kronecker)
        vecRQR <- matrix(tcrossprod(R_stationary %*% result$Q, R_stationary))
        vecPstar <- Tinv %*% vecRQR
        result$P_star <- BlockMatrix(
          matrix(
            0,
            sum(sarima_spec$s * sarima_spec$i) * n_sarima,
            sum(sarima_spec$s * sarima_spec$i) * n_sarima
          ),
          matrix(vecPstar, dim(T_stationary)[[1]], dim(T_stationary)[[2]])
        )
      }
    }
  }

  return(result)
}
