#' Kalman Filter with Univariate Treatment
#'
#' Applies the kalman filter with univariate treatment after the
#' initialisation stages to calculate the estimated state and
#' corresponding variance for the next step.
#'
#' @param y Scalar containing the dependent variable.
#' @param a Column vector containing the state parameters.
#' @param P Variance - covariance matrix of the state vector.
#' @param Z Z system matrix of the State Space model.
#' @param Tmat T system matrix of the State Space model.
#' @param R R system matrix of the State Space model.
#' @param Q Q system matrix of the State Space model.
#' @param timestep Boolean indicating whether a transition to the next
#'   timepoint should be made.
#'
#' @return
#' A list containing:
#' * `a`: The estimated state for the next step.
#' * `P`: The corresponding variance - covariance matrix of the estimated state
#'   for the next step.
#' * `loglik`: The loglikelihood for the current step.
#' * `a_fil`: The filtered state for the current step.
#' * `P_fil`: The corresponding variance - covariance matrix of the filtered
#'   state for the current step.
#'
#' @noRd
KalmanUT <- function(y, a, P, Z, Tmat = NULL, R = NULL, Q = NULL, timestep) {

  # Initialise list to return
  result <- list()

  # If y is missing, use kalman filter formulae with Z = 0
  if (is.na(y)) {

    # Estimated state vector and corresponding variance - covariance matrix
    # for the next step
    a_new <- a
    P_new <- P

    # If transition to next timepoint, do some additional computations
    if (timestep) {

      # Return filtered state and variance
      result$a_fil <- a_new
      result$P_fil <- P_new

      # Estimated state vector and corresponding variance - covariance matrix
      # for the next step
      a_new <- Tmat %*% a_new
      P_new <- Tmat %*% P_new %*% t(Tmat) + R %*% Q %*% t(R)
    }

    # Loglik not available
    loglik <- NA

    # Adding predicted state, variance and loglikelihood to the list that
    # will be returned
    result$a <- a_new
    result$P <- P_new
    result$loglik <- loglik

    # Return the list
    return(result)
  }

  # PZ' as in Kalman formulae
  PtZ <- P %*% t(Z)

  # Variance matrix of the current residual/fitted value
  Fmat <- c(Z %*% PtZ)

  # Check if Fmat is nearly 0
  if (Fmat < 1e-7 & timestep) {

    # No new information available in this step
    a_new <- a
    P_new <- P
    loglik <- NA

  } else {

    # Inverse of Fmat
    Finv <- 1 / Fmat

    # Current residual
    v <- y - c(Z %*% a)

    # Loglikelihood
    loglik <- -0.5 * (log(2*pi) + log(Fmat) + v^2 * Finv)

    # Kernel matrix
    K <- PtZ * Finv

    # Estimated state vector and corresponding variance - covariance matrix
    # for the next step
    a_new <- a + K * v
    P_new <- P - K %*% t(K) * Fmat
  }

  # If transition to next timepoint, do some additional multiplications
  if (timestep) {

    # Return filtered state and variance
    result$a_fil <- a_new
    result$P_fil <- P_new

    # Estimated state vector and corresponding variance - covariance matrix
    # for the next step
    a_new <- Tmat %*% a_new
    P_new <- Tmat %*% P_new %*% t(Tmat) + R %*% Q %*% t(R)
  }

  # Adding predicted state, variance and loglikelihood to the list that
  # will be returned
  result$a <- a_new
  result$P <- P_new
  result$loglik <- loglik

  # Return the list
  return(result)
}

#' Kalman Filter with Exact Initialisation
#'
#' Applies an exact initialisation during the initialisation stages
#' for the kalman filter with univariate treatment to calculate the
#' estimated state and corresponding variance for the next step.
#'
#' @param P_inf Diffuse part of the variance - covariance matrix of the state
#'   vector.
#' @param P_star Stationary part of the variance - covariance matrix of the
#'   state vector.
#' @inheritParams KalmanUT
#'
#' @return
#' A list containing:
#' * `a`: The estimated state for the next step.
#' * `P_inf`: The corresponding diffuse part of the variance - covariance
#'   matrix of the estimated state for the next step.
#' * `P_star`: The corresponding stationary part of the variance - covariance
#'   matrix of the estimated state for the next step.
#' * `loglik`: The loglikelihood for the current step.
#' * `a_fil`: The filtered state for the current step.
#' * `P_inf_fil`: The corresponding diffuse part of the variance - covariance
#'   matrix of the filtered state for the current step.
#' * `P_star_fil`: The corresponding stationary part of the variance -
#'   covariance matrix of the filtered state for the current step.
#'
#' @noRd
KalmanEI <- function(y, a, P_inf, P_star, Z,
                     Tmat = NULL, R = NULL, Q = NULL, timestep) {

  # Initialise list to return
  result <- list()

  # If y is missing, use kalman filter formulae with Z = 0
  if (is.na(y)) {

    # Estimated state vector and corresponding variance - covariance matrix
    # for the next step
    a_new <- a
    P_inf_new <- P_inf
    P_star_new <- P_star

    # If transition to next timepoint, do some additional computations
    if (timestep) {

      # Return filtered state and variance
      result$a_fil <- a_new
      result$P_inf_fil <- P_inf_new
      result$P_star_fil <- P_star_new

      # Estimated state vector and corresponding variance - covariance matrix
      # for the next step
      a_new <- Tmat %*% a_new
      P_inf_new <- Tmat %*% P_inf_new %*% t(Tmat)
      P_star_new <- Tmat %*% P_star_new %*% t(Tmat) + R %*% Q %*% t(R)
    }

    # Loglik not available
    loglik <- NA

    # Adding predicted state, variance and loglikelihood to the list that
    # will be returned
    result$a <- a_new
    result$P_inf <- P_inf_new
    result$P_star <- P_star_new
    result$loglik <- loglik

    # Return the list
    return(result)
  }

  # PZ' as in Kalman formulae
  M_inf <- P_inf %*% t(Z)
  M_star <- P_star %*% t(Z)

  # Variance matrix of the current residual/fitted value
  F_inf <- c(Z %*% M_inf)
  F_star <- c(Z %*% M_star)

  # Check if F_inf is nearly 0
  if (F_inf < 1e-7) {

    # Check if F_star is nearly 0
    if (F_star < 1e-7 & timestep) {

      # No new information available in this step
      a_new <- a
      P_inf_new <- P_inf
      P_star_new <- P_star
      loglik <- NA

    } else {

      # Inverse of Fmat
      F_1 <- 1 / F_star

      # Current residual
      v <- y - c(Z %*% a)

      # Auxiliary matrices
      K_0 <- M_star * F_1
      L_0 <- diag(length(Z)) - K_0 %*% Z

      # Estimated state vector and corresponding variance - covariance matrix
      # for the next step
      a_new <- a + K_0 * v
      P_inf_new <- P_inf
      P_star_new <- P_star %*% t(L_0)

      # Loglikelihood
      loglik <- -0.5 * (log(2*pi) + log(F_star) + v^2 * F_1)
    }

  } else {

    # Inverse of Fmat
    F_1 <- 1 / F_inf
    F_2 <- -F_1 * F_star * F_1

    # Current residual
    v <- y - c(Z %*% a)

    # Auxiliary matrices
    K_0 <- M_inf * F_1
    L_0 <- diag(length(Z)) - K_0 %*% Z
    K_1 <- M_star * F_1 + M_inf * F_2
    L_1 <- -K_1 %*% Z

    # Estimated state vector and corresponding variance - covariance matrix
    # for the next step
    a_new <- a + K_0 * v
    P_inf_new <- P_inf %*% t(L_0)
    P_star_new <- P_inf %*% t(L_1) + P_star %*% t(L_0)

    # Loglikelihood
    loglik <- -0.5 * (log(2*pi) + log(F_inf))
  }

  # If transition to next timepoint, do some additional multiplications
  if (timestep) {

    # Return filtered state and variance
    result$a_fil <- a_new
    result$P_inf_fil <- P_inf_new
    result$P_star_fil <- P_star_new

    # Estimated state vector and corresponding variance - covariance matrix
    # for the next step
    a_new <- Tmat %*% a_new
    P_inf_new <- Tmat %*% P_inf_new %*% t(Tmat)
    P_star_new <- Tmat %*% P_star_new %*% t(Tmat) + R %*% Q %*% t(R)
  }

  # Adding predicted state, variance and loglikelihood to the list that
  # will be returned
  result$a <- a_new
  result$P_inf <- P_inf_new
  result$P_star <- P_star_new
  result$loglik <- loglik

  # Return the list
  return(result)
}
