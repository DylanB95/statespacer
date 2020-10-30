#' Generating Random Samples using the Simulation Smoother
#'
#' Draws random samples of the specified model conditional
#' on the observed data.
#'
#' @param nsim Number of random samples to draw. Defaults to `1`.
#' @param components Boolean indicating whether the components of
#'   the model should be extracted in each of the random samples.
#' @inheritParams predict.statespacer
#'
#' @return
#' A list containing the simulated state parameters and disturbances.
#' In addition, it returns the components as specified
#' by the State Space model.
#'
#' @author Dylan Beijers, \email{dylanbeijers@@gmail.com}
#'
#' @examples
#' # Fits a local level model for the Nile data
#' library(datasets)
#' y <- matrix(Nile)
#' fit <- statespacer(initial = 10, y = y, local_level_ind = TRUE)
#'
#' # Obtain random sample using the fitted model
#' #sim <- SimSmoother(fit, nsim = 1, components = TRUE)
#'
#' # Plot the simulated level
#' #plot(1:10, sim$level[ , 1, 1], type = "l")
#' @export
SimSmoother <- function(object,
                        nsim = 1,
                        components = TRUE) {

  # Initialising list to return
  result <- list()

  # Number of observations
  N <- dim(object$function_call$y)[[1]]

  # Number of dependent variables
  p <- dim(object$function_call$y)[[2]]

  # Number of state parameters
  m <- dim(object$smoothed$a)[[2]]

  # Number of state disturbances
  r <- dim(object$smoothed$eta)[[2]]

  # Initialising arrays
  y <- array(0, dim = c(N, p, nsim))
  eta <- array(0, dim = c(N, r, nsim))
  a <- array(0, dim = c(N, m, nsim))

  # Current last filled column indices of a and eta
  eta_index <- 0
  a_index <- 0

  #### Local Level ####
  if (object$function_call$local_level_ind && !object$function_call$slope_ind &&
      is.null(object$function_call$level_addvar_list)) {

    # Model components
    Z <- object$system_matrices$Z$level
    Tmat <- object$system_matrices$T$level
    R <- object$system_matrices$R$level
    Q <- object$system_matrices$Q$level
    a1 <- object$system_matrices$a1$level
    dim(Z) <- c(dim(Z), 1)
    dim(Tmat) <- c(dim(Tmat), 1)
    dim(R) <- c(dim(R), 1)
    dim(Q) <- c(dim(Q), 1)

    # Indices of eta and a components
    eta_num <- dim(R)[[2]]
    a_num <- dim(Z)[[2]]
    eta_indices <- eta_index + 1:eta_num
    a_indices <- a_index + 1:a_num
    eta_index <- eta_index + eta_num
    a_index <- a_index + a_num

    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = N,
      a = a1,
      Z = Z,
      T = Tmat,
      R = R,
      Q = Q,
      P_star = matrix(0),
      draw_initial = FALSE,
      eta_only = FALSE
    )
    y <- y + sim$y
    eta[ , eta_indices, ] <- sim$eta
    a[ , a_indices, ] <- sim$a
  }

  #### Local Level + Slope ####
  if (object$function_call$slope_ind &&
      is.null(object$function_call$level_addvar_list)) {

    # Model components
    Z <- object$system_matrices$Z$level
    Tmat <- object$system_matrices$T$level
    R <- object$system_matrices$R$level
    Q <- BlockMatrix(
      object$system_matrices$Q$level,
      object$system_matrices$Q$slope
    )
    a1 <- object$system_matrices$a1$level
    dim(Z) <- c(dim(Z), 1)
    dim(Tmat) <- c(dim(Tmat), 1)
    dim(R) <- c(dim(R), 1)
    dim(Q) <- c(dim(Q), 1)

    # Indices of eta and a components
    eta_num <- dim(R)[[2]]
    a_num <- dim(Z)[[2]]
    eta_indices <- eta_index + 1:eta_num
    a_indices <- a_index + 1:a_num
    eta_index <- eta_index + eta_num
    a_index <- a_index + a_num

    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = N,
      a = a1,
      Z = Z,
      T = Tmat,
      R = R,
      Q = Q,
      P_star = matrix(0),
      draw_initial = FALSE,
      eta_only = FALSE
    )
    y <- y + sim$y
    eta[ , eta_indices, ] <- sim$eta
    a[ , a_indices, ] <- sim$a
  }

  #### BSM ####
  if (length(object$function_call$BSM_vec) > 0) {
    for (s in object$function_call$BSM_vec) {

      # Model components
      Z <- object$system_matrices$Z[[paste0("BSM", s)]]
      Tmat <- object$system_matrices$T[[paste0("BSM", s)]]
      R <- object$system_matrices$R[[paste0("BSM", s)]]
      Q <- object$system_matrices$Q[[paste0("BSM", s)]]
      a1 <- object$system_matrices$a1[[paste0("BSM", s)]]
      dim(Z) <- c(dim(Z), 1)
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)

      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num

      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = s - 1,
        N = N,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = matrix(0),
        draw_initial = FALSE,
        eta_only = FALSE
      )
      y <- y + sim$y
      eta[ , eta_indices, ] <- sim$eta
      a[ , a_indices, ] <- sim$a
    }
  }

  #### Explanatory Variables ####
  if (!is.null(object$function_call$addvar_list)) {

    # Model components
    Z <- object$system_matrices$Z$addvar
    Tmat <- object$system_matrices$T$addvar
    R <- object$system_matrices$R$addvar
    Q <- object$system_matrices$Q$addvar
    a1 <- object$system_matrices$a1$addvar
    dim(Tmat) <- c(dim(Tmat), 1)
    dim(R) <- c(dim(R), 1)
    dim(Q) <- c(dim(Q), 1)

    # Indices of eta and a components
    eta_num <- dim(R)[[2]]
    a_num <- dim(Z)[[2]]
    eta_indices <- eta_index + 1:eta_num
    a_indices <- a_index + 1:a_num
    eta_index <- eta_index + eta_num
    a_index <- a_index + a_num

    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = N,
      a = a1,
      Z = Z,
      T = Tmat,
      R = R,
      Q = Q,
      P_star = matrix(0),
      draw_initial = FALSE,
      eta_only = FALSE
    )
    y <- y + sim$y
    eta[ , eta_indices, ] <- sim$eta
    a[ , a_indices, ] <- sim$a
  }

  #### Local Level + Explanatory Variables ####
  if (!is.null(object$function_call$level_addvar_list) &&
      !object$function_call$slope_ind) {

    # Model components
    Z <- object$system_matrices$Z$level
    Tmat <- object$system_matrices$T$level
    R <- object$system_matrices$R$level
    Q <- BlockMatrix(
      object$system_matrices$Q$level,
      object$system_matrices$Q$level_addvar
    )
    a1 <- object$system_matrices$a1$level
    dim(Z) <- c(dim(Z), 1)
    dim(R) <- c(dim(R), 1)
    dim(Q) <- c(dim(Q), 1)

    # Indices of eta and a components
    eta_num <- dim(R)[[2]]
    a_num <- dim(Z)[[2]]
    eta_indices <- eta_index + 1:eta_num
    a_indices <- a_index + 1:a_num
    eta_index <- eta_index + eta_num
    a_index <- a_index + a_num

    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = N,
      a = a1,
      Z = Z,
      T = Tmat,
      R = R,
      Q = Q,
      P_star = matrix(0),
      draw_initial = FALSE,
      eta_only = FALSE
    )
    y <- y + sim$y
    eta[ , eta_indices, ] <- sim$eta
    a[ , a_indices, ] <- sim$a
  }

  #### Local Level + Explanatory Variables + Slope ####
  if (!is.null(object$function_call$level_addvar_list) &&
      object$function_call$slope_ind) {

    # Model components
    Z <- object$system_matrices$Z$level
    Tmat <- object$system_matrices$T$level
    R <- object$system_matrices$R$level
    Q <- BlockMatrix(
      object$system_matrices$Q$level,
      object$system_matrices$Q$slope,
      object$system_matrices$Q$level_addvar
    )
    a1 <- object$system_matrices$a1$level
    dim(Z) <- c(dim(Z), 1)
    dim(R) <- c(dim(R), 1)
    dim(Q) <- c(dim(Q), 1)

    # Indices of eta and a components
    eta_num <- dim(R)[[2]]
    a_num <- dim(Z)[[2]]
    eta_indices <- eta_index + 1:eta_num
    a_indices <- a_index + 1:a_num
    eta_index <- eta_index + eta_num
    a_index <- a_index + a_num

    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = N,
      a = a1,
      Z = Z,
      T = Tmat,
      R = R,
      Q = Q,
      P_star = matrix(0),
      draw_initial = FALSE,
      eta_only = FALSE
    )
    y <- y + sim$y
    eta[ , eta_indices, ] <- sim$eta
    a[ , a_indices, ] <- sim$a
  }

  #### Cycle ####
  if (object$function_call$cycle_ind) {
    for (i in seq_along(object$function_call$format_cycle_list)) {

      # Model components
      Z <- object$system_matrices$Z[[paste0("Cycle", i)]]
      Tmat <- object$system_matrices$T[[paste0("Cycle", i)]]
      R <- object$system_matrices$R[[paste0("Cycle", i)]]
      Q <- object$system_matrices$Q[[paste0("Cycle", i)]]
      a1 <- object$system_matrices$a1[[paste0("Cycle", i)]]
      P_star <- object$system_matrices$P_star[[paste0("Cycle", i)]]
      dim(Z) <- c(dim(Z), 1)
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)
      draw_initial <- any(P_star > 0)

      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num

      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 2,
        N = N,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P_star,
        draw_initial = draw_initial,
        eta_only = FALSE
      )
      y <- y + sim$y
      eta[ , eta_indices, ] <- sim$eta
      a[ , a_indices, ] <- sim$a
    }
  }

  #### ARIMA ####
  if (!is.null(object$function_call$arima_list)) {
    for (i in seq_along(object$function_call$arima_list)) {

      # Model components
      Z <- object$system_matrices$Z[[paste0("ARIMA", i)]]
      Tmat <- object$system_matrices$T[[paste0("ARIMA", i)]]
      R <- object$system_matrices$R[[paste0("ARIMA", i)]]
      Q <- object$system_matrices$Q[[paste0("ARIMA", i)]]
      a1 <- object$system_matrices$a1[[paste0("ARIMA", i)]]
      P_star <- object$system_matrices$P_star[[paste0("ARIMA", i)]]
      dim(Z) <- c(dim(Z), 1)
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)

      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num

      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = N,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P_star,
        draw_initial = TRUE,
        eta_only = FALSE
      )
      y <- y + sim$y
      eta[ , eta_indices, ] <- sim$eta
      a[ , a_indices, ] <- sim$a
    }
  }

  #### SARIMA ####
  if (!is.null(object$function_call$sarima_list)) {
    for (i in seq_along(object$function_call$sarima_list)) {

      # Model components
      Z <- object$system_matrices$Z[[paste0("SARIMA", i)]]
      Tmat <- object$system_matrices$T[[paste0("SARIMA", i)]]
      R <- object$system_matrices$R[[paste0("SARIMA", i)]]
      Q <- object$system_matrices$Q[[paste0("SARIMA", i)]]
      a1 <- object$system_matrices$a1[[paste0("SARIMA", i)]]
      P_star <- object$system_matrices$P_star[[paste0("SARIMA", i)]]
      dim(Z) <- c(dim(Z), 1)
      dim(Tmat) <- c(dim(Tmat), 1)
      dim(R) <- c(dim(R), 1)
      dim(Q) <- c(dim(Q), 1)

      # Indices of eta and a components
      eta_num <- dim(R)[[2]]
      a_num <- dim(Z)[[2]]
      eta_indices <- eta_index + 1:eta_num
      a_indices <- a_index + 1:a_num
      eta_index <- eta_index + eta_num
      a_index <- a_index + a_num

      # Obtain random samples of component
      sim <- SimulateC(
        nsim = nsim,
        repeat_Q = 1,
        N = N,
        a = a1,
        Z = Z,
        T = Tmat,
        R = R,
        Q = Q,
        P_star = P_star,
        draw_initial = TRUE,
        eta_only = FALSE
      )
      y <- y + sim$y
      eta[ , eta_indices, ] <- sim$eta
      a[ , a_indices, ] <- sim$a
    }
  }

  #### Self Specified ####
  if (!is.null(object$function_call$self_spec_list)) {

    # Model components
    Z <- object$system_matrices$Z$self_spec
    Tmat <- object$system_matrices$T$self_spec
    R <- object$system_matrices$R$self_spec
    Q <- object$system_matrices$Q$self_spec
    a1 <- object$system_matrices$a1$self_spec
    P_star <- object$system_matrices$P_star$self_spec
    if (is.matrix(Z)) {
      dim(Z) <- c(dim(Z), 1)
    }
    if (is.matrix(Tmat)) {
      dim(Tmat) <- c(dim(Tmat), 1)
    }
    if (is.matrix(R)) {
      dim(R) <- c(dim(R), 1)
    }
    if (is.matrix(Q)) {
      dim(Q) <- c(dim(Q), 1)
    }
    draw_initial <- any(P_star > 0)

    # Indices of eta and a components
    eta_num <- dim(R)[[2]]
    a_num <- dim(Z)[[2]]
    eta_indices <- eta_index + 1:eta_num
    a_indices <- a_index + 1:a_num
    eta_index <- eta_index + eta_num
    a_index <- a_index + a_num

    # Obtain random samples of component
    sim <- SimulateC(
      nsim = nsim,
      repeat_Q = 1,
      N = N,
      a = a1,
      Z = Z,
      T = Tmat,
      R = R,
      Q = Q,
      P_star = P_star,
      draw_initial = draw_initial,
      eta_only = FALSE
    )
    y <- y + sim$y
    eta[ , eta_indices, ] <- sim$eta
    a[ , a_indices, ] <- sim$a
  }

  #### Residuals ####

  # Model components
  Q <- object$system_matrices$H$H
  if (is.matrix(Q)) {
    dim(Q) <- c(dim(Q), 1)
  }

  # Obtain random samples of component
  sim <- SimulateC(
    nsim = nsim,
    repeat_Q = 1,
    N = N,
    a = matrix(0),
    Z = array(0, dim = c(1, 1, 1)),
    T = array(0, dim = c(1, 1, 1)),
    R = array(0, dim = c(1, 1, 1)),
    Q = Q,
    P_star = matrix(0),
    draw_initial = FALSE,
    eta_only = TRUE
  )
  epsilon <- sim$eta
  y <- y + epsilon

  #### Obtain smoothed state and disturbances of random samples ####

  # Model components
  Z <- object$system_matrices$Z$full
  Tmat <- object$system_matrices$T$full
  R <- object$system_matrices$R$full
  Q <- object$system_matrices$Q$full
  a1 <- object$system_matrices$a1$full
  P_star <- object$system_matrices$P_star$full
  P_inf <- object$system_matrices$P_inf$full
  H <- object$system_matrices$H$H

  # Add residuals to system matrices for univariate treatment
  Z <- CombineZ(diag(1, p, p), Z)
  Tmat <- CombineTRQ(matrix(0, p, p), Tmat)
  R <- CombineTRQ(diag(1, p, p), R)
  if (is.matrix(H)) {
    Q <- CombineTRQ(H, Q)
    P_star <- BlockMatrix(H, P_star)
  } else {
    Q <- CombineTRQ(H[, , c(2:dim(H)[[3]], 1)], Q)
    P_star <- BlockMatrix(H[, , 1], P_star)
  }
  a1 <- rbind(matrix(0, p, 1), a1)
  P_inf <- BlockMatrix(matrix(0, p, p), P_inf)

  # Augment dimensions of system matrices for compatibility with FastSmoother
  if (is.matrix(Z)) {
    dim(Z) <- c(dim(Z), 1)
  }
  if (is.matrix(Tmat)) {
    dim(Tmat) <- c(dim(Tmat), 1)
  }
  if (is.matrix(R)) {
    dim(R) <- c(dim(R), 1)
  }
  if (is.matrix(Q)) {
    dim(Q) <- c(dim(Q), 1)
  }

  # Smoothed state and disturbances
  fast_smoother <- FastSmootherC(
    y = y,
    a = a1,
    P_inf = P_inf,
    P_star = P_star,
    Z = Z,
    T = Tmat,
    R = R,
    Q = Q,
    initialisation_steps = object$diagnostics$initialisation_steps,
    transposed_state = components
  )
  a_smooth <- fast_smoother$a_smooth[ , -(1:p), , drop = FALSE]
  epsilon_smooth <- fast_smoother$a_smooth[ , 1:p, , drop = FALSE]
  eta_smooth <- fast_smoother$eta[ , -(1:p), , drop = FALSE]
  if (components) {
    a_t <- fast_smoother$a_t[-(1:p), , , drop = FALSE]
  }

  # Adjust random samples by mean corrections
  a <- a - a_smooth + array(object$smoothed$a, dim = c(N, m, nsim))
  epsilon <- epsilon - epsilon_smooth +
    array(object$smoothed$epsilon, dim = c(N, p, nsim))
  eta <- eta - eta_smooth + array(object$smoothed$eta, dim = c(N, r, nsim))

  # Return the list containing the simulated samples
  result$a <- a
  result$epsilon <- epsilon
  result$eta <- eta
  return(result)
}
