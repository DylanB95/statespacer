#' Transform arbitrary matrices into partial autocorrelation matrices
#' 
#' Creates valid partial autocorrelation matrices.
#' 
#' @param A An array of arbitrary square matrices in the multivariate case,
#'   or a vector of arbitrary numbers in the univariate case.
#' 
#' @return An array of partial autocorrelation matrices in the multivariate case,
#'   or a vector of partial autocorrelations in the univariate case.
#' 
#' @noRd
PACMat <- function(A) {
  
  # If univariate, computations are less expensive
  if (is.vector(A)) {
    return(A / sqrt(1 + A^2))
  }
  
  # Multivariate, for a single partial autocorrelation matrix
  if (is.matrix(A)) {
    return(solve(t(chol(diag(dim(A)[1]) + A %*% t(A)))) %*% A)
  }
  
  # Multivariate, for multiple partial autocorrelation matrices
  return(array(apply(A, 3, PACMat), dim = dim(A)))
}

#' Transform partial autocorrelation matrices into coefficient matrices
#' 
#' @description
#' Creates coefficient matrices for which the characteristic polynomial
#' corresponds to a stationary process. Note, this function covers a 
#' subset of the space of coefficient matrices that correspond to a 
#' stationary process.
#' 
#' @param P An array of partial autocorrelation matrices in the multivariate case,
#'   or a vector of partial autocorrelations in the univariate case.
#' 
#' @return
#' Multivariate case:
#'   A list containing:
#'     An array of coefficient matrices.
#'     Variance - covariance matrix of the noise, in order to transform 
#'       coefficients further.
#' Univariate case:
#'   A vector of coefficients.
#' 
#' @noRd
TransformPAC <- function(P) {
  
  coeff_new <- P

  # If univariate, computations are less expensive
  if (is.vector(P)) {
    if (length(P) > 1) {
      for (i in 1:(length(P) - 1)) {
        coeff_old <- coeff_new
        for (j in 1:i) {
          coeff_new[j] <- coeff_old[j] - coeff_new[i + 1] * coeff_old[i - j + 1]
        }
      }
    }
    return(coeff_new)
  }
  
  # Multivariate
  # Initialise with i = 0
  sigma_new <- diag(dim(P)[1]) - coeff_new[,,1] %*% t(coeff_new[,,1])
  
  if (dim(P)[3] > 1) {
    
    # i = 0
    coeff_star_new <- P
    coeff_star_new[,,1] <- t(P[,,1])
    sigma_star_new <- diag(dim(P)[1]) - coeff_star_new[,,1] %*% t(coeff_star_new[,,1])
    L <- t(chol(sigma_new))
    L_star <- t(chol(sigma_star_new))
    
    for (i in 1:(dim(P)[3] - 1)) {
  
      # Storing former values
      coeff_old <- coeff_new
      coeff_star_old <- coeff_star_new
      sigma_old <- sigma_new
      sigma_star_old <- sigma_star_new
      
      # Calculating new values
      coeff_new[,,i + 1] <- L %*% P[,,i + 1] %*% solve(L_star)
      sigma_new <- sigma_old - coeff_new[,,i + 1] %*% sigma_star_old %*% t(coeff_new[,,i + 1])
      if (i < (dim(P)[3] - 1)) {
        coeff_star_new[,,i + 1] <- L_star %*% t(P[,,i + 1]) %*% solve(L)
        sigma_star_new <- sigma_star_old - coeff_star_new[,,i + 1] %*% sigma_old %*% t(coeff_star_new[,,i + 1])
        L_star <- t(chol(sigma_star_new))
        L <- t(chol(sigma_new))
      }
      for (j in 1:i) {
        coeff_new[,,j] <- coeff_old[,,j] - coeff_new[,,i + 1] %*% coeff_star_old[,,i - j + 1]
        if (i < (dim(P)[3] - 1)) {
          coeff_star_new[,,j] <- coeff_star_old[,,j] - coeff_star_new[,,i + 1] %*% coeff_old[,,i - j + 1]
        }
      }
    }
  }
  return(list(coeff = coeff_new, variance = sigma_new))
}

#' Transform arbitrary matrices into ARMA coefficient matrices
#' 
#' @description
#' Creates coefficient matrices for which the characteristic polynomial
#' corresponds to a stationary process.
#' 
#' @param A An array of arbitrary square matrices in the multivariate case,
#'   or a vector of arbitrary numbers in the univariate case.
#' @param variance A variance - covariance matrix in the multivariate case,
#'   or the variance in the univariate case.
#' 
#' @return A list containing:
#'   Multivariate case:
#'     An array of coefficient matrices for the AR part.
#'     An array of coefficient matrices for the MA part.
#'   Univariate case:
#'     A vector of coefficients for the AR part.
#'     A vector of coefficients for the MA part.
#' 
#' @noRd
CoeffARMA <- function(A, variance = NULL, ar = 1, ma = 0) {
  
  # Check for erroneous input
  if (ar < 0) {
    stop("The order of the autoregressive part must be >= 0.")
  }
  if (ma < 0) {
    stop("The order of the moving average part must be >= 0.")
  }
  if (ar + ma == 0) {
    stop("At least one of the orders of the AR and MA parts must be positive.")
  }
  
  # Initialise list to return
  result <- list()
  
  # Obtain partial autocorrelation matrices
  P <- PACMat(A)
  
  # If univariate, coefficients are readily returned by TransformPAC
  if (is.vector(P)) {
    if (ar > 0) {
      P_ar <- P[1:ar]
      result$ar <- TransformPAC(P_ar)
    }
    if (ma > 0) {
      P_ma <- P[(ar + 1):(ar + ma)]
      result$ma <- -TransformPAC(P_ma) # Note the minus sign
    }
    return(result)
  }
  
  # Cholesky decomposition of variance - covariance matrix
  L <- t(chol(variance))
  L_inv <- solve(L)
  
  # Obtain coefficient matrices for AR part
  if (ar > 0) {
    P_ar <- P[,,1:ar, drop = FALSE]
    ar_part <- TransformPAC(P_ar)
    L_ar <- t(chol(ar_part$variance))
    L_ar_inv <- solve(L_ar)
    result$ar <- array(
      apply(ar_part$coeff, 3, function(x) L %*% L_ar_inv %*% x %*% L_ar %*% L_inv),
      dim = dim(ar_part$coeff)
    )
  }
  
  # Obtain coefficient matrices for MA part
  if (ma > 0) {
    P_ma <- P[,,(ar + 1):(ar + ma), drop = FALSE]
    ma_part <- TransformPAC(P_ma)
    L_ma <- t(chol(ma_part$variance))
    L_ma_inv <- solve(L_ma)
    result$ma <- -array( # Note the minus sign
      apply(ma_part$coeff, 3, function(x) L %*% L_ma_inv %*% x %*% L_ma %*% L_inv),
      dim = dim(ma_part$coeff)
    )
  }
  return(result)
}