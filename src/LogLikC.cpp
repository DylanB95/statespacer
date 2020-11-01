#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]

//' Employ the Kalman Filter to calculate the loglikelihood
//'
//' @param y Matrix of observations.
//' @param y_isna Matrix indicating which observations are missing.
//' @param a Initial values of the state vector.
//' @param P_inf Diffuse part of the variance - covariance matrix of the
//'   state vector.
//' @param P_star Stationary part of the variance - covariance matrix of the
//'   state vector.
//' @param Z Z system matrix of the State Space model.
//' @param T T system matrix of the State Space model.
//' @param R R system matrix of the State Space model.
//' @param Q Q system matrix of the State Space model.
//'
//' @noRd
// [[Rcpp::export]]
double LogLikC(const Rcpp::NumericMatrix& y,
               const Rcpp::LogicalMatrix& y_isna,
               arma::colvec a,
               arma::mat P_inf,
               arma::mat P_star,
               const arma::cube& Z,
               const arma::cube& T,
               const arma::cube& R,
               const arma::cube& Q) {

  // Number of observations, dependent variables, and state parameters
  int N = y.nrow(), p = y.ncol(), m = a.n_rows;

  // Keep track of the limits of the indices
  int N_min1 = N - 1, p_min1 = p - 1;

  // Check which system matrices are time-varying
  bool Z_tv = Z.n_slices > 1, T_tv = T.n_slices > 1,
       R_tv = R.n_slices > 1, Q_tv = Q.n_slices > 1;

  // Initial system matrices
  arma::mat Z_mat = Z.slice(0), T_mat = T.slice(0),
            R_mat = R.slice(0), Q_mat = Q.slice(0);
  arma::rowvec Z_row = Z_mat.row(0);

  // Indicator for whether the first row should be assigned
  bool row_assign = Z_tv || p > 1;

  // Check if P_inf is already 0
  bool initialisation = !arma::all(arma::vectorise(arma::abs(P_inf)) < 1e-7);

  // Initialise number of non-initialisation steps
  // and number of occurrences that F = 0
  int non_init = 0, F_eq0 = 0;

  // Initialise loglikelihood
  double loglik = 0;

  // Initialise objects used in computations
  arma::colvec M_inf(m), M_star(m), K_0(m), K_1(m);
  arma::mat L_0(m, m), L_1(m, m);
  double F_inf, F_star, F_1, F_2, v;
  double constant = -log(2 * M_PI) / 2;

  // Iterators
  int i, j;

  // Loop over time points
  for (i = 0; i < N; i++) {

    // Get system matrices of current time point
    if (Z_tv && i > 0) {
      Z_mat = Z.slice(i);
    }
    if (T_tv && i > 0) {
      T_mat = T.slice(i);
    }
    if (R_tv && i > 0) {
      R_mat = R.slice(i);
    }
    if (Q_tv && i > 0) {
      Q_mat = Q.slice(i);
    }

    // Loop over dependent variables
    for (j = 0; j < p; j++) {

      // Check for missing value
      if (y_isna(i, j)) {
        continue;
      }

      // Retrieve row of Z
      if (j > 0 || (i > 0 && row_assign)) {
        Z_row = Z_mat.row(j);
      }

      // Exact Kalman filter in initialisation steps
      if (initialisation) {

        // PZ' as in Kalman formulae
        M_inf = P_inf * Z_row.t();
        M_star = P_star * Z_row.t();

        // Variance matrix of the current residual/fitted value
        F_inf = arma::as_scalar(Z_row * M_inf);
        F_star = arma::as_scalar(Z_row * M_star);

        // Check if F_inf is nearly 0
        if (F_inf < 1e-7) {

          // Check if F_star is nearly 0
          if (F_star < 1e-7) {
            continue;
          } else {

            // Inverse of Fmat
            F_1 = 1 / F_star;

            // Current residual
            v = y(i, j) - arma::as_scalar(Z_row * a);

            // Auxiliary matrices
            K_0 = M_star * F_1;
            L_0 = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;

            // Estimated state vector and corresponding variance - covariance
            // matrix for the next step
            a = a + K_0 * v;
            P_star = P_star * L_0.t();
            loglik += constant - (log(F_star) + pow(v, 2.0) * F_1) / 2;
            continue;
          }
        } else {

          // Inverse of Fmat
          F_1 = 1 / F_inf;
          F_2 = -pow(F_1, 2.0) * F_star;

          // Current residual
          v = y(i, j) - arma::as_scalar(Z_row * a);

          // Auxiliary matrices
          K_0 = M_inf * F_1;
          L_0 = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;
          K_1 = M_star * F_1 + M_inf * F_2;
          L_1 = -K_1 * Z_row;

          // Estimated state vector and corresponding variance - covariance
          // matrix for the next step
          a = a + K_0 * v;
          P_star = P_inf * L_1.t() + P_star * L_0.t();
          P_inf = P_inf * L_0.t();
          loglik += constant - log(F_inf) / 2;
        }

        // Check if P_inf converged to 0
        if (j < p_min1) {
          initialisation = !arma::all(arma::vectorise(arma::abs(P_inf)) < 1e-7);
        }

      } else {

        // Increment non_init
        non_init++;

        // PZ' as in Kalman formulae
        M_star = P_star * Z_row.t();

        // Variance matrix of the current residual/fitted value
        F_star = arma::as_scalar(Z_row * M_star);

        // Check if F_star is nearly 0
        if (F_star < 1e-7) {

          // Increment number of times F = 0
          F_eq0++;

        } else {

          // Inverse of Fmat
          F_1 = 1 / F_star;

          // Current residual
          v = y(i, j) - arma::as_scalar(Z_row * a);

          // Auxiliary matrices
          K_0 = M_star * F_1;
          L_0 = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;

          // Estimated state vector and corresponding variance - covariance
          // matrix for the next step
          a = a + K_0 * v;
          P_star = P_star * L_0.t();
          loglik += constant - (log(F_star) + pow(v, 2.0) * F_1) / 2;
        }
      }
    }

    // Perform computations for the next time point
    if (i < N_min1) {
      a = T_mat * a;
      P_star = T_mat * P_star * T_mat.t() + R_mat * Q_mat * R_mat.t();
      if (initialisation) {
        P_inf = T_mat * P_inf * T_mat.t();
        initialisation = !arma::all(arma::vectorise(arma::abs(P_inf)) < 1e-7);
      }
    }
  }

  // Return NA if F equals 0 during all of the non-initialisation steps
  if (non_init == F_eq0) {
    return NA_REAL;
  }
  return loglik / N;
}
