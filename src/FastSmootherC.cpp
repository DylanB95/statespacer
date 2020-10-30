#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]

//' Employ the Fast Kalman Smoother on multiple samples
//'
//' @param y Array of samples.
//' @param a Initial values of the state vector.
//' @param P_inf Diffuse part of the variance - covariance matrix of the
//'   state vector.
//' @param P_star Stationary part of the variance - covariance matrix of the
//'   state vector.
//' @param Z Z system matrix of the State Space model.
//' @param T T system matrix of the State Space model.
//' @param R R system matrix of the State Space model.
//' @param Q Q system matrix of the State Space model.
//' @param initialisation_steps Number of steps that were
//'   needed during initialisation.
//' @param transposed_state Boolean indicating whether a
//'   transposed variant of the state should be returned.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List FastSmootherC(const arma::cube& y,
                         const arma::colvec& a,
                         const arma::mat& P_inf,
                         const arma::mat& P_star,
                         const arma::cube& Z,
                         const arma::cube& T,
                         const arma::cube& R,
                         const arma::cube& Q,
                         const int& initialisation_steps,
                         const bool& transposed_state) {

  // Number of observations, dependent variables, state parameters,
  // and state disturbances
  int N = y.n_rows, p = y.n_cols, m = a.n_rows, r = R.n_cols, nsim = y.n_slices;

  // Keep track of the current index and limits of indices
  int index = 0, Np = N * p, Np_min1 = Np - 1, N_min1 = N - 1, p_min1 = p - 1;

  // Check which system matrices are time-varying
  bool Z_tv = Z.n_slices > 1, T_tv = T.n_slices > 1,
       R_tv = R.n_slices > 1, Q_tv = Q.n_slices > 1;

  // Initial system matrices
  arma::mat Z_mat = Z.slice(0), T_mat = T.slice(0),
            R_mat = R.slice(0), Q_mat = Q.slice(0);
  arma::rowvec Z_row = Z_mat.row(0);

  // Indicator for whether the first row should be assigned
  bool row_assign = Z_tv || p > 1;

  // Initialise objects used in computations
  arma::colvec M_inf(m), M_star(m), K_0(m), K_1(m);
  arma::cube L_0(m, m, Np), L_1(m, m, Np);
  Rcpp::NumericVector F_inf(Np), F_star(Np),
                      F_1(Np), v_UT(Np);
  double F_2;

  // Initialise state and corresponding variance
  arma::colvec a_vec;
  arma::mat P_inf_mat, P_star_mat;

  // Initialise r for smoother
  arma::cube r_UT(m, nsim, Np + 1, arma::fill::zeros),
             r_vec(m, nsim, N, arma::fill::zeros);

  // Initialise helpers r_1, QtR, and tT
  arma::colvec r_1(m, arma::fill::zeros);
  arma::mat r_1_mat(m, nsim);
  arma::cube QtR(r, m, N), tT(m, m, N);
  arma::mat QtR_mat = Q_mat * R_mat.t(), tT_mat = T_mat.t();
  QtR.slice(0) = QtR_mat;
  tT.slice(0) = tT_mat;

  // Indicator for whether QtR is time-varying
  bool QtR_tv = Q_tv || R_tv;

  // Initialising smoothed state and its transposed variant, and residuals
  arma::cube a_smooth(N, m, nsim), a_t(m, nsim, N),
             eta(N, r, nsim, arma::fill::zeros);
  arma::mat a_temp(m, nsim), eta_temp(r, nsim);

  // Iterators
  int sim, i, j;

  // Loop over random samples
  for (sim = 0; sim < nsim; sim++) {

    // Set index to the first index
    index = 0;

    // Reset state and corresponding variance, and r_1
    if (sim > 0) {
      a_vec = a;
      P_inf_mat = P_inf;
      P_star_mat = P_star;
      r_1.zeros();
    }

    // Kalman Filter
    // Loop over timepoints
    for (i = 0; i < N; i++) {

      // Get system matrices of current timepoint
      if (Z_tv && i > 0) {
        Z_mat = Z.slice(i);
      }
      if (T_tv && (i > 0 || sim > 0)) {
        T_mat = T.slice(i);
      }
      if (R_tv && (i > 0 || sim > 0)) {
        R_mat = R.slice(i);
      }
      if (Q_tv && (i > 0 || sim > 0)) {
        Q_mat = Q.slice(i);
      }

      // These matrices only need to be computed once
      if (sim == 0 && i > 0) {
        if (QtR_tv) {
          QtR_mat = Q_mat * R_mat.t();
          QtR.slice(i) = QtR_mat;
        }
        if (T_tv) {
          tT_mat = T_mat.t();
          tT.slice(i) = tT_mat;
        }
      } else if (sim > 0) {
        if (QtR_tv) {
          QtR_mat = QtR.slice(i);
        }
        if (T_tv && i > 0) {
          tT_mat = tT.slice(i);
        }
      }

      // Loop over dependent variables
      for (j = 0; j < p; j++, index++) {

        // Retrieve row of Z
        if (j > 0 || (i > 0 && row_assign)) {
          Z_row = Z_mat.row(j);
        }

        // Exact Kalman filter in initialisation steps
        if (index < initialisation_steps) {

          // PZ' as in Kalman formulae
          M_inf = P_inf_mat * Z_row.t();
          M_star = P_star_mat * Z_row.t();

          // Variance matrix of the current residual/fitted value
          F_inf(index) = arma::as_scalar(Z_row * M_inf);
          F_star(index) = arma::as_scalar(Z_row * M_star);

          // Check if F_inf is nearly 0
          if (F_inf(index) < 1e-7) {

            // Check if F_star is nearly 0
            if (F_star(index) < 1e-7) {

              // No new information
              continue;
            } else {

              // Inverse of Fmat
              F_1(index) = 1 / F_star(index);

              // Current residual
              v_UT(index) = y(i, j, sim) - arma::as_scalar(Z_row * a_vec);

              // Auxiliary matrices
              K_0 = M_star * F_1(index);
              L_0.slice(index) = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;

              // Estimated state vector and corresponding variance - covariance
              // matrix for the next step
              if (index < Np_min1) {
                a_vec = a_vec + K_0 * v_UT(index);
                P_star_mat = P_star_mat * L_0.slice(index).t();
              }
              continue;
            }
          } else {

            // Inverse of Fmat
            F_1(index) = 1 / F_inf(index);
            F_2 = -pow(F_1(index), 2.0) * F_star(index);

            // Current residual
            v_UT(index) = y(i, j, sim) - arma::as_scalar(Z_row * a_vec);

            // Auxiliary matrices
            K_0 = M_inf * F_1(index);
            L_0.slice(index) = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;
            K_1 = M_star * F_1(index) + M_inf * F_2;
            L_1.slice(index) = -K_1 * Z_row;

            // Estimated state vector and corresponding variance - covariance
            // matrix for the next step
            if (index < Np_min1) {
              a_vec = a_vec + K_0 * v_UT(index);
              P_star_mat = P_inf_mat * L_1.slice(index).t() +
                P_star_mat * L_0.slice(index).t();
              P_inf_mat = P_inf_mat * L_0.slice(index).t();
            }
          }
        } else {

          // PZ' as in Kalman formulae
          M_star = P_star_mat * Z_row.t();

          // Variance matrix of the current residual/fitted value
          F_star(index) = arma::as_scalar(Z_row * M_star);

          // Check if F_star is nearly 0
          if (F_star(index) < 1e-7) {

            // No new information
            continue;
          } else {

            // Inverse of Fmat
            F_1(index) = 1 / F_star(index);

            // Current residual
            v_UT(index) = y(i, j, sim) - arma::as_scalar(Z_row * a_vec);

            // Auxiliary matrices
            K_0 = M_star * F_1(index);
            L_0.slice(index) = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;

            // Estimated state vector and corresponding variance - covariance
            // matrix for the next step
            if (index < Np_min1) {
              a_vec = a_vec + K_0 * v_UT(index);
              P_star_mat = P_star_mat * L_0.slice(index).t();
            }
          }
        }
      }

      // Transition to the next timepoint
      if (i < N_min1) {
        a_vec = T_mat * a_vec;
        P_star_mat = T_mat * P_star_mat * tT_mat + R_mat * QtR_mat;
        if ((index - 1) < initialisation_steps) {
          P_inf_mat = T_mat * P_inf_mat * tT_mat;
        }
      }
    }

    // Set index to the last index
    index = Np_min1;

    // Kalman Smoother
    // Loop backwards over timepoints
    for (i = N_min1; i >= 0; i--) {

      // Get system matrix of current timepoint for Z
      if (Z_tv && i < N_min1) {
        Z_mat = Z.slice(i);
      }

      // Get system matrix of previous timepoint for T
      if (T_tv && i > 0) {
        tT_mat = tT.slice(i - 1);
      }

      // Loop backwards over dependent variables
      for (j = p_min1; j >= 0; j--, index--) {

        // Retrieve row of Z
        if (j < p_min1 || (i < N_min1 && row_assign)) {
          Z_row = Z_mat.row(j);
        }

        // Exact Kalman smoother in initialisation steps
        if (index < initialisation_steps) {

          // Check if F_inf is nearly 0
          if (F_inf(index) < 1e-7) {

            // Check if F_star is nearly 0
            if (F_star(index) < 1e-7) {

              // No new information
              r_UT.slice(index).col(sim) = r_UT.slice(index + 1).col(sim);
            } else {

              // New r
              r_UT.slice(index).col(sim) = Z_row.t() * F_1(index) *
                v_UT(index) + r_UT.slice(index + 1).col(sim) *
                L_0.slice(index);
            }
          } else {

            // New r
            r_UT.slice(index).col(sim) =
              r_UT.slice(index + 1).col(sim) * L_0.slice(index);
            r_1 = Z_row.t() * F_1(index) * v_UT(index) + r_1 *
              L_0.slice(index) + r_UT.slice(index + 1).col(sim) *
              L_1.slice(index);
          }
        } else {

          // Check if F_star is nearly 0
          if (F_star(index) < 1e-7) {

            // No new information
            r_UT.slice(index).col(sim) = r_UT.slice(index + 1).col(sim);
          } else {

            // New r
            r_UT.slice(index).col(sim) = Z_row.t() * F_1(index) *
              v_UT(index) + r_UT.slice(index + 1).col(sim) *
              L_0.slice(index);
          }
        }
      }

      // Save r for each timepoint
      r_vec.slice(i).col(sim) = r_UT.slice(index + 1).col(sim);
      r_1_mat.col(sim) = r_1;

      // r and N for the previous timepoint, not valid for i = 0
      if (i > 0) {
        r_UT.slice(index + 1).col(sim) =
          tT_mat * r_UT.slice(index + 1).col(sim);
        if ((index + 1) < initialisation_steps) {
          r_1 = tT_mat * r_1;
        }
      }
    }
  }

  // Initial smoothed state
  a_temp = arma::repmat(a, 1, nsim) + P_star * r_vec.slice(0) + P_inf * r_1_mat;
  a_smooth.row(0) = a_temp;
  if (transposed_state) {
    a_t.slice(0) = a_temp;
  }

  // Calculate smoothed eta and state
  for (i = 1; i < N; i++) {

    // Get system matrices of current timepoint
    if (T_tv) {
      T_mat = T.slice(i - 1);
    }
    if (R_tv) {
      R_mat = R.slice(i - 1);
    }
    if (QtR_tv) {
      QtR_mat = QtR.slice(i - 1);
    }

    // Calculate smoothed state disturbance
    eta_temp = QtR_mat * r_vec.slice(i);
    eta.row(i - 1) = eta_temp;

    // Calculate smoothed state
    a_temp = T_mat * a_temp + R_mat * eta_temp;
    a_smooth.row(i) = a_temp;
    if (transposed_state) {
      a_t.slice(i) = a_temp;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("a_smooth") = a_smooth,
    Rcpp::Named("a_t") = a_t,
    Rcpp::Named("eta") = eta
  );
}
