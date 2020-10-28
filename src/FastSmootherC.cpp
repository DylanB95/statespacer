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
                         const int& initialisation_steps) {

  // Number of observations, dependent variables, state parameters,
  // and state disturbances
  int N = y.n_rows, p = y.n_cols, m = a.n_rows, r = R.n_cols, nsim = y.n_slices;

  // Keep track of the current index and limits of indices
  int index = 0, Np = N * p, Np_min1 = Np - 1, N_min1 = N - 1, p_min1 = p - 1;

  // Check which system matrices are time-varying
  bool Z_tv = Z.n_slices > 1, T_tv = T.n_slices > 1,
       R_tv = R.n_slices > 1, Q_tv = Q.n_slices > 1;

  // Indicator for whether the first row should be assigned
  bool row_assign = Z_tv || p > 1;

  // Initial system matrices
  arma::mat Z_mat = Z.slice(0), T_mat = T.slice(0),
            R_mat = R.slice(0), Q_mat = Q.slice(0);
  arma::rowvec Z_row = Z_mat.row(0);

  // Initialise objects used in computations
  arma::colvec M_inf(m), M_star(m), K_0(m), K_1(m);
  arma::cube L_0(m, m, Np), L_1(m, m, Np);
  Rcpp::NumericMatrix F_inf(Np, nsim), F_star(Np, nsim),
                      F_1(Np, nsim), v_UT(Np, nsim);
  double F_2;

  // Initialise state and corresponding variance
  arma::colvec a_vec;
  arma::mat P_inf_mat, P_star_mat;

  // Initialising smoothed state and residuals
  arma::cube a_smooth(N, m, nsim), eta(N, r, nsim);
  arma::mat a_temp(m, nsim), eta_temp(r, nsim);

  // Initialise r for smoother
  arma::cube r_UT(m, nsim, Np + 1, arma::fill::zeros),
             r_vec(m, nsim, N + 1, arma::fill::zeros);

  // Initialise helpers r_1, QtR, and tT
  arma::colvec r_1(m, arma::fill::zeros);
  arma::mat r_1_vec(m, nsim);
  arma::cube QtR(r, m, N), tT(m, m, N);
  arma::mat QtR_mat = Q_mat * R_mat.t(), tT_mat = T_mat.t();
  QtR.slice(0) = QtR_mat;
  tT.slice(0) = tT_mat;

  // Iterators
  int sim, i, j;

  // Loop over random samples
  for (sim = 0; sim < nsim; sim++) {

    // Set index to the first index
    index = 0;

    // Reset state and corresponding variance
    a_vec = a;
    P_inf_mat = P_inf;
    P_star_mat = P_star;

    // Kalman Filter
    // Loop over timepoints
    for (i = 0; i < N; i++) {

      // Get system matrices of current timepoint
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

      // These matrices only need to be computed once
      if (sim == 0 && i > 0) {
        if (Q_tv || R_tv) {
          QtR_mat = Q_mat * R_mat.t();
          QtR.slice(i) = QtR_mat;
        }
        if (T_tv) {
          tT_mat = T_mat.t();
          tT.slice(i) = tT_mat;
        }
      }
      if (sim > 0 && i > 0) {
        if (Q_tv || R_tv) {
          QtR_mat = QtR.slice(i);
        }
        if (T_tv) {
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
          F_inf(index, sim) = arma::as_scalar(Z_row * M_inf);
          F_star(index, sim) = arma::as_scalar(Z_row * M_star);

          // Check if F_inf is nearly 0
          if (F_inf(index, sim) < 1e-7) {

            // Check if F_star is nearly 0
            if (F_star(index, sim) < 1e-7) {

              // No new information
              continue;
            } else {

              // Inverse of Fmat
              F_1(index, sim) = 1 / F_star(index, sim);

              // Current residual
              v_UT(index, sim) = y(i, j, sim) - arma::as_scalar(Z_row * a_vec);

              // Auxiliary matrices
              K_0 = M_star * F_1(index, sim);
              L_0.slice(index) = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;

              // Estimated state vector and corresponding variance - covariance
              // matrix for the next step
              if (index < Np_min1) {
                a_vec = a_vec + K_0 * v_UT(index, sim);
                P_star_mat = P_star_mat * L_0.slice(index).t();
              }
              continue;
            }
          } else {

            // Inverse of Fmat
            F_1(index, sim) = 1 / F_inf(index, sim);
            F_2 = -pow(F_1(index, sim), 2.0) * F_star(index, sim);

            // Current residual
            v_UT(index, sim) = y(i, j, sim) - arma::as_scalar(Z_row * a_vec);

            // Auxiliary matrices
            K_0 = M_inf * F_1(index, sim);
            L_0.slice(index) = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;
            K_1 = M_star * F_1(index, sim) + M_inf * F_2;
            L_1.slice(index) = -K_1 * Z_row;

            // Estimated state vector and corresponding variance - covariance
            // matrix for the next step
            if (index < Np_min1) {
              a_vec = a_vec + K_0 * v_UT(index, sim);
              P_star_mat = P_inf_mat * L_1.slice(index).t() +
                P_star_mat * L_0.slice(index).t();
              P_inf_mat = P_inf_mat * L_0.slice(index).t();
            }
          }
        } else {

          // PZ' as in Kalman formulae
          M_star = P_star_mat * Z_row.t();

          // Variance matrix of the current residual/fitted value
          F_star(index, sim) = arma::as_scalar(Z_row * M_star);

          // Check if F_star is nearly 0
          if (F_star(index, sim) < 1e-7) {

            // No new information
            continue;
          } else {

            // Inverse of Fmat
            F_1(index, sim) = 1 / F_star(index, sim);

            // Current residual
            v_UT(index, sim) = y(i, j, sim) - arma::as_scalar(Z_row * a_vec);

            // Auxiliary matrices
            K_0 = M_star * F_1(index, sim);
            L_0.slice(index) = arma::mat(m, m, arma::fill::eye) - K_0 * Z_row;

            // Estimated state vector and corresponding variance - covariance
            // matrix for the next step
            if (index < Np_min1) {
              a_vec = a_vec + K_0 * v_UT(index, sim);
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
    for (int i = N_min1; i >= 0; i--) {

      // Get system matrix of current timepoint for Z
      if (Z_tv && i < N_min1) {
        Z_mat = Z.slice(i);
      }
      if ((Q_tv || R_tv) && i < N_min1) {
        QtR_mat = QtR.slice(i);
      }
      if (T_tv && i < N_min1) {
        tT_mat = tT.slice(i);
      }

      // Loop backwards over dependent variables
      for (int j = p_min1; j >= 0; j--, index--) {

        // Retrieve row of Z
        if (j < p_min1 || (i < N_min1 && row_assign)) {
          Z_row = Z_mat.row(j);
        }

        // Exact Kalman smoother in initialisation steps
        if (index < initialisation_steps) {

          // Check if F_inf is nearly 0
          if (F_inf(index, sim) < 1e-7) {

            // Check if F_star is nearly 0
            if (F_star(index, sim) < 1e-7) {

              // No new information
              r_UT.slice(index).col(sim) = r_UT.slice(index + 1).col(sim);
            } else {

              // New r
              r_UT.slice(index).col(sim) = Z_row.t() * F_1(index, sim) *
                v_UT(index, sim) + r_UT.slice(index + 1).col(sim) *
                L_0.slice(index);
            }
          } else {

            // New r
            r_UT.slice(index).col(sim) =
              r_UT.slice(index + 1).col(sim) * L_0.slice(index);
            r_1 = Z_row.t() * F_1(index, sim) * v_UT(index, sim) + r_1 *
              L_0.slice(index) + r_UT.slice(index + 1).col(sim) *
              L_1.slice(index);
          }
        } else {

          // Check if F_star is nearly 0
          if (F_star(index, sim) < 1e-7) {

            // No new information
            r_UT.slice(index).col(sim) = r_UT.slice(index + 1).col(sim);
          } else {

            // New r
            r_UT.slice(index).col(sim) = Z_row.t() * F_1(index, sim) *
              v_UT(index, sim) + r_UT.slice(index + 1).col(sim) *
              L_0.slice(index);
          }
        }
      }

      // Save r for each timepoint
      r_vec.slice(i).col(sim) = r_UT.slice(index + 1).col(sim);
      r_1_vec.col(sim) = r_1;

      // r and N for the previous timepoint, not valid for i = 0
      if (i > 0) {
        r_UT.slice(index + 1).col(sim) = tT.slice(i - 1) * r_UT.slice(index + 1).col(sim);
        if ((index + 1) < initialisation_steps) {
          r_1 = tT.slice(i - 1) * r_1;
        }
      }
    }
  }

  // Calculate smoothed eta and state
  for (int i = 0; i < N; i++) {

  }

  // List to return
  nested["initialisation_steps"] = initialisation_steps;
  nested["loglik"] = loglik;
  nested["a_pred"] = a_pred;
  nested["a_fil"] = a_fil;
  nested["a_smooth"] = a_smooth;
  nested["P_pred"] = P_pred;
  nested["P_fil"] = P_fil;
  nested["V"] = V;
  nested["P_inf_pred"] = P_inf_pred;
  return Rcpp::List::create(
    Rcpp::Named("nested") = nested,
    Rcpp::Named("P_star_pred") = P_star_pred,
    Rcpp::Named("P_inf_fil") = P_inf_fil,
    Rcpp::Named("P_star_fil") = P_star_fil,
    Rcpp::Named("yfit") = yfit,
    Rcpp::Named("v") = v,
    Rcpp::Named("v_norm") = v_norm,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("e") = e,
    Rcpp::Named("Fmat") = Fmat,
    Rcpp::Named("eta_var") = eta_var,
    Rcpp::Named("D") = D,
    Rcpp::Named("a_fc") = a_fc,
    Rcpp::Named("P_inf_fc") = P_inf_fc,
    Rcpp::Named("P_star_fc") = P_star_fc,
    Rcpp::Named("P_fc") = P_fc,
    Rcpp::Named("Tstat_observation") = Tstat_observation,
    Rcpp::Named("Tstat_state") = Tstat_state,
    Rcpp::Named("r_vec") = r_vec,
    Rcpp::Named("Nmat") = Nmat
  );
}
