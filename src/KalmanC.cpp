#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]

//' Employ the Kalman Filter and Smoother
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
//' @param diagnostics Boolean indicating whether diagnostics should be computed.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List KalmanC(const arma::mat& y,
                   const Rcpp::LogicalMatrix& y_isna,
                   const arma::colvec& a,
                   const arma::mat& P_inf,
                   const arma::mat& P_star,
                   const arma::cube& Z,
                   const arma::cube& T,
                   const arma::cube& R,
                   const arma::cube& Q,
                   const bool& diagnostics) {

  // Number of observations, dependent variables, state parameters,
  // and state disturbances
  int N = y.n_rows, p = y.n_cols, m = a.n_rows, r = R.n_cols;

  // Keep track of the current index and limits of indices
  int index = 0, Np = N * p, Np_min1 = Np - 1, N_min1 = N - 1, p_min1 = p - 1;

  // Check which system matrices are time-varying
  bool Z_tv = Z.n_slices > 1, T_tv = T.n_slices > 1,
       R_tv = R.n_slices > 1, Q_tv = Q.n_slices > 1;

  // Initial system matrices
  arma::mat Z_mat = Z.slice(0), T_mat = T.slice(0),
            R_mat = R.slice(0), Q_mat = Q.slice(0);
  arma::rowvec Z_row = Z_mat.row(0);

  // Check if P_inf is already 0
  bool initialisation = !arma::all(arma::vectorise(arma::abs(P_inf)) < 1e-7);

  // Initialise number of initialisation steps
  int initialisation_steps = 0;

  // Initialise objects used in computations
  arma::mat M_inf(m, Np), M_star(m, Np), K_0(m, Np), K_1(m, Np);
  arma::cube L_0(m, m, Np), L_1(m, m, Np);
  Rcpp::NumericVector F_inf(Np), F_star(Np),
                      F_1(Np), F_2(Np), v_UT(Np);
  double constant = -log(2 * M_PI) / 2;
  double kappa = 1.0e7;

  // Initialise loglikelihood
  Rcpp::NumericVector loglik(Np);

  // Initialise (filtered and smoothed) state and corresponding variance
  arma::mat a_mat(Np, m), a_pred(N, m), a_fil(N, m), a_smooth(N, m);
  arma::cube P_inf_cube(m, m, Np, arma::fill::zeros),
             P_star_cube(m, m, Np), P_pred(m, m, N),
             P_fil(m, m, N), V(m, m, N);
  a_mat.row(0) = a.t();
  P_inf_cube.slice(0) = P_inf;
  P_star_cube.slice(0) = P_star;

  // Initialising fitted values, (smoothed) residuals and variance
  arma::mat yfit(N, p), v(N, p), v_norm(N, p), eta(N, r), e(N, p);
  arma::cube Fmat(p, p, N), eta_var(r, r, N), D(p, p, N);
  v_norm.fill(NA_REAL);

  // Helpers for calculation of v_norm
  bool calculate_vnorm = false;
  arma::cube Finv(p, p, N);
  arma::mat U_Finv(p, p), V_Finv(p, p), Finv_root(p, p), v_0na(N, p);
  arma::colvec s_Finv(p);
  arma::rowvec v_norm_row(p);

  // One step ahead forecast of the state and corresponding uncertainty
  arma::colvec a_fc(m);
  arma::mat P_inf_fc(m, m, arma::fill::zeros), P_star_fc(m, m), P_fc(m, m);

  // Initialise Tstats
  arma::mat Tstat_observation(N, p), Tstat_state(N + 1, m);

  // Inititalise r and N for smoother
  arma::mat r_UT(Np + 1, m, arma::fill::zeros),
            r_vec(N + 1, m, arma::fill::zeros);
  arma::cube N_UT(m, m, Np + 1, arma::fill::zeros),
             Nmat(m, m, N + 1, arma::fill::zeros);

  // Initialise helpers r_1, N_1, N_2, tZZ, K, e_row, RtQ for smoother
  arma::rowvec r_1(m, arma::fill::zeros), e_row(p);
  arma::mat N_1(m, m, arma::fill::zeros), N_2(m, m, arma::fill::zeros),
            tZZ(m, m), K(m, p), RtQ(m, r);

  // Kalman Filter
  // Loop over timepoints
  for (int i = 0; i < N; i++) {

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

    // Predicted state and corresponding uncertainty
    a_pred.row(i) = a_mat.row(index);
    if (initialisation) {
      P_pred.slice(i) =
        kappa * P_inf_cube.slice(index) + P_star_cube.slice(index);
      if (diagnostics) {
        calculate_vnorm = arma::all(arma::vectorise(arma::abs(
          Z_mat * P_inf_cube.slice(index) * Z_mat.t())) < 1e-7);
      }
    } else {
      P_pred.slice(i) = P_star_cube.slice(index);
    }

    // Predicted values
    yfit.row(i) = arma::trans(Z_mat * a_pred.row(i).t());

    // Residuals and corresponding uncertainty
    v.row(i) = y.row(i) - yfit.row(i);
    Fmat.slice(i) = Z_mat * P_pred.slice(i) * Z_mat.t();

    // Normalised residuals
    if (diagnostics && (!initialisation || calculate_vnorm)) {

      // Inverse of Fmat
      Finv.slice(i) = arma::inv_sympd(Fmat.slice(i));

      // Singular Value decomposition of the inverse of Fmat
      arma::svd(U_Finv, s_Finv, V_Finv, Finv.slice(i));

      // Symmetric square root of Finv
      Finv_root = U_Finv * arma::diagmat(arma::sqrt(s_Finv)) * U_Finv.t();

      // Normalised prediction error
      v_0na.row(i) = v.row(i);
      v_0na.row(i).replace(NA_REAL, 0);
      v_norm_row = v_0na.row(i) * Finv_root;
      v_norm_row.elem(arma::find_nonfinite(v.row(i))).fill(NA_REAL);
      v_norm.row(i) = v_norm_row;
    }

    // Loop over dependent variables
    for (int j = 0; j < p; j++, index++) {

      // Check for missing value
      if (y_isna(i, j)) {
        loglik(index) = NA_REAL;
        if (index < Np_min1) {
          a_mat.row(index + 1) = a_mat.row(index);
          P_star_cube.slice(index + 1) = P_star_cube.slice(index);
          if (initialisation) {
            P_inf_cube.slice(index + 1) = P_inf_cube.slice(index);
          }
        } else {
          a_fc = a_mat.row(index).t();
          P_star_fc = P_star_cube.slice(index);
          if (initialisation) {
            P_inf_fc = P_inf_cube.slice(index);
          }
        }
        continue;
      }

      // Retrieve row of Z
      if (i > 0 || j > 0) {
        Z_row = Z_mat.row(j);
      }

      // Exact Kalman filter in initialisation steps
      if (initialisation) {

        // Increment initialisation_steps
        initialisation_steps++;

        // PZ' as in Kalman formulae
        M_inf.col(index) = P_inf_cube.slice(index) * Z_row.t();
        M_star.col(index) = P_star_cube.slice(index) * Z_row.t();

        // Variance matrix of the current residual/fitted value
        F_inf(index) = arma::as_scalar(Z_row * M_inf.col(index));
        F_star(index) = arma::as_scalar(Z_row * M_star.col(index));

        // Check if F_inf is nearly 0
        if (F_inf(index) < 1e-7) {

          // Check if F_star is nearly 0
          if (F_star(index) < 1e-7) {

            // No new information
            loglik(index) = NA_REAL;
            if (index < Np_min1) {
              a_mat.row(index + 1) = a_mat.row(index);
              P_star_cube.slice(index + 1) = P_star_cube.slice(index);
              P_inf_cube.slice(index + 1) = P_inf_cube.slice(index);
            } else {
              a_fc = a_mat.row(index).t();
              P_star_fc = P_star_cube.slice(index);
              P_inf_fc = P_inf_cube.slice(index);
            }
            continue;
          } else {

            // Inverse of Fmat
            F_1(index) = 1 / F_star(index);

            // Current residual
            v_UT(index) = y(i, j) - arma::as_scalar(Z_row * a_mat.row(index).t());

            // Auxiliary matrices
            K_0.col(index) = M_star.col(index) * F_1(index);
            L_0.slice(index) =
              arma::mat(m, m, arma::fill::eye) - K_0.col(index) * Z_row;

            // Estimated state vector and corresponding variance - covariance
            // matrix for the next step
            if (index < Np_min1) {
              a_mat.row(index + 1) =
                a_mat.row(index) + K_0.col(index).t() * v_UT(index);
              P_star_cube.slice(index + 1) =
                P_star_cube.slice(index) * L_0.slice(index).t();
              P_inf_cube.slice(index + 1) = P_inf_cube.slice(index);
            } else {
              a_fc = a_mat.row(index).t() + K_0.col(index) * v_UT(index);
              P_star_fc = P_star_cube.slice(index) * L_0.slice(index).t();
              P_inf_fc = P_inf_cube.slice(index);
            }
            loglik(index) = constant -
              (log(F_star(index)) + pow(v_UT(index), 2.0) * F_1(index)) / 2;
            continue;
          }
        } else {

          // Inverse of Fmat
          F_1(index) = 1 / F_inf(index);
          F_2(index) = -pow(F_1(index), 2.0) * F_star(index);

          // Current residual
          v_UT(index) = y(i, j) - arma::as_scalar(Z_row * a_mat.row(index).t());

          // Auxiliary matrices
          K_0.col(index) = M_inf.col(index) * F_1(index);
          L_0.slice(index) =
            arma::mat(m, m, arma::fill::eye) - K_0.col(index) * Z_row;
          K_1.col(index) =
            M_star.col(index) * F_1(index) + M_inf.col(index) * F_2(index);
          L_1.slice(index) = -K_1.col(index) * Z_row;

          // Estimated state vector and corresponding variance - covariance
          // matrix for the next step
          if (index < Np_min1) {
            a_mat.row(index + 1) =
              a_mat.row(index) + K_0.col(index).t() * v_UT(index);
            P_star_cube.slice(index + 1) = P_inf_cube.slice(index) *
              L_1.slice(index).t() +
              P_star_cube.slice(index) * L_0.slice(index).t();
            P_inf_cube.slice(index + 1) =
              P_inf_cube.slice(index) * L_0.slice(index).t();
          } else {
            a_fc = a_mat.row(index).t() + K_0.col(index) * v_UT(index);
            P_star_fc = P_inf_cube.slice(index) *
              L_1.slice(index).t() +
              P_star_cube.slice(index) * L_0.slice(index).t();
            P_inf_fc = P_inf_cube.slice(index) * L_0.slice(index).t();
          }
          loglik(index) = constant - log(F_inf(index)) / 2;
        }

        // Check if P_inf converged to 0
        if (index < Np_min1 && j < p_min1) {
          initialisation = !arma::all(
            arma::vectorise(arma::abs(P_inf_cube.slice(index + 1))) < 1e-7
          );
        }
      } else {

        // PZ' as in Kalman formulae
        M_star.col(index) = P_star_cube.slice(index) * Z_row.t();

        // Variance matrix of the current residual/fitted value
        F_star(index) = arma::as_scalar(Z_row * M_star.col(index));

        // Check if F_star is nearly 0
        if (F_star(index) < 1e-7) {

          // No new information
          loglik(index) = NA_REAL;
          if (index < Np_min1) {
            a_mat.row(index + 1) = a_mat.row(index);
            P_star_cube.slice(index + 1) = P_star_cube.slice(index);
          } else {
            a_fc = a_mat.row(index).t();
            P_star_fc = P_star_cube.slice(index);
          }
        } else {

          // Inverse of Fmat
          F_1(index) = 1 / F_star(index);

          // Current residual
          v_UT(index) = y(i, j) - arma::as_scalar(Z_row * a_mat.row(index).t());

          // Auxiliary matrices
          K_0.col(index) = M_star.col(index) * F_1(index);
          L_0.slice(index) =
            arma::mat(m, m, arma::fill::eye) - K_0.col(index) * Z_row;

          // Estimated state vector and corresponding variance - covariance
          // matrix for the next step
          if (index < Np_min1) {
            a_mat.row(index + 1) =
              a_mat.row(index) + K_0.col(index).t() * v_UT(index);
            P_star_cube.slice(index + 1) =
              P_star_cube.slice(index) * L_0.slice(index).t();
          } else {
            a_fc = a_mat.row(index).t() + K_0.col(index) * v_UT(index);
            P_star_fc = P_star_cube.slice(index) * L_0.slice(index).t();
          }
          loglik(index) = constant -
            (log(F_star(index)) + pow(v_UT(index), 2.0) * F_1(index)) / 2;
        }
      }
    }

    // Filtered state and corresponding uncertainty
    // Followed by transition to the next timepoint
    if (index < Np) {
      a_fil.row(i) = a_mat.row(index);
      if (initialisation) {
        P_fil.slice(i) =
          kappa * P_inf_cube.slice(index) + P_star_cube.slice(index);
      } else {
        P_fil.slice(i) = P_star_cube.slice(index);
      }
      a_mat.row(index) = arma::trans(T_mat * a_mat.row(index).t());
      P_star_cube.slice(index) =
        T_mat * P_star_cube.slice(index) * T_mat.t() + R_mat * Q_mat * R_mat.t();
      if (initialisation) {
        P_inf_cube.slice(index) = T_mat * P_inf_cube.slice(index) * T_mat.t();
        // Check if P_inf converged to 0
        initialisation = !arma::all(
          arma::vectorise(arma::abs(P_inf_cube.slice(index))) < 1e-7
        );
      }
    } else {
      a_fil.row(i) = a_fc.t();
      if (initialisation) {
        P_fil.slice(i) = kappa * P_inf_fc + P_star_fc;
      } else {
        P_fil.slice(i) = P_star_fc;
      }
      a_fc = T_mat * a_fc;
      P_star_fc = T_mat * P_star_fc * T_mat.t() + R_mat * Q_mat * R_mat.t();
      if (initialisation) {
        P_inf_fc = T_mat * P_inf_fc * T_mat.t();
        P_fc = kappa * P_inf_fc + P_star_fc;
      } else {
        P_fc = P_star_fc;
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

    // Get system matrices of previous timepoint for T, R, and Q
    if (T_tv && i > 0) {
      T_mat = T.slice(i - 1);
    }
    if (R_tv && i > 0) {
      R_mat = R.slice(i - 1);
    }
    if (Q_tv && i > 0) {
      Q_mat = Q.slice(i - 1);
    }

    // Loop backwards over dependent variables
    for (int j = p_min1; j >= 0; j--, index--) {

      // Check for missing value
      if (y_isna(i, j)) {
        r_UT.row(index) = r_UT.row(index + 1);
        N_UT.slice(index) = N_UT.slice(index + 1);
        continue;
      }

      // Retrieve row of Z
      if (i < N_min1 || j < p_min1) {
        Z_row = Z_mat.row(j);
      }

      // Exact Kalman smoother in initialisation steps
      if (index < initialisation_steps) {

        // Check if F_inf is nearly 0
        if (F_inf(index) < 1e-7) {

          // Check if F_star is nearly 0
          if (F_star(index) < 1e-7) {

            // No new information
            r_UT.row(index) = r_UT.row(index + 1);
            N_UT.slice(index) = N_UT.slice(index + 1);
          } else {

            // New r and N
            r_UT.row(index) = Z_row * F_1(index) * v_UT(index) +
              r_UT.row(index + 1) * L_0.slice(index);
            N_UT.slice(index) = Z_row.t() * Z_row * F_1(index) +
              L_0.slice(index).t() * N_UT.slice(index + 1) * L_0.slice(index);
            N_1 = N_1 * L_0.slice(index);
          }
        } else {

          // New r and N
          r_UT.row(index) = r_UT.row(index + 1) * L_0.slice(index);
          r_1 = Z_row * F_1(index) * v_UT(index) +
            r_1 * L_0.slice(index) + r_UT.row(index + 1) * L_1.slice(index);
          N_UT.slice(index) =
            L_0.slice(index).t() * N_UT.slice(index + 1) * L_0.slice(index);
          tZZ = Z_row.t() * Z_row;
          N_1 = tZZ * F_1(index) +
            L_0.slice(index).t() * N_1 * L_0.slice(index) +
            L_1.slice(index).t() * N_UT.slice(index + 1) * L_0.slice(index) +
            L_0.slice(index).t() * N_UT.slice(index + 1) * L_1.slice(index);
          N_2 = tZZ * F_2(index) +
            L_0.slice(index).t() * N_2 * L_0.slice(index) +
            L_1.slice(index).t() * N_1.t() * L_0.slice(index) +
            L_1.slice(index).t() * N_UT.slice(index + 1) * L_1.slice(index);
        }
      } else {

        // Check if F_star is nearly 0
        if (F_star(index) < 1e-7) {

          // No new information
          r_UT.row(index) = r_UT.row(index + 1);
          N_UT.slice(index) = N_UT.slice(index + 1);
        } else {

          // New r and N
          r_UT.row(index) = Z_row * F_1(index) * v_UT(index) +
            r_UT.row(index + 1) * L_0.slice(index);
          N_UT.slice(index) = Z_row.t() * Z_row * F_1(index) +
            L_0.slice(index).t() * N_UT.slice(index + 1) * L_0.slice(index);
        }
      }
    }

    // Save r and N for each timepoint
    r_vec.row(i) = r_UT.row(index + 1);
    Nmat.slice(i) = N_UT.slice(index + 1);

    // Compute smoothed state and variance
    if ((index + 1) < initialisation_steps) {
      a_smooth.row(i) = a_pred.row(i) +
        r_vec.row(i) * P_star_cube.slice(index + 1) +
        r_1 * P_inf_cube.slice(index + 1);
      V.slice(i) = P_star_cube.slice(index + 1) -
        P_star_cube.slice(index + 1) * Nmat.slice(i) * P_star_cube.slice(index + 1) -
        P_star_cube.slice(index + 1) * N_1.t() * P_inf_cube.slice(index + 1) -
        P_inf_cube.slice(index + 1) * N_1 * P_star_cube.slice(index + 1) -
        P_inf_cube.slice(index + 1) * N_2 * P_inf_cube.slice(index + 1);
    } else {
      a_smooth.row(i) = a_pred.row(i) +
        r_vec.row(i) * P_pred.slice(i);
      V.slice(i) = P_pred.slice(i) -
        P_pred.slice(i) * Nmat.slice(i) * P_pred.slice(i);
    }

    if (diagnostics) {

      // Tstat for the state equation
      Tstat_state.row(i + 1) = r_vec.row(i + 1) /
        arma::sqrt(arma::diagvec(Nmat.slice(i + 1)).t());

      // Smoothing error e and corresponding variance D
      K = T_mat * P_pred.slice(i) * Z_mat.t() * Finv.slice(i);
      e_row = v_0na.row(i) * Finv.slice(i) - r_vec.row(i + 1) * K;
      e_row.elem(arma::find_nonfinite(v.row(i))).fill(NA_REAL);
      e.row(i) = e_row;
      D.slice(i) = Finv.slice(i) + K.t() * Nmat.slice(i + 1) * K;

      // Tstat for the observation equation
      Tstat_observation.row(i + 1) = e_row /
        arma::sqrt(arma::diagvec(D.slice(i)).t());
    }

    // Compute smoothed state disturbance and corresponding variance
    RtQ = R_mat * Q_mat.t();
    eta.row(i) = r_vec.row(i + 1) * RtQ;
    eta_var.slice(i) = Q_mat - RtQ.t() * Nmat.slice(i + 1) * RtQ;

    // r and N for the previous timepoint, not valid for i = 0
    if (i > 0) {
      r_UT.row(index + 1) = r_UT.row(index + 1) * T_mat;
      N_UT.slice(index + 1) = T_mat.t() * N_UT.slice(index + 1) * T_mat;
      if ((index + 1) < initialisation_steps) {
        r_1 = r_1 * T_mat;
        N_1 = T_mat.t() * N_1 * T_mat;
        N_2 = T_mat.t() * N_2 * T_mat;
      }
    }
  }

  // List to return
  return Rcpp::List::create(
    Rcpp::Named("initialisation_steps") = initialisation_steps,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("a_pred") = a_pred,
    Rcpp::Named("a_fil") = a_fil,
    Rcpp::Named("a_smooth") = a_smooth,
    Rcpp::Named("P_pred") = P_pred,
    Rcpp::Named("P_fil") = P_fil,
    Rcpp::Named("V") = V,
    Rcpp::Named("yfit") = yfit,
    Rcpp::Named("v") = v,
    Rcpp::Named("v_norm") = v_norm,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("e") = e,
    Rcpp::Named("Fmat") = Fmat,
    Rcpp::Named("eta_var") = eta_var,
    Rcpp::Named("D") = D,
    Rcpp::Named("a_fc") = a_fc,
    Rcpp::Named("P_fc") = P_fc,
    Rcpp::Named("Tstat_observation") = Tstat_observation,
    Rcpp::Named("Tstat_state") = Tstat_state,
    Rcpp::Named("r_vec") = r_vec,
    Rcpp::Named("Nmat") = Nmat
  );
}
