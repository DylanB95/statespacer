#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]

//' Draw Random Samples of a State Space Component
//'
//' @param nsim Number of random samples to draw.
//' @param repeat_Q Number of times the drawing of random samples
//'   using Q should be repeated.
//' @param N Number of timepoints.
//' @param a Initial values of the state vector of the component.
//' @param Z Z system matrix of the State Space model component.
//' @param T T system matrix of the State Space model component.
//' @param R R system matrix of the State Space model component.
//' @param Q Q system matrix of the State Space model component.
//' @param P_star P_star system matrix of the State Space model component.
//'   Only used when draw_initial = true.
//' @param draw_initial Boolean indicating whether the initial state
//'   should be simulated.
//' @param eta_only Boolean indicating whether only the
//'   disturbances should be simulated.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List SimulateC(const int& nsim,
                     const int& repeat_Q,
                     const int& N,
                     const arma::colvec& a,
                     const arma::cube& Z,
                     const arma::cube& T,
                     const arma::cube& R,
                     const arma::cube& Q,
                     const arma::mat& P_star,
                     const bool& draw_initial,
                     const bool& eta_only) {

  // Number of dependent variables, state parameters and state disturbances
  int p = Z.n_rows, m = a.n_rows, r = Q.n_rows, r_tot = r * repeat_Q,
      nsim_rep = nsim * repeat_Q, total = r_tot * nsim, m_nsim = m * nsim;

  // The last time point
  int N_min1 = N - 1;

  // Check which system matrices are time-varying
  bool Z_tv = Z.n_slices > 1, T_tv = T.n_slices > 1,
       R_tv = R.n_slices > 1, Q_tv = Q.n_slices > 1;

  // Check if repeat_Q is greater than 1
  bool rep_g1 = repeat_Q > 1;

  // Initial system matrices
  arma::mat Z_mat = Z.slice(0), T_mat = T.slice(0),
            R_mat = R.slice(0), Q_mat = Q.slice(0);

  // Initialise simulated disturbances, state, and dependent
  arma::cube eta(N, r_tot, nsim), a_cube(N, m, nsim), y(N, p, nsim);
  arma::mat eta_sim(r, nsim_rep), eta_temp(r_tot, nsim);

  // Initialise vector of random draws
  Rcpp::NumericVector draw(total), a_draw(m_nsim);

  // Helpers for calculating root of Q and P_star
  arma::mat Q_root(r, r), U_Q(r, r), V_Q(r, r),
            P_root(m, m), U_P(m, m), V_P(m, m);
  arma::colvec s_Q(r), s_P(m);

  // Initial state
  arma::mat a_temp = arma::repmat(a, 1, nsim);
  if (draw_initial) {
    arma::svd(U_P, s_P, V_P, P_star);
    P_root = U_P * arma::diagmat(arma::sqrt(s_P)) * U_P.t();
    a_draw = Rcpp::rnorm(m_nsim);
    a_temp = a_temp + P_root * arma::mat(a_draw.begin(), m, nsim);
  }
  a_cube.row(0) = a_temp;

  // Initial root of Q
  arma::svd(U_Q, s_Q, V_Q, Q_mat);
  Q_root = U_Q * arma::diagmat(arma::sqrt(s_Q)) * U_Q.t();

  // Loop over timepoints
  for (int i = 0; i < N; i++) {

    // Get system matrices of current timepoint
    if (Z_tv && !eta_only && i > 0) {
      Z_mat = Z.slice(i);
    }
    if (T_tv && !eta_only && i > 0) {
      T_mat = T.slice(i);
    }
    if (R_tv && !eta_only && i > 0) {
      R_mat = R.slice(i);
    }
    if (Q_tv && i > 0) {
      Q_mat = Q.slice(i);
      arma::svd(U_Q, s_Q, V_Q, Q_mat);
      Q_root = U_Q * arma::diagmat(arma::sqrt(s_Q)) * U_Q.t();
    }

    // Simulated y component
    if (!eta_only) {
      y.row(i) = Z_mat * a_temp;
    }

    // Draw random numbers following the standard normal distribution
    draw = Rcpp::rnorm(total);

    // Transform such that the random numbers have variance Q
    if (rep_g1) {
      eta_sim = Q_root * arma::mat(draw.begin(), r, nsim_rep);
      eta_temp = arma::mat(eta_sim.begin(), r_tot, nsim);
    } else {
      eta_temp = Q_root * arma::mat(draw.begin(), r_tot, nsim);
    }
    eta.row(i) = eta_temp;

    // Simulated state
    if (!eta_only && i < N_min1) {
      a_temp = T_mat * a_temp + R_mat * eta_temp;
      a_cube.row(i + 1) = a_temp;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("a") = a_cube,
    Rcpp::Named("eta") = eta
  );
}
