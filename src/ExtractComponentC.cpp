#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]

//' Extract Component of Random Samples
//'
//' @param a Array containing the state vectors for each time point
//'   and random sample.
//' @param Z Z system matrix of the State Space model component.
//'
//' @noRd
// [[Rcpp::export]]
arma::cube SimulateC(const arma::cube& a,
                     const arma::cube& Z) {

  // Number of time points, dependent variables, and random samples
  int N = a.n_slices, p = Z.n_rows, nsim = a.n_cols;

  // Check whether Z is time-varying
  bool Z_tv = Z.n_slices > 1;

  // Initial Z system matrix
  arma::mat Z_mat = Z.slice(0);

  // Initialise component
  arma::cube component(N, p, nsim);

  // Loop over time points
  for (int i = 0; i < N; i++) {

    // Get system matrices of current time point
    if (Z_tv && i > 0) {
      Z_mat = Z.slice(i);
    }

    // Simulated component
    component.row(i) = Z_mat * a.slice(i);
  }

  return component;
}
