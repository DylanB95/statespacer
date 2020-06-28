// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LogLikC
double LogLikC(const Rcpp::NumericMatrix& y, const Rcpp::LogicalMatrix& y_isna, arma::colvec a, arma::mat P_inf, arma::mat P_star, const arma::cube& Z, const arma::cube& T, const arma::cube& R, const arma::cube& Q);
RcppExport SEXP _statespacer_LogLikC(SEXP ySEXP, SEXP y_isnaSEXP, SEXP aSEXP, SEXP P_infSEXP, SEXP P_starSEXP, SEXP ZSEXP, SEXP TSEXP, SEXP RSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalMatrix& >::type y_isna(y_isnaSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_inf(P_infSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P_star(P_starSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(LogLikC(y, y_isna, a, P_inf, P_star, Z, T, R, Q));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_statespacer_LogLikC", (DL_FUNC) &_statespacer_LogLikC, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_statespacer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}