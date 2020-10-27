// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// KalmanC
Rcpp::List KalmanC(const arma::mat& y, const Rcpp::LogicalMatrix& y_isna, const arma::colvec& a, const arma::mat& P_inf, const arma::mat& P_star, const arma::cube& Z, const arma::cube& T, const arma::cube& R, const arma::cube& Q, const bool& diagnostics);
RcppExport SEXP _statespacer_KalmanC(SEXP ySEXP, SEXP y_isnaSEXP, SEXP aSEXP, SEXP P_infSEXP, SEXP P_starSEXP, SEXP ZSEXP, SEXP TSEXP, SEXP RSEXP, SEXP QSEXP, SEXP diagnosticsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalMatrix& >::type y_isna(y_isnaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P_inf(P_infSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P_star(P_starSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const bool& >::type diagnostics(diagnosticsSEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanC(y, y_isna, a, P_inf, P_star, Z, T, R, Q, diagnostics));
    return rcpp_result_gen;
END_RCPP
}
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
// SimulateC
Rcpp::List SimulateC(const int& nsim, const int& repeat_Q, const int& N, const arma::colvec& a, const arma::cube& Z, const arma::cube& T, const arma::cube& R, const arma::cube& Q, const arma::mat& P_star, const bool& draw_initial, const bool& eta_only);
RcppExport SEXP _statespacer_SimulateC(SEXP nsimSEXP, SEXP repeat_QSEXP, SEXP NSEXP, SEXP aSEXP, SEXP ZSEXP, SEXP TSEXP, SEXP RSEXP, SEXP QSEXP, SEXP P_starSEXP, SEXP draw_initialSEXP, SEXP eta_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< const int& >::type repeat_Q(repeat_QSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P_star(P_starSEXP);
    Rcpp::traits::input_parameter< const bool& >::type draw_initial(draw_initialSEXP);
    Rcpp::traits::input_parameter< const bool& >::type eta_only(eta_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(SimulateC(nsim, repeat_Q, N, a, Z, T, R, Q, P_star, draw_initial, eta_only));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_statespacer_KalmanC", (DL_FUNC) &_statespacer_KalmanC, 10},
    {"_statespacer_LogLikC", (DL_FUNC) &_statespacer_LogLikC, 9},
    {"_statespacer_SimulateC", (DL_FUNC) &_statespacer_SimulateC, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_statespacer(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
