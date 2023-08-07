// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// im_equ
Rcpp::IntegerVector im_equ(Eigen::MatrixXd A, Eigen::MatrixXd b);
RcppExport SEXP _polysf_im_equ(SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(im_equ(A, b));
    return rcpp_result_gen;
END_RCPP
}
// lasso_fr
Rcpp::IntegerVector lasso_fr(Eigen::MatrixXd Omega, Eigen::MatrixXd pi);
RcppExport SEXP _polysf_lasso_fr(SEXP OmegaSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_fr(Omega, pi));
    return rcpp_result_gen;
END_RCPP
}
// ind_rows
Rcpp::List ind_rows(Eigen::MatrixXd M);
RcppExport SEXP _polysf_ind_rows(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(ind_rows(M));
    return rcpp_result_gen;
END_RCPP
}
// add_to_full
Eigen::MatrixXd add_to_full(Eigen::MatrixXd N);
RcppExport SEXP _polysf_add_to_full(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(add_to_full(N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polysf_im_equ", (DL_FUNC) &_polysf_im_equ, 2},
    {"_polysf_lasso_fr", (DL_FUNC) &_polysf_lasso_fr, 2},
    {"_polysf_ind_rows", (DL_FUNC) &_polysf_ind_rows, 1},
    {"_polysf_add_to_full", (DL_FUNC) &_polysf_add_to_full, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_polysf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}