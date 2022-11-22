// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cxxMixSquarem
List cxxMixSquarem(NumericMatrix matrix_lik, NumericVector prior, NumericVector pi_init, List control);
RcppExport SEXP _ashr_cxxMixSquarem(SEXP matrix_likSEXP, SEXP priorSEXP, SEXP pi_initSEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type matrix_lik(matrix_likSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pi_init(pi_initSEXP);
    Rcpp::traits::input_parameter< List >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(cxxMixSquarem(matrix_lik, prior, pi_init, control));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ashr_cxxMixSquarem", (DL_FUNC) &_ashr_cxxMixSquarem, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ashr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
