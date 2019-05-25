// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// timesThree
int timesThree(int x);
RcppExport SEXP _q2e_timesThree(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesThree(x));
    return rcpp_result_gen;
END_RCPP
}
// rq2e
int rq2e(Rcpp::StringVector argVector, Rcpp::StringVector fnVector, int Rnreplicates);
RcppExport SEXP _q2e_rq2e(SEXP argVectorSEXP, SEXP fnVectorSEXP, SEXP RnreplicatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type argVector(argVectorSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type fnVector(fnVectorSEXP);
    Rcpp::traits::input_parameter< int >::type Rnreplicates(RnreplicatesSEXP);
    rcpp_result_gen = Rcpp::wrap(rq2e(argVector, fnVector, Rnreplicates));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_q2e_timesThree", (DL_FUNC) &_q2e_timesThree, 1},
    {"_q2e_rq2e", (DL_FUNC) &_q2e_rq2e, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_q2e(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
