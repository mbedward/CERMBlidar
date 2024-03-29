// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_overlap_midpoints
NumericMatrix get_overlap_midpoints(NumericMatrix xyf, double res, bool NSorientation);
RcppExport SEXP _CERMBlidar_get_overlap_midpoints(SEXP xyfSEXP, SEXP resSEXP, SEXP NSorientationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xyf(xyfSEXP);
    Rcpp::traits::input_parameter< double >::type res(resSEXP);
    Rcpp::traits::input_parameter< bool >::type NSorientation(NSorientationSEXP);
    rcpp_result_gen = Rcpp::wrap(get_overlap_midpoints(xyf, res, NSorientation));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CERMBlidar_get_overlap_midpoints", (DL_FUNC) &_CERMBlidar_get_overlap_midpoints, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CERMBlidar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
