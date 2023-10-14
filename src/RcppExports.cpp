// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// spmaC
Rcpp::List spmaC(double alpha, double epi_abs, double epi_rel, double rho, double lambda, int iter, NumericMatrix covF1, NumericMatrix covF2, NumericMatrix covF3, arma::mat augIMat, int p, int M, IntegerVector setN, arma::mat thetaF);
RcppExport SEXP _SPMA_spmaC(SEXP alphaSEXP, SEXP epi_absSEXP, SEXP epi_relSEXP, SEXP rhoSEXP, SEXP lambdaSEXP, SEXP iterSEXP, SEXP covF1SEXP, SEXP covF2SEXP, SEXP covF3SEXP, SEXP augIMatSEXP, SEXP pSEXP, SEXP MSEXP, SEXP setNSEXP, SEXP thetaFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type epi_abs(epi_absSEXP);
    Rcpp::traits::input_parameter< double >::type epi_rel(epi_relSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covF1(covF1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covF2(covF2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covF3(covF3SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type augIMat(augIMatSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type setN(setNSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type thetaF(thetaFSEXP);
    rcpp_result_gen = Rcpp::wrap(spmaC(alpha, epi_abs, epi_rel, rho, lambda, iter, covF1, covF2, covF3, augIMat, p, M, setN, thetaF));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SPMA_spmaC", (DL_FUNC) &_SPMA_spmaC, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_SPMA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
