// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matCov
SEXP matCov(const arma::mat A, const arma::mat B, const bool corMat, const double eps);
RcppExport SEXP _MGMM_matCov(SEXP ASEXP, SEXP BSEXP, SEXP corMatSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< const bool >::type corMat(corMatSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(matCov(A, B, corMat, eps));
    return rcpp_result_gen;
END_RCPP
}
// eigSym
SEXP eigSym(const arma::mat A);
RcppExport SEXP _MGMM_eigSym(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(eigSym(A));
    return rcpp_result_gen;
END_RCPP
}
// matDet
SEXP matDet(const arma::mat A, const bool logDet);
RcppExport SEXP _MGMM_matDet(SEXP ASEXP, SEXP logDetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const bool >::type logDet(logDetSEXP);
    rcpp_result_gen = Rcpp::wrap(matDet(A, logDet));
    return rcpp_result_gen;
END_RCPP
}
// matInv
SEXP matInv(const arma::mat A);
RcppExport SEXP _MGMM_matInv(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matInv(A));
    return rcpp_result_gen;
END_RCPP
}
// matIP
SEXP matIP(const arma::mat A, const arma::mat B);
RcppExport SEXP _MGMM_matIP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matIP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// MMP
SEXP MMP(const arma::mat A, const arma::mat B);
RcppExport SEXP _MGMM_MMP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MMP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matOP
SEXP matOP(const arma::mat A, const arma::mat B);
RcppExport SEXP _MGMM_matOP(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matOP(A, B));
    return rcpp_result_gen;
END_RCPP
}
// matQF
SEXP matQF(const arma::mat X, const arma::mat A);
RcppExport SEXP _MGMM_matQF(SEXP XSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(matQF(X, A));
    return rcpp_result_gen;
END_RCPP
}
// SchurC
SEXP SchurC(const arma::mat Ibb, const arma::mat Iaa, const arma::mat Iba);
RcppExport SEXP _MGMM_SchurC(SEXP IbbSEXP, SEXP IaaSEXP, SEXP IbaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Ibb(IbbSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Iaa(IaaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Iba(IbaSEXP);
    rcpp_result_gen = Rcpp::wrap(SchurC(Ibb, Iaa, Iba));
    return rcpp_result_gen;
END_RCPP
}
// tr
SEXP tr(const arma::mat A);
RcppExport SEXP _MGMM_tr(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(tr(A));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MGMM_matCov", (DL_FUNC) &_MGMM_matCov, 4},
    {"_MGMM_eigSym", (DL_FUNC) &_MGMM_eigSym, 1},
    {"_MGMM_matDet", (DL_FUNC) &_MGMM_matDet, 2},
    {"_MGMM_matInv", (DL_FUNC) &_MGMM_matInv, 1},
    {"_MGMM_matIP", (DL_FUNC) &_MGMM_matIP, 2},
    {"_MGMM_MMP", (DL_FUNC) &_MGMM_MMP, 2},
    {"_MGMM_matOP", (DL_FUNC) &_MGMM_matOP, 2},
    {"_MGMM_matQF", (DL_FUNC) &_MGMM_matQF, 2},
    {"_MGMM_SchurC", (DL_FUNC) &_MGMM_SchurC, 3},
    {"_MGMM_tr", (DL_FUNC) &_MGMM_tr, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MGMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
