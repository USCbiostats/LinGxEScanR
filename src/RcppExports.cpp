// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// increment
int increment(arma::ivec& n);
RcppExport SEXP _LinGxEScanR_increment(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(increment(n));
    return rcpp_result_gen;
END_RCPP
}
// initlslinreg
int initlslinreg(const arma::vec& y, const arma::mat& xl, arma::mat& xtx, arma::mat& xtxinv0, arma::mat& ql, arma::mat& rtl, arma::vec& k, arma::vec& zt, arma::vec& resids, arma::vec& sigma2, arma::vec& s2, arma::vec& loglike);
RcppExport SEXP _LinGxEScanR_initlslinreg(SEXP ySEXP, SEXP xlSEXP, SEXP xtxSEXP, SEXP xtxinv0SEXP, SEXP qlSEXP, SEXP rtlSEXP, SEXP kSEXP, SEXP ztSEXP, SEXP residsSEXP, SEXP sigma2SEXP, SEXP s2SEXP, SEXP loglikeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xtx(xtxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xtxinv0(xtxinv0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type zt(ztSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type loglike(loglikeSEXP);
    rcpp_result_gen = Rcpp::wrap(initlslinreg(y, xl, xtx, xtxinv0, ql, rtl, k, zt, resids, sigma2, s2, loglike));
    return rcpp_result_gen;
END_RCPP
}
// lslinregfit
int lslinregfit(const arma::vec& dosage, const arma::vec& y, const arma::mat& xl, arma::mat& xr, arma::mat& bt, arma::mat& bb, const arma::mat& ql, arma::mat& qr, const arma::mat& rtl, arma::mat& rtr, arma::mat& rbr, arma::mat& h, const arma::mat& k, arma::mat& t, arma::mat& zb);
RcppExport SEXP _LinGxEScanR_lslinregfit(SEXP dosageSEXP, SEXP ySEXP, SEXP xlSEXP, SEXP xrSEXP, SEXP btSEXP, SEXP bbSEXP, SEXP qlSEXP, SEXP qrSEXP, SEXP rtlSEXP, SEXP rtrSEXP, SEXP rbrSEXP, SEXP hSEXP, SEXP kSEXP, SEXP tSEXP, SEXP zbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type bt(btSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtr(rtrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rbr(rbrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type zb(zbSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinregfit(dosage, y, xl, xr, bt, bb, ql, qr, rtl, rtr, rbr, h, k, t, zb));
    return rcpp_result_gen;
END_RCPP
}
// lslinregresiduals
int lslinregresiduals(const arma::vec& y, const arma::mat& xl, const arma::mat& xr, const arma::mat& bt, const arma::mat& bb, arma::vec& resids);
RcppExport SEXP _LinGxEScanR_lslinregresiduals(SEXP ySEXP, SEXP xlSEXP, SEXP xrSEXP, SEXP btSEXP, SEXP bbSEXP, SEXP residsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bt(btSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type resids(residsSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinregresiduals(y, xl, xr, bt, bb, resids));
    return rcpp_result_gen;
END_RCPP
}
// lslinregsigma2
int lslinregsigma2(const arma::mat& xl, const arma::mat& xr, const arma::vec& resids, arma::vec& sigma2, arma::vec& s2, arma::vec& loglike);
RcppExport SEXP _LinGxEScanR_lslinregsigma2(SEXP xlSEXP, SEXP xrSEXP, SEXP residsSEXP, SEXP sigma2SEXP, SEXP s2SEXP, SEXP loglikeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type loglike(loglikeSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinregsigma2(xl, xr, resids, sigma2, s2, loglike));
    return rcpp_result_gen;
END_RCPP
}
// lslinregxtx
int lslinregxtx(const arma::mat& xl, const arma::mat& xr, arma::mat& xtx);
RcppExport SEXP _LinGxEScanR_lslinregxtx(SEXP xlSEXP, SEXP xrSEXP, SEXP xtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xtx(xtxSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinregxtx(xl, xr, xtx));
    return rcpp_result_gen;
END_RCPP
}
// lslinregxtxinv
int lslinregxtxinv(const arma::mat& xtx, arma::mat& xtxinv);
RcppExport SEXP _LinGxEScanR_lslinregxtxinv(SEXP xtxSEXP, SEXP xtxinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xtx(xtxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xtxinv(xtxinvSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinregxtxinv(xtx, xtxinv));
    return rcpp_result_gen;
END_RCPP
}
// lslinregwaldtest
int lslinregwaldtest(const arma::mat& xl, const arma::mat& xr, const arma::mat& bb, const arma::vec& s2, const arma::mat& xtxinv, arma::vec& std_err, arma::mat& xrs2, arma::vec& chi2);
RcppExport SEXP _LinGxEScanR_lslinregwaldtest(SEXP xlSEXP, SEXP xrSEXP, SEXP bbSEXP, SEXP s2SEXP, SEXP xtxinvSEXP, SEXP std_errSEXP, SEXP xrs2SEXP, SEXP chi2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xtxinv(xtxinvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type std_err(std_errSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xrs2(xrs2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type chi2(chi2SEXP);
    rcpp_result_gen = Rcpp::wrap(lslinregwaldtest(xl, xr, bb, s2, xtxinv, std_err, xrs2, chi2));
    return rcpp_result_gen;
END_RCPP
}
// lslinreghwtest
int lslinreghwtest(const arma::mat& xl, const arma::mat& xr, const arma::vec& resids, const arma::mat& xtxinv, arma::vec& s2a, arma::mat& hws2);
RcppExport SEXP _LinGxEScanR_lslinreghwtest(SEXP xlSEXP, SEXP xrSEXP, SEXP residsSEXP, SEXP xtxinvSEXP, SEXP s2aSEXP, SEXP hws2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xtxinv(xtxinvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s2a(s2aSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hws2(hws2SEXP);
    rcpp_result_gen = Rcpp::wrap(lslinreghwtest(xl, xr, resids, xtxinv, s2a, hws2));
    return rcpp_result_gen;
END_RCPP
}
// lslinreg
int lslinreg(const arma::vec& dosage, const arma::vec& y, const arma::mat& xl, arma::mat& xr, arma::mat& xtx, arma::mat& bt, arma::mat& bb, const arma::mat& ql, arma::mat& qr, const arma::mat& rtl, arma::mat& rtr, arma::mat& rbr, arma::mat& h, const arma::mat& k, arma::mat& t, arma::mat& zb, arma::vec& resids, arma::vec& sigma2, arma::vec& s2, arma::mat& hws2, arma::mat& xtxinv, arma::vec& std_err, arma::mat& xrs2, arma::vec& chi2, arma::vec& loglike);
RcppExport SEXP _LinGxEScanR_lslinreg(SEXP dosageSEXP, SEXP ySEXP, SEXP xlSEXP, SEXP xrSEXP, SEXP xtxSEXP, SEXP btSEXP, SEXP bbSEXP, SEXP qlSEXP, SEXP qrSEXP, SEXP rtlSEXP, SEXP rtrSEXP, SEXP rbrSEXP, SEXP hSEXP, SEXP kSEXP, SEXP tSEXP, SEXP zbSEXP, SEXP residsSEXP, SEXP sigma2SEXP, SEXP s2SEXP, SEXP hws2SEXP, SEXP xtxinvSEXP, SEXP std_errSEXP, SEXP xrs2SEXP, SEXP chi2SEXP, SEXP loglikeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xl(xlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xtx(xtxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type bt(btSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type bb(bbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ql(qlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type qr(qrSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rtl(rtlSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rtr(rtrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type rbr(rbrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type zb(zbSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hws2(hws2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xtxinv(xtxinvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type std_err(std_errSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type xrs2(xrs2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type chi2(chi2SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type loglike(loglikeSEXP);
    rcpp_result_gen = Rcpp::wrap(lslinreg(dosage, y, xl, xr, xtx, bt, bb, ql, qr, rtl, rtr, rbr, h, k, t, zb, resids, sigma2, s2, hws2, xtxinv, std_err, xrs2, chi2, loglike));
    return rcpp_result_gen;
END_RCPP
}
// initreg
Rcpp::List initreg(const arma::mat& X, const arma::colvec& y);
RcppExport SEXP _LinGxEScanR_initreg(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(initreg(X, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LinGxEScanR_increment", (DL_FUNC) &_LinGxEScanR_increment, 1},
    {"_LinGxEScanR_initlslinreg", (DL_FUNC) &_LinGxEScanR_initlslinreg, 12},
    {"_LinGxEScanR_lslinregfit", (DL_FUNC) &_LinGxEScanR_lslinregfit, 15},
    {"_LinGxEScanR_lslinregresiduals", (DL_FUNC) &_LinGxEScanR_lslinregresiduals, 6},
    {"_LinGxEScanR_lslinregsigma2", (DL_FUNC) &_LinGxEScanR_lslinregsigma2, 6},
    {"_LinGxEScanR_lslinregxtx", (DL_FUNC) &_LinGxEScanR_lslinregxtx, 3},
    {"_LinGxEScanR_lslinregxtxinv", (DL_FUNC) &_LinGxEScanR_lslinregxtxinv, 2},
    {"_LinGxEScanR_lslinregwaldtest", (DL_FUNC) &_LinGxEScanR_lslinregwaldtest, 8},
    {"_LinGxEScanR_lslinreghwtest", (DL_FUNC) &_LinGxEScanR_lslinreghwtest, 6},
    {"_LinGxEScanR_lslinreg", (DL_FUNC) &_LinGxEScanR_lslinreg, 25},
    {"_LinGxEScanR_initreg", (DL_FUNC) &_LinGxEScanR_initreg, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_LinGxEScanR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
