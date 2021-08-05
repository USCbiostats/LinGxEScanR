// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

const double LOG2PIP1 = log(M_PI + M_PI) + 1.;

// [[Rcpp::export]]
int increment(arma::ivec &n) {
  n[0]++;
  return 0;
}

// [[Rcpp::export]]
int initlslinreg(const arma::vec &y,
                 const arma::mat &xl,
                 arma::mat &xtx,
                 arma::mat &xtxinv0,
                 arma::mat &ql,
                 arma::mat &rtl,
                 arma::vec &k,
                 arma::vec &zt,
                 arma::vec &resids,
                 arma::vec &sigma2,
                 arma::vec &s2,
                 arma::vec &loglike) {
  qr_econ(ql, rtl, xl);
  zt = ql.t() * y;
  solve(k, rtl, zt);
  resids = y - xl * k;
  sigma2[0] = arma::dot(resids, resids) / xl.n_rows;
  s2[0] = sigma2[0] * xl.n_rows / (xl.n_rows - xl.n_cols);
  xtx.submat(0, 0, xl.n_cols - 1, xl.n_cols - 1) = trans(xl) * xl;
  xtxinv0 = pinv(xtx.submat(0, 0, xl.n_cols - 1, xl.n_cols - 1));
  loglike[0] = -0.5 * xl.n_rows * (log(sigma2[0]) + LOG2PIP1);

  return 0;
}

// [[Rcpp::export]]
int lslinregfit(const arma::vec &dosage,
                const arma::vec &y,
                const arma::mat &xl,
                arma::mat &xr,
                arma::mat &bt,
                arma::mat &bb,
                const arma::mat &ql,
                arma::mat &qr,
                const arma::mat &rtl,
                arma::mat &rtr,
                arma::mat &rbr,
                arma::mat &h,
                const arma::mat &k,
                arma::mat &t,
                arma::mat &zb) {
  xr.col(0) = dosage;
  if (xr.n_cols == 2)
    xr.col(1) = xr.col(0) % xl.col(xl.n_cols - 1);
  rtr = trans(ql) * xr;
  t = xr - ql * rtr;
  qr_econ(qr, rbr, t);
  solve(h, rtl, rtr);
  zb = trans(qr) * y;
  arma::mat bbt = bb.t();
  solve(bbt, rbr, zb);
  bt = k - h*bbt;
  bb = bbt.t();

  return 0;
}

// [[Rcpp::export]]
int lslinregresiduals(const arma::vec &y,
                      const arma::mat &xl,
                      const arma::mat &xr,
                      const arma::mat &bt,
                      const arma::mat &bb,
                      arma::vec &resids) {
  resids = y - xl*bt - xr*bb;
  
  return 0;
}

// [[Rcpp::export]]
int lslinregsigma2(const arma::mat &xl,
                   const arma::mat &xr,
                   const arma::vec &resids,
                   arma::vec &sigma2,
                   arma::vec &s2,
                   arma::vec &loglike) {
  sigma2[1] = std::inner_product(resids.begin(), resids.end(), resids.begin(), 0.0);
  s2[1] = sigma2[1] / (xl.n_rows - xl.n_cols - xr.n_cols);
  sigma2[1] /= xl.n_rows;
  loglike[1] = -0.5 * xl.n_rows * (log(sigma2[1]) + LOG2PIP1);
  
  return 0;
}

// [[Rcpp::export]]
int lslinregxtx(const arma::mat &xl,
                const arma::mat &xr,
                arma::mat &xtx) {
  int c1, c2, r1, r2;

  r1 = xl.n_cols;
  r2 = xtx.n_cols - 1;
  c1 = 0;
  c2 = xl.n_cols - 1;
  
  xtx.submat(r1, c1, r2, c2) = xr.t() * xl;
  xtx.submat(c1, r1, c2, r2) = xtx.submat(r1, c1, r2, c2).t();
  xtx.submat(r1, r1, r2, r2) = xr.t() * xr;
  
  return 0;
}

// [[Rcpp::export]]
int lslinregxtxinv(const arma::mat &xtx,
                   arma::mat &xtxinv) {
  xtxinv = pinv(xtx);
  
  return 0;
}

// [[Rcpp::export]]
int lslinregwaldtest(const arma::mat &xl,
                     const arma::mat &xr,
                     const arma::mat &bb,
                     const arma::vec &s2,
                     const arma::mat &xtxinv,
                     arma::vec &std_err,
                     arma::mat &xrs2,
                     arma::vec &chi2) {
  int r1, r2;
  
  r1 = xl.n_cols;
  r2 = xtxinv.n_cols - 1;

  std_err = arma::sqrt(s2[1] * arma::diagvec(xtxinv));
  xrs2 = s2[1] * xtxinv.submat(r1,r1,r2,r2);
  chi2 = bb.t() * pinv(xrs2) * bb;

  return 0;
}

// [[Rcpp::export]]
int lslinreghwtest(const arma::mat &xl,
                   const arma::mat &xr,
                   const arma::vec &resids,
                   const arma::mat &xtxinv,
                   arma::mat &s2a,
                   arma::mat &hws2) {
  int r1, r2, r3;
  
  r1 = xl.n_cols - 1;
  r2 = xl.n_cols;
  r3 = xtxinv.n_cols - 1;
  
  s2a.diag() = resids % resids;
  
  hws2.submat(0, 0, r1, r1) = xl.t() * s2a * xl;
  hws2.submat(0, r2, r1, r3) = xl.t() * s2a * xr;
  hws2.submat(r2, 0, r3, r1) = hws2.submat(0, r2, r1, r3).t();
  hws2.submat(r2, r2, r3, r3) = xr.t() * s2a * xr;
  hws2.submat(0, 0, r3, r3) = xtxinv * hws2 * xtxinv;
  
  return 0;
}

// [[Rcpp::export]]
int lslinreg(const arma::vec &dosage,
             const arma::vec &y,
             const arma::mat &xl,
             arma::mat &xr,
             arma::mat &xtx,
             arma::mat &bt,
             arma::mat &bb,
             const arma::mat &ql,
             arma::mat &qr,
             const arma::mat &rtl,
             arma::mat &rtr,
             arma::mat &rbr,
             arma::mat &h,
             const arma::mat &k,
             arma::mat &t,
             arma::mat &zb,
             arma::vec &resids,
             arma::vec &sigma2,
             arma::vec &s2,
             arma::mat &hws2,
             arma::mat &xtxinv,
             arma::vec &std_err,
             arma::mat &xrs2,
             arma::vec &chi2,
             arma::vec &loglike) {
  int c1, c2, r1, r2;

  xr.col(0) = dosage;
  if (xr.n_cols == 2)
    xr.col(1) = xr.col(0) % xl.col(xl.n_cols - 1);
  rtr = trans(ql) * xr;
  t = xr - ql * rtr;
  qr_econ(qr, rbr, t);
  solve(h, rtl, rtr);
  zb = trans(qr) * y;
  arma::mat bbt = bb;
  solve(bbt, rbr, zb);
  bt = k - h*bbt;
  bb = bbt;

  resids = y - xl*bt - xr*bb;
  sigma2[1] = std::inner_product(resids.begin(), resids.end(), resids.begin(), 0.0);
  s2[1] = sigma2[1] / (xl.n_rows - xl.n_cols - xr.n_cols);
  sigma2[1] /= xl.n_rows;

  loglike[1] = -0.5 * xl.n_rows * (log(sigma2[1]) + LOG2PIP1);
  
  // Building the xtx matrix  
  r1 = xl.n_cols;
  r2 = xtx.n_cols - 1;
  c1 = 0;
  c2 = xl.n_cols - 1;
  
  xtx.submat(r1, c1, r2, c2) = trans(xr) * xl;
  xtx.submat(c1, r1, c2, r2) = trans(xtx.submat(r1, c1, r2, c2));
  xtx.submat(r1, r1, r2, r2) = trans(xr) * xr;
  
  xtxinv = pinv(xtx);
  std_err = arma::sqrt(s2[1] * arma::diagvec(xtxinv));
  xrs2 = s2[1] * xtxinv.submat(r1,r1,r2,r2);
  chi2 = bb.t() * pinv(xrs2) * bb;
  
  hws2.submat(0, 0, xl.n_cols - 1, xl.n_cols - 1)
    = xl.t() * arma::diagmat(resids % resids) * xl;
  hws2.submat(0, xl.n_cols, xl.n_cols - 1, xtx.n_cols - 1)
    = xl.t() * arma::diagmat(resids % resids) * xr;
  hws2.submat(xl.n_cols, 0, xtx.n_cols - 1, xl.n_cols - 1)
    = hws2.submat(0, xl.n_cols, xl.n_cols - 1, xtx.n_cols - 1).t();
  hws2.submat(xl.n_cols, xl.n_cols, xtx.n_cols - 1, xtx.n_cols - 1)
    = xr.t() * arma::diagmat(resids % resids) * xr;
  hws2.submat(0, 0, xtx.n_cols - 1, xtx.n_cols - 1) = xtxinv * hws2 * xtxinv;

  return 0;
}

// [[Rcpp::export]]
Rcpp::List initreg(const arma::mat &X, const arma::colvec &y) {
  int n = X.n_rows, k = X.n_cols;
  
  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals
  
  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
  
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));  
  
  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = std_err,
                            Rcpp::Named("df.residual")  = n - k,
                            Rcpp::Named("s2") = s2);
}
