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
                 arma::mat &ql,
                 arma::mat &rtl,
                 arma::vec &k,
                 arma::vec &zt,
                 arma::vec &resids,
                 arma::vec &s2) {
  double sigma2;
  
  qr_econ(ql, rtl, xl);
  zt = ql.t() * y;
  solve(k, rtl, zt);
  resids = y - xl * k;
  sigma2 = arma::dot(resids, resids) / xl.n_rows;
  s2[0] = sigma2 * xl.n_rows / (xl.n_rows - xl.n_cols);
  xtx.submat(0, 0, xl.n_cols - 1, xl.n_cols - 1) = trans(xl) * xl;

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
             arma::vec &s2,
             arma::mat &xtxinv,
             arma::vec &std_err,
             arma::mat &xrs2,
             arma::vec &chi2) {
  int c1, c2, r1, r2;

  xr.col(0) = dosage;
  if (xr.n_cols == 2)
    xr.col(1) = xr.col(0) % xl.col(xl.n_cols - 1);
  rtr = trans(ql) * xr;
  t = xr - ql * rtr;
  qr_econ(qr, rbr, t);
  solve(h, rtl, rtr);
  zb = trans(qr) * y;
  solve(bb, rbr, zb);
  bt = k - h*bb;
  
  resids = y - xl*bt - xr*bb;
  s2[0] = std::inner_product(resids.begin(), resids.end(), resids.begin(), 0.0);
  s2[0] /= (xl.n_rows - xl.n_cols - xr.n_cols);
  
  // Building the xtx matrix  
  r1 = xl.n_cols;
  r2 = xtx.n_cols - 1;
  c1 = 0;
  c2 = xl.n_cols - 1;
  
  xtx.submat(r1, c1, r2, c2) = trans(xr) * xl;
  xtx.submat(c1, r1, c2, r2) = trans(xtx.submat(r1, c1, r2, c2));
  xtx.submat(r1, r1, r2, r2) = trans(xr) * xr;
  
  xtxinv = s2[0] * pinv(xtx);
  std_err = arma::sqrt(arma::diagvec(xtxinv));  
  xrs2 = xtxinv.submat(r1,r1,r2,r2);
  chi2 = trans(bb) * pinv(xrs2) * bb;
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
