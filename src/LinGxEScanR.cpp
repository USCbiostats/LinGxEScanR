// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

const double LOG2PIP1 = log(M_PI + M_PI) + 1.;

// [[Rcpp::export]]
Rcpp::List initlslinreg(const arma::vec &y,
                        const arma::mat &x) {
  arma::mat ql, rtl, xtx;
  arma::vec residuals, zt, k;
  double sigma2, s2, loglike;
  
  qr_econ(ql, rtl, x);
  zt = ql.t() * y;
  solve(k, rtl, zt);
  residuals = y - x * k;
  sigma2 = arma::dot(residuals, residuals) / x.n_rows;
  s2 = sigma2 * x.n_rows / (x.n_rows - x.n_cols);
  xtx = trans(x) * x;
  loglike = -(y.n_elem * (LOG2PIP1 + log(sigma2)) / 2.);
  
  return Rcpp::List::create(Rcpp::Named("residuals") = residuals,
                            Rcpp::Named("ql") = ql,
                            Rcpp::Named("rtl") = rtl,
                            Rcpp::Named("zt") = zt,
                            Rcpp::Named("k") = k,
                            Rcpp::Named("s2") = s2,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("xtx") = xtx,
                            Rcpp::Named("loglike") = loglike);
}

// [[Rcpp::export]]
int lslinreg(const arma::vec &y,
             const arma::mat &xl,
             const arma::mat &xr,
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
             arma::vec &res,
             arma::vec &s2,
             arma::mat &xtxinv,
             arma::vec &std_err,
             arma::mat &xrS2,
             arma::vec &xrChi2) {
  int c1, c2, r1, r2;

  rtr = trans(ql) * xr;
  t = xr - ql * rtr;
  qr_econ(qr, rbr, t);
  solve(h, rtl, rtr);
  zb = trans(qr) * y;
  solve(bb, rbr, zb);
  bt = k - h*bb;
  
  res = y - xl*bt - xr*bb;
  s2[0] = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
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
  xrS2 = xtxinv.submat(r1,r1,r2,r2);
  xrChi2 = trans(bb) * pinv(xrS2) * bb;
  return 0;
}

// and we can use Rcpp::List to return both at the same time
//
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
