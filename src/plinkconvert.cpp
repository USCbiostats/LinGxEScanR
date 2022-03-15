#include <Rcpp.h>

// [[Rcpp::export]]
int snprawtoint(Rcpp::RawVector &r, Rcpp::IntegerVector &d) {
  Rcpp::RawVector::iterator rit;
  Rcpp::IntegerVector::iterator dit;
  int i, j;
  unsigned char mask[4] = {0x03, 0x0c, 0x30, 0xc0};
  int rawtointmap[4] = {0, NA_INTEGER, 1, 2};
  
  rit = r.begin();
  for (dit = d.begin(); dit != d.end(); ++rit) {
    for (i = 0, j = 0; i < 4 && dit != d.end(); ++i, j += 2, ++dit) {
      *dit = rawtointmap[(*rit & mask[i]) >> j];
    }
  }
  return 0;
}

// [[Rcpp::export]]
int snpinttoraw(Rcpp::IntegerVector &d, Rcpp::RawVector &r) {
  Rcpp::RawVector::iterator rit;
  Rcpp::IntegerVector::iterator dit;
  int i, j;
  int inttorawmap[3] = {0x00, 0x01, 0x03};
  
  rit = r.begin();
  for (dit = d.begin(); dit != d.end(); ++rit) {
    *rit = 0x00;
    for (i = 0, j = 0; i < 4 && dit != d.end(); ++i, j += 2, ++dit) {
      if (*dit == NA_INTEGER)
        *rit |= (0x10 << j);
      else
        *rit |= (inttorawmap[*dit] << j);
    }
  }
  return 0;
}
