#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double pBH_cutoff(vec &p, double q);
double lfdr_cutoff(vec &lfdr, double q);
  
// [[Rcpp::export]]
SEXP cutoff(SEXP P, SEXP rLIS, SEXP Lfdr, SEXP covLfdr, SEXP carLIS, SEXP q_in) {
  vec p = as<arma::vec>(P);
  vec lis = as<arma::vec>(rLIS);
  vec lfdr = as<arma::vec>(Lfdr);
  vec covlfdr = as<arma::vec>(covLfdr);
  vec carlis = as<arma::vec>(carLIS);
  double q = Rcpp::as<double>(q_in);
  
  double pBH_thr, rLIS_thr, Lfdr_thr, covLfdr_thr, carLIS_thr;
  pBH_thr = pBH_cutoff(p, q);
  rLIS_thr = lfdr_cutoff(lis, q);
  Lfdr_thr = lfdr_cutoff(lfdr, q);
  covLfdr_thr = lfdr_cutoff(covlfdr, q);
  carLIS_thr = lfdr_cutoff(carlis, q);
  
  return Rcpp::List::create(Rcpp::Named("pBH_thr") = pBH_thr,
                            Rcpp::Named("rLIS_thr") = rLIS_thr,
                            Rcpp::Named("Lfdr_thr") = Lfdr_thr,
                            Rcpp::Named("covLfdr_thr") = covLfdr_thr,
                            Rcpp::Named("carLIS_thr") = carLIS_thr);
}

double pBH_cutoff(vec &p, double q){
  int J = p.size();
  vec sort_p = sort(p);
  
  double thr = 0, fdp = 1;
  for(int i = J; i > 0; i--){
    thr = sort_p(i-1);
    fdp = thr * J / i;
    if(fdp <= q)
      return thr;
  }
  return 0;
}

double lfdr_cutoff(vec &lfdr, double q){
  vec sort_lfdr = sort(lfdr);
  vec cum_lfdr = cumsum(lfdr);

  double thr = 0, fdp = 1;
  for(int i = lfdr.size(); i > 0; i--){
    // std::cout << i << std::endl;
    thr = sort_lfdr(i-1);
    // fdp = sum(sort_rLIS(span(0, i-1))) / i;
    fdp = cum_lfdr(i-1) / i;
    if(fdp <= q)
      return thr;
  }
  return 0;
}
