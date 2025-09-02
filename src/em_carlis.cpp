#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

vec my_pava(vec &values, vec &weight, bool decreasing);
vec carlis(vec &pi, mat &A, vec &f1, vec &f2, vec &f3);
vec fdr(vec &carLIS);

// [[Rcpp::export]]
SEXP em_carlis(SEXP pa_in, SEXP pb_in, SEXP x_in, SEXP pi0a_in, SEXP pi0b_in, SEXP pi0x_in) {
  try{
    const int maxIter = 200;
    const double tol = 1e-3;
    const double pvalue_cutoff = 1e-15;
    const double f_cutoff = 1e-15;
    
    vec pa = as<arma::vec>(pa_in), pb = as<arma::vec>(pb_in), x = as<arma::vec>(x_in);
    const double pi0_pa = Rcpp::as<double>(pi0a_in), pi0_pb = Rcpp::as<double>(pi0b_in), pi0_x = Rcpp::as<double>(pi0x_in);
    
    int m = pa.size();
    double min_a = pa.elem(find(pa>0)).min(), min_b = pb.elem(find(pb>0)).min(), min_x = x.elem(find(x>0)).min();
    pa.elem(find(pa<=0)).fill(pvalue_cutoff < min_a ? pvalue_cutoff : min_a);
    pb.elem(find(pb<=0)).fill(pvalue_cutoff < min_b ? pvalue_cutoff : min_b);
    x.elem(find(x<=0)).fill(pvalue_cutoff < min_x ? pvalue_cutoff : min_x);
    
    vec p1 = sort(pa), p2 = sort(pb), p3 = sort(x);
    uvec ix1 = sort_index(pa), ix2 = sort_index(pb), ix3 = sort_index(x);
    
    vec p1_diff(m), p2_diff(m), p3_diff(m);
    p1_diff(0) = p1(0);
    p2_diff(0) = p2(0);
    p3_diff(0) = p3(0);
    for(int i = 1; i<m; ++i){
      p1_diff(i) = p1(i) - p1(i-1);
      p2_diff(i) = p2(i) - p2(i-1);
      p3_diff(i) = p3(i) - p3(i-1);
    }
    
    // Initialization
    vec pi(8);
    pi(0) = pi0_pa * pi0_pb * pi0_x; 
    pi(1) = pi0_pa * pi0_pb * (1 - pi0_x);
    pi(2) = pi0_pa * (1 - pi0_pb) * pi0_x; 
    pi(3) = pi0_pa * (1 - pi0_pb) * (1 - pi0_x);
    pi(4) = (1 - pi0_pa) * pi0_pb * pi0_x;
    pi(5) = (1 - pi0_pa) * pi0_pb * (1 - pi0_x);
    pi(6) = (1 - pi0_pa) * (1 - pi0_pb) * pi0_x;
    pi(7) = (1 - pi0_pa) * (1 - pi0_pb) * (1 - pi0_x);
    
    mat A = {{0.72, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04},
    {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04},
    {0.04, 0.04, 0.72, 0.04, 0.04, 0.04, 0.04, 0.04},
    {0.04, 0.04, 0.04, 0.72, 0.04, 0.04, 0.04, 0.04},
    {0.04, 0.04, 0.04, 0.04, 0.72, 0.04, 0.04, 0.04},
    {0.04, 0.04, 0.04, 0.04, 0.04, 0.72, 0.04, 0.04},
    {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.72, 0.04},
    {0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.72}};
    
    vec f0 = ones(m,1), f1(m), f2(m), f3(m);
    f1 = 1 - pa;
    f2 = 1 - pb;
    f3 = 1 - x;
    
    mat f = join_rows(join_rows(f0, f3, f2, f2%f3), join_rows(f1, f1%f3, f1%f2, f1%f2%f3));
    
    vec loglik(maxIter);
    loglik(0) = -datum::inf;
    
    mat alpha = zeros(m, 8), beta = zeros(m, 8), gamma = zeros(m, 8);
    cube xi = zeros(8, 8, m-1);
    
    vec xi_kl(m-1), Q1(m), Q2(m), Q3(m), y1(m), y2(m), y3(m), res1(m), res2(m), res3(m);
    mat xi_k(8,m-1);
    mat f_jp1(8,8);
    double sub_loglik, loglik_delta;
    
    std::cout << "EM begins:" << std::endl;
    
    for (int i = 1; i < maxIter; i++){
      // E-step
      // calculate the forward, backward and posterior probabilities based on current HMM
      alpha.row(0) = pi.t() % f.row(0);
      alpha.row(0) = alpha.row(0)/sum(alpha.row(0));
      
      for (int j = 1; j < m; j++){
        alpha(j,0) = sum(alpha.row(j-1).t() % A.col(0)) * f(j,0);
        alpha(j,1) = sum(alpha.row(j-1).t() % A.col(1)) * f(j,1);
        alpha(j,2) = sum(alpha.row(j-1).t() % A.col(2)) * f(j,2);
        alpha(j,3) = sum(alpha.row(j-1).t() % A.col(3)) * f(j,3);
        alpha(j,4) = sum(alpha.row(j-1).t() % A.col(4)) * f(j,4);
        alpha(j,5) = sum(alpha.row(j-1).t() % A.col(5)) * f(j,5);
        alpha(j,6) = sum(alpha.row(j-1).t() % A.col(6)) * f(j,6);
        alpha(j,7) = sum(alpha.row(j-1).t() % A.col(7)) * f(j,7);
        
        alpha.row(j) = alpha.row(j)/sum(alpha.row(j));
      }
      
      beta.row(m-1).fill(0.125);
      for(int j = m-2; j >=0; j--){
        beta(j,0) = sum(beta.row(j+1) % A.row(0) % f.row(j+1));
        beta(j,1) = sum(beta.row(j+1) % A.row(1) % f.row(j+1));
        beta(j,2) = sum(beta.row(j+1) % A.row(2) % f.row(j+1));
        beta(j,3) = sum(beta.row(j+1) % A.row(3) % f.row(j+1));
        beta(j,4) = sum(beta.row(j+1) % A.row(4) % f.row(j+1));
        beta(j,5) = sum(beta.row(j+1) % A.row(5) % f.row(j+1));
        beta(j,6) = sum(beta.row(j+1) % A.row(6) % f.row(j+1));
        beta(j,7) = sum(beta.row(j+1) % A.row(7) % f.row(j+1));
        
        beta.row(j) = beta.row(j)/sum(beta.row(j));
      }
      
      for(int j = 0; j < m; j++){
        gamma.row(j) = alpha.row(j) % beta.row(j) / sum(alpha.row(j) % beta.row(j));
      }
      
      for(int j = 0; j < m-1; j++){
        f_jp1 = repmat(f.row(j+1), 8, 1);
        xi.slice(j) = alpha.row(j).t() * beta.row(j+1) % A % f_jp1/
          accu(alpha.row(j).t() * beta.row(j+1) % A % f_jp1);
      }
      
      // M-step
      // update the parameters pi and A
      pi = gamma.row(0).t();
      for(int k = 0; k < 8; k++){
        for(int l = 0; l < 8; l++){
          xi_kl = xi.tube(k,l);
          xi_k = xi.row(k);
          A(k,l) = sum(xi_kl)/accu(xi_k);
        }
      }
      
      // update f1 and f2
      Q1 = gamma.col(4) + gamma.col(5) + gamma.col(6) + gamma.col(7);
      Q2 = gamma.col(2) + gamma.col(3) + gamma.col(6) + gamma.col(7);
      Q3 = gamma.col(1) + gamma.col(3) + gamma.col(5) + gamma.col(7);
      Q1 = Q1(ix1);
      Q2 = Q2(ix2);
      Q3 = Q3(ix3);
      
      y1 = - p1_diff * sum(Q1) / Q1;
      y2 = - p2_diff * sum(Q2) / Q2;
      y3 = - p3_diff * sum(Q3) / Q3;
      
      y1.elem(find_nonfinite(y1)).fill(y1.elem(find_finite(y1)).min());
      y2.elem(find_nonfinite(y2)).fill(y2.elem(find_finite(y2)).min());
      y3.elem(find_nonfinite(y3)).fill(y3.elem(find_finite(y3)).min());
      
      res1 = my_pava(y1, Q1, true);
      res2 = my_pava(y2, Q2, true);
      res3 = my_pava(y3, Q3, true);
      
      f1 = -1 / res1;
      f1 = f1 / sum(f1 % p1_diff);
      f1(ix1) = f1;
      f1.elem(find_nan(f1)).fill(f1.min());
      
      f2 = -1 / res2;
      f2 = f2 / sum(f2 % p2_diff);
      f2(ix2) = f2;
      f2.elem(find_nan(f2)).fill(f2.min());

      f3 = -1 / res3;
      f3 = f3 / sum(f3 % p3_diff);
      f3(ix3) = f3;
      f3.elem(find_nan(f3)).fill(f3.min());
      
      double min_f1 = f1.elem(find(f1>0)).min(), min_f2 = f2.elem(find(f2>0)).min(), min_f3 = f3.elem(find(f3>0)).min();
      f1.elem(find(f1<=0)).fill(f_cutoff < min_f1 ? f_cutoff : min_f1);
      f2.elem(find(f2<=0)).fill(f_cutoff < min_f2 ? f_cutoff : min_f2);
      f3.elem(find(f3<=0)).fill(f_cutoff < min_f3 ? f_cutoff : min_f3);
      
      f = join_rows(join_rows(f0, f3, f2, f2%f3), join_rows(f1, f1%f3, f1%f2, f1%f2%f3));
      
      // calculate the updated log-likelihood
      sub_loglik = 0;
      for(int j = 0; j < m - 1; j++){
        sub_loglik = sub_loglik + accu(log(A) % xi.slice(j));
      }
      
      loglik(i) = sum(log(pi).t() % gamma.row(1)) + sub_loglik + accu(gamma % log(f));
      loglik_delta = abs((loglik(i) - loglik(i-1))/loglik(i-1));
      
      std::cout<<i<<". "<< loglik(i) << ", delta = " << loglik_delta << std::endl;
      
      if(loglik_delta < tol){
        break;
      }
    }
    
    vec carLIS = carlis(pi, A, f1, f2, f3);
    vec caradj = fdr(carLIS);
    
    return Rcpp::List::create(Rcpp::Named("carLIS") = carLIS,
                              Rcpp::Named("fdr") = caradj,
                              Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("pi") = pi,
                              Rcpp::Named("A") = A,
                              Rcpp::Named("f1") = f1,
                              Rcpp::Named("f2") = f2,
                              Rcpp::Named("f3") = f3);
  } catch( std::exception &ex ) {
    forward_exception_to_r(ex);
    return Rcpp::List::create();
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
    return Rcpp::List::create();
  }
}

double na_rm(vec &p){
  double _min = std::numeric_limits<double>::max();
  
  p.for_each([&_min](double &val) { if(val > 0 && val < _min) _min = val; });
  
  return _min;
}

vec carlis(vec &pi, mat &A, vec &f1, vec &f2, vec &f3){
  int m = f1.size();
  vec f0 = ones(m,1);
  mat f = join_rows(join_rows(f0, f3, f2, f2%f3), join_rows(f1, f1%f3, f1%f2, f1%f2%f3));
  
  mat alpha = zeros(m, 8), beta = zeros(m, 8), f_jp1(8,8);
  
  alpha.row(0) = pi.t() % f.row(0);
  alpha.row(0) = alpha.row(0)/sum(alpha.row(0));
  
  for (int j = 1; j < m; j++){
    alpha(j,0) = sum(alpha.row(j-1).t() % A.col(0)) * f(j,0);
    alpha(j,1) = sum(alpha.row(j-1).t() % A.col(1)) * f(j,1);
    alpha(j,2) = sum(alpha.row(j-1).t() % A.col(2)) * f(j,2);
    alpha(j,3) = sum(alpha.row(j-1).t() % A.col(3)) * f(j,3);
    alpha(j,4) = sum(alpha.row(j-1).t() % A.col(4)) * f(j,4);
    alpha(j,5) = sum(alpha.row(j-1).t() % A.col(5)) * f(j,5);
    alpha(j,6) = sum(alpha.row(j-1).t() % A.col(6)) * f(j,6);
    alpha(j,7) = sum(alpha.row(j-1).t() % A.col(7)) * f(j,7);
    
    alpha.row(j) = alpha.row(j)/sum(alpha.row(j));
  }
  
  beta.row(m-1).fill(0.125);
  for(int j = m-2; j >=0; j--){
    beta(j,0) = sum(beta.row(j+1) % A.row(0) % f.row(j+1));
    beta(j,1) = sum(beta.row(j+1) % A.row(1) % f.row(j+1));
    beta(j,2) = sum(beta.row(j+1) % A.row(2) % f.row(j+1));
    beta(j,3) = sum(beta.row(j+1) % A.row(3) % f.row(j+1));
    beta(j,4) = sum(beta.row(j+1) % A.row(4) % f.row(j+1));
    beta(j,5) = sum(beta.row(j+1) % A.row(5) % f.row(j+1));
    beta(j,6) = sum(beta.row(j+1) % A.row(6) % f.row(j+1));
    beta(j,7) = sum(beta.row(j+1) % A.row(7) % f.row(j+1));
    
    beta.row(j) = beta.row(j)/sum(beta.row(j));
  }
  
  vec carLIS(m);
  for(int j = 0; j < m; j++){
    carLIS(j) = sum(alpha.row(j).head(5) % beta.row(j).head(5))/
      sum(alpha.row(j) % beta.row(j));
  }
  
  return carLIS;
}

vec fdr(vec &carLIS){
  int m = carLIS.size();
  
  vec ordered_lis = sort(carLIS), s = linspace(1,m,m);
  uvec ix_lis = sort_index(carLIS);
  
  vec caradj = cumsum(ordered_lis)/s;
  caradj(ix_lis) = caradj;
  
  return caradj;
}

vec my_pava(vec &values, vec &weight, bool decreasing)
{
  if(decreasing){
    values = reverse(values);
    weight = reverse(weight);
  }
  vec w(values.size(), fill::zeros);
  vec x(values.size(), fill::zeros);
  x[0] = values[0];
  w[0] = weight[0];
  unsigned j = 0;
  vec s(values.size(), fill::zeros);
  
  for (unsigned i = 1; i < values.size(); i++) {
    j += 1;
    x[j] = values[i];
    w[j] = weight[i];
    while (j > 0 && x[j - 1]>x[j]) {
      x[j - 1] = (w[j] * x[j] + w[j - 1] * x[j - 1]) / (w[j] + w[j - 1]);
      w[j - 1] += w[j];
      j -= 1;
    }
    s[j + 1] = i + 1;
  }
  
  vec ww(values.size(), fill::zeros);
  vec xx(values.size(), fill::zeros);
  for (unsigned k = 0; k < j + 1; k++) {
    for (unsigned i = s[k]; i < s[k + 1]; i++) {
      ww[i] = w[k];
      xx[i] = x[k];
    }
  }
  
  if(decreasing){
    xx = reverse(xx);
  }
  
  return xx;
}
