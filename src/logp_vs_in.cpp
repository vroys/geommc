#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Helper function
inline double logp_vs_core(const arma::mat& x_g,
                           const arma::vec& Xty_sub,
                           int p0,
                           double yty,
                           double mult_c,
                           double add_c,
                           double lam,
                           double logw) {
  
  if(p0 == 0) {
    return -mult_c * std::log(add_c + yty);
  }
  
  arma::mat xtx = x_g.t() * x_g;
  xtx.diag() += lam;
  
  arma::mat R = arma::chol(xtx);
  arma::vec z = arma::solve(arma::trimatl(R.t()), Xty_sub);
  
  double logp = 0.5 * p0 * std::log(lam) - 
    arma::sum(arma::log(R.diag())) - 
    mult_c * std::log(yty - arma::dot(z, z) + add_c) + 
    p0 * logw;
  
  return logp;
}

// Version for DENSE matrices
// [[Rcpp::export]]
double logp_vs_in_dense(IntegerVector model,
                        const arma::mat& X,
                        double yty,
                        const arma::vec& Xty,
                        double mult_c,
                        double add_c,
                        double lam,
                        double logw) {
  
  int p0 = model.size();
  
  if(p0 == 0) {
    return -mult_c * std::log(add_c + yty);
  }
  
  arma::mat x_g(X.n_rows, p0);
  arma::vec Xty_sub(p0);
  
  for(int j = 0; j < p0; j++) {
    int idx = model[j] - 1;
    arma::vec col = X.col(idx);
    
    double col_mean = arma::mean(col);
    double col_sd = arma::stddev(col, 0);
    x_g.col(j) = (col - col_mean) / col_sd;
    
    Xty_sub(j) = Xty(idx);
  }
  
  return logp_vs_core(x_g, Xty_sub, p0, yty, mult_c, add_c, lam, logw);
}

// Version for SPARSE matrices
// [[Rcpp::export]]
double logp_vs_in_sparse(IntegerVector model,
                         const arma::sp_mat& X,
                         double yty,
                         const arma::vec& Xty,
                         double mult_c,
                         double add_c,
                         double lam,
                         double logw) {
  
  int p0 = model.size();
  
  if(p0 == 0) {
    return -mult_c * std::log(add_c + yty);
  }
  
  arma::mat x_g(X.n_rows, p0);
  arma::vec Xty_sub(p0);
  
  for(int j = 0; j < p0; j++) {
    int idx = model[j] - 1;
    arma::vec col = arma::vec(X.col(idx));  // Convert sparse to dense
    
    double col_mean = arma::mean(col);
    double col_sd = arma::stddev(col, 0);
    x_g.col(j) = (col - col_mean) / col_sd;
    
    Xty_sub(j) = Xty(idx);
  }
  
  return logp_vs_core(x_g, Xty_sub, p0, yty, mult_c, add_c, lam, logw);
}