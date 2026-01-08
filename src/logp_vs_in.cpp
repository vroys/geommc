#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]
double logp_vs_in(IntegerVector model,
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
  
  // Pre-allocate
  arma::mat x_g(X.n_rows, p0);
  arma::vec Xty_sub(p0);
  
  // Extract columns and scale (combined loop)
  for(int j = 0; j < p0; j++) {
    int idx = model[j] - 1;  // Convert to 0-based
    arma::vec col = X.col(idx);
    
    double col_mean = arma::mean(col);
    double col_sd = arma::stddev(col, 0);
    x_g.col(j) = (col - col_mean) / col_sd;
    
    Xty_sub(j) = Xty(idx);
  }
  
  // X'X + lambda*I
  arma::mat xtx = x_g.t() * x_g;
  xtx.diag() += lam;
  
  // Cholesky
  arma::mat R = arma::chol(xtx);
  
  // Solve triangular system
  arma::vec z = arma::solve(arma::trimatl(R.t()), Xty_sub);
  
  // Log probability
  double logp = 0.5 * p0 * std::log(lam) - 
    arma::sum(arma::log(R.diag())) - 
    mult_c * std::log(yty - arma::dot(z, z) + add_c) + 
    p0 * logw;
  
  return logp;
}