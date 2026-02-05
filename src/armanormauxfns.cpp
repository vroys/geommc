#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @noRd
// [[Rcpp::export]]
double ldens_mvnorm(SEXP y_sexp, SEXP mean_sexp, SEXP Sigma_sexp) {
  arma::vec y;
  if (Rf_isMatrix(y_sexp) || Rf_length(y_sexp) > 1) {
    y = Rcpp::as<arma::vec>(y_sexp);
  } else {
    y = arma::vec(1);
    y(0) = Rcpp::as<double>(y_sexp);
  }
  
  arma::vec mean;
  if (Rf_isMatrix(mean_sexp) || Rf_length(mean_sexp) > 1) {
    mean = Rcpp::as<arma::vec>(mean_sexp);
  } else {
    mean = arma::vec(1);
    mean(0) = Rcpp::as<double>(mean_sexp);
  }
  
  arma::mat Sigma;
  if (Rf_isMatrix(Sigma_sexp)) {
    Sigma = Rcpp::as<arma::mat>(Sigma_sexp);
  } else {
    Sigma = arma::mat(1, 1);
    Sigma(0, 0) = Rcpp::as<double>(Sigma_sexp);
  }
  
  int d = y.n_elem;
  
  if (mean.n_elem != d) {
    Rcpp::stop("Dimension mismatch: y and mean must have the same length");
  }
  if (Sigma.n_rows != d || Sigma.n_cols != d) {
    Rcpp::stop("Dimension mismatch: Sigma must be d x d matrix");
  }
  
  arma::vec diff = y - mean;
  double logdet, quad_form;
  
  bool is_diagonal = arma::approx_equal(Sigma, arma::diagmat(Sigma), "absdiff", 1e-16);
  
  if (is_diagonal) {
    logdet = 0.0;
    quad_form = 0.0;
    
    for (int i = 0; i < d; ++i) {
      double var = Sigma(i, i);
      
      if (var <= 0.0) {
        return -arma::datum::inf;
      }
      
      logdet += std::log(var);
      quad_form += (diff(i) * diff(i)) / var;
    }
   logdet=0.5*logdet; 
  } else {
    arma::mat R;
    bool success = arma::chol(R, Sigma);
    
    if (!success) {
      return -arma::datum::inf;
    }
    
    // Solve R^T * z = diff
    arma::vec z = arma::solve(arma::trimatl(R.t()), diff);
    
    logdet = arma::sum(arma::log(R.diag()));
    
    quad_form = arma::dot(z, z);
  }
  
  double result = -logdet -0.5 * (quad_form + d * std::log(2.0 * M_PI));
  
  return result;
}
 //' Log density of multivariate normal distribution (Cholesky parameterization)
 // [[Rcpp::export]]
 double ldens_mvnormchol(SEXP y_sexp, SEXP mean_sexp, SEXP chol_sig_sexp) {
   // Convert y to vector
   arma::vec y;
   if (Rf_isMatrix(y_sexp) || Rf_length(y_sexp) > 1) {
     y = Rcpp::as<arma::vec>(y_sexp);
   } else {
     y = arma::vec(1);
     y(0) = Rcpp::as<double>(y_sexp);
   }
   
   arma::vec mean;
   if (Rf_isMatrix(mean_sexp) || Rf_length(mean_sexp) > 1) {
     mean = Rcpp::as<arma::vec>(mean_sexp);
   } else {
     mean = arma::vec(1);
     mean(0) = Rcpp::as<double>(mean_sexp);
   }
   
   arma::mat chol_sig;
   if (Rf_isMatrix(chol_sig_sexp)) {
     chol_sig = Rcpp::as<arma::mat>(chol_sig_sexp);
     
     if (chol_sig.n_rows != chol_sig.n_cols) {
       return -arma::datum::inf;
     }
   } else {
     chol_sig = arma::mat(1, 1);
     chol_sig(0, 0) = Rcpp::as<double>(chol_sig_sexp);
   }
   
   int d = y.n_elem;
   
   // Solve R^T * z = (y - mean) where R is upper triangular
   arma::vec diff = y - mean;
   arma::vec z = arma::solve(arma::trimatl(chol_sig.t()), diff);
   
   double logdet = arma::sum(arma::log(chol_sig.diag()));
   
   double result = -logdet - 0.5 * arma::dot(z, z) - 0.5 * d * std::log(2.0 * M_PI);
   
   return result;
 }
 //' Generate random sample from multivariate normal distribution
// [[Rcpp::export]]
arma::vec rmvnorm(SEXP mu_sexp, SEXP Sigma_sexp) {
  arma::vec mu;
  if (Rf_isMatrix(mu_sexp) || Rf_length(mu_sexp) > 1) {
    mu = Rcpp::as<arma::vec>(mu_sexp);
  } else {
    mu = arma::vec(1);
    mu(0) = Rcpp::as<double>(mu_sexp);
  }
  
  int d = mu.n_elem;
  
  arma::vec z = arma::randn<arma::vec>(d);
  
  if (d == 1) {
    double sigma = Rf_isMatrix(Sigma_sexp) ? 
    Rcpp::as<arma::mat>(Sigma_sexp)(0, 0) : 
    Rcpp::as<double>(Sigma_sexp);
    
    if (sigma <= 0.0) {
      Rcpp::stop("Variance must be positive");
    }
    
    mu(0) += std::sqrt(sigma) * z(0);
    return mu;
  }
  
  arma::mat Sigma;
  if (Rf_isMatrix(Sigma_sexp)) {
    Sigma = Rcpp::as<arma::mat>(Sigma_sexp);
  } else {
    Sigma = arma::mat(d, d, arma::fill::eye);
    Sigma *= Rcpp::as<double>(Sigma_sexp);
  }
  
  if (Sigma.n_rows != d || Sigma.n_cols != d) {
    Rcpp::stop("Dimension mismatch: Sigma must be d x d matrix");
  }
  
  bool is_diagonal = arma::approx_equal(Sigma, arma::diagmat(Sigma), "absdiff", 1e-16);
  
  if (is_diagonal) {
    arma::vec result = mu;
    
    for (int i = 0; i < d; ++i) {
      double var = Sigma(i, i);
      
      if (var <= 0.0) {
        Rcpp::stop("All variances must be positive");
      }
      
      result(i) += std::sqrt(var) * z(i);
    }
    
    return result;
    
  } else {
    arma::mat R;
    bool success = arma::chol(R, Sigma);
    
    if (!success) {
      Rcpp::stop("Sigma must be positive definite");
    }
    
    return mu + R.t() * z;
  }
}

 //' Bhattacharya coefficient
 // [[Rcpp::export]]
 arma::mat bc(SEXP mu1_sexp, SEXP mu2_sexp, SEXP sig1_sexp, SEXP sig2_sexp, bool diag_var) {
   arma::vec mu1;
   if (Rf_isMatrix(mu1_sexp) || Rf_length(mu1_sexp) > 1) {
     mu1 = Rcpp::as<arma::vec>(mu1_sexp);
   } else {
     mu1 = arma::vec(1);
     mu1(0) = Rcpp::as<double>(mu1_sexp);
   }
   
   arma::vec mu2;
   if (Rf_isMatrix(mu2_sexp) || Rf_length(mu2_sexp) > 1) {
     mu2 = Rcpp::as<arma::vec>(mu2_sexp);
   } else {
     mu2 = arma::vec(1);
     mu2(0) = Rcpp::as<double>(mu2_sexp);
   }
   
   int ddd = mu1.n_elem;
   double prod;
   
   if (ddd > 1) {
     if (!diag_var) {
       arma::mat sig1 = Rcpp::as<arma::mat>(sig1_sexp);
       arma::mat sig2 = Rcpp::as<arma::mat>(sig2_sexp);
       arma::mat sig = 0.5 * (sig1 + sig2);
       
       arma::mat R1 = arma::chol(sig1);
       arma::mat R2 = arma::chol(sig2);
       arma::mat R = arma::chol(sig);
       
       arma::vec diff = mu1 - mu2;
       arma::vec z = arma::solve(arma::trimatl(R.t()), diff);
       
       double sum_z2 = arma::dot(z, z);
       double log_det_R = arma::sum(arma::log(R.diag()));
       double log_det_R1 = arma::sum(arma::log(R1.diag()));
       double log_det_R2 = arma::sum(arma::log(R2.diag()));
       
       prod = std::exp(-sum_z2 / 8.0 - (log_det_R - 0.5 * (log_det_R1 + log_det_R2)));
     } else {
       arma::mat sig1_mat = Rcpp::as<arma::mat>(sig1_sexp);
       arma::mat sig2_mat = Rcpp::as<arma::mat>(sig2_sexp);
       
       arma::vec sig1_diag = sig1_mat.diag();
       arma::vec sig2_diag = sig2_mat.diag();
       arma::vec sig = 0.5 * (sig1_diag + sig2_diag);
       
       arma::vec diff = mu1 - mu2;
       arma::vec diff_sq = arma::square(diff);
       
       double sum_term = arma::sum(diff_sq / sig);
       
       double sum_log_sig = arma::sum(arma::log(sig));
       double sum_log_sig1_sig2 = arma::sum(arma::log(sig1_diag) + arma::log(sig2_diag));
       
       prod = std::exp(-sum_term / 8.0 - 0.5 * (sum_log_sig - 0.5 * sum_log_sig1_sig2));
     }
   } else {
     double sig1_val = Rf_isMatrix(sig1_sexp) ? 
     Rcpp::as<arma::mat>(sig1_sexp)(0, 0) : 
     Rcpp::as<double>(sig1_sexp);
     double sig2_val = Rf_isMatrix(sig2_sexp) ? 
     Rcpp::as<arma::mat>(sig2_sexp)(0, 0) : 
       Rcpp::as<double>(sig2_sexp);
     
     double diff = mu1(0) - mu2(0);
     double sig_sum = sig1_val + sig2_val;
     
     prod = std::exp(-diff * diff / (4.0 * sig_sum) - 
       0.5 * (std::log(sig_sum) - std::log(2.0 * std::sqrt(sig1_val * sig2_val))));
   }
   
   arma::mat result(1, 2);
   result(0, 0) = prod;
   
   double bounded_prod = std::max(std::numeric_limits<double>::epsilon(),
                                  std::min(prod, 1.0 - 1e-16));
   result(0, 1) = std::acos(bounded_prod);
   
   return result;
 }