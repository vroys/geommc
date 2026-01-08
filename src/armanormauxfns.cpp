#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @noRd
// [[Rcpp::export]]
 double ldens_mvnorm(SEXP y_sexp, SEXP mean_sexp, SEXP Sigma_sexp) {
   // Convert y to vector
   arma::vec y;
   if (Rf_isMatrix(y_sexp) || Rf_length(y_sexp) > 1) {
     y = Rcpp::as<arma::vec>(y_sexp);
   } else {
     y = arma::vec(1);
     y(0) = Rcpp::as<double>(y_sexp);
   }
   
   // Convert mean to vector
   arma::vec mean;
   if (Rf_isMatrix(mean_sexp) || Rf_length(mean_sexp) > 1) {
     mean = Rcpp::as<arma::vec>(mean_sexp);
   } else {
     mean = arma::vec(1);
     mean(0) = Rcpp::as<double>(mean_sexp);
   }
   
   // Convert Sigma to matrix
   arma::mat Sigma;
   if (Rf_isMatrix(Sigma_sexp)) {
     Sigma = Rcpp::as<arma::mat>(Sigma_sexp);
   } else {
     Sigma = arma::mat(1, 1);
     Sigma(0, 0) = Rcpp::as<double>(Sigma_sexp);
   }
   
   int d = y.n_elem;
   
   // Cholesky decomposition
   arma::mat R;
   bool success = arma::chol(R, Sigma);
   
   if (!success) {
     return -arma::datum::inf;
   }
   
   // Solve R^T * z = (y - mean)
   arma::vec diff = y - mean;
   arma::vec z = arma::solve(arma::trimatl(R.t()), diff);
   
   // Log determinant
   double logdet = arma::sum(arma::log(R.diag()));
   
   // Compute log density
   double result = -logdet - 0.5 * arma::dot(z, z) - 0.5 * d * std::log(2.0 * M_PI);
   
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
   
   // Convert mean to vector
   arma::vec mean;
   if (Rf_isMatrix(mean_sexp) || Rf_length(mean_sexp) > 1) {
     mean = Rcpp::as<arma::vec>(mean_sexp);
   } else {
     mean = arma::vec(1);
     mean(0) = Rcpp::as<double>(mean_sexp);
   }
   
   // Convert chol.sig to matrix
   arma::mat chol_sig;
   if (Rf_isMatrix(chol_sig_sexp)) {
     chol_sig = Rcpp::as<arma::mat>(chol_sig_sexp);
     
     // Check if square matrix
     if (chol_sig.n_rows != chol_sig.n_cols) {
       return -arma::datum::inf;
     }
   } else {
     // Scalar case: create 1x1 matrix
     chol_sig = arma::mat(1, 1);
     chol_sig(0, 0) = Rcpp::as<double>(chol_sig_sexp);
   }
   
   int d = y.n_elem;
   
   // Solve R^T * z = (y - mean) where R is upper triangular
   arma::vec diff = y - mean;
   arma::vec z = arma::solve(arma::trimatl(chol_sig.t()), diff);
   
   // Log determinant: sum of log of diagonal elements
   double logdet = arma::sum(arma::log(chol_sig.diag()));
   
   // Compute log density
   double result = -logdet - 0.5 * arma::dot(z, z) - 0.5 * d * std::log(2.0 * M_PI);
   
   return result;
 }
 //' Generate random sample from multivariate normal distribution
 // [[Rcpp::export]]
 arma::vec rmvnorm(SEXP mu_sexp, SEXP Sigma_sexp) {
   // Convert mu to vector
   arma::vec mu;
   if (Rf_isMatrix(mu_sexp) || Rf_length(mu_sexp) > 1) {
     mu = Rcpp::as<arma::vec>(mu_sexp);
   } else {
     mu = arma::vec(1);
     mu(0) = Rcpp::as<double>(mu_sexp);
   }
   
   int d = mu.n_elem;
   
   // Generate standard normal random variables
   arma::vec z = arma::randn<arma::vec>(d);
   
   // Handle scalar case
   if (d == 1) {
     double sigma = Rf_isMatrix(Sigma_sexp) ? 
     Rcpp::as<arma::mat>(Sigma_sexp)(0, 0) : 
     Rcpp::as<double>(Sigma_sexp);
     mu(0) += std::sqrt(sigma) * z(0);
     return mu;
   }
   
   // Multivariate case
   arma::mat Sigma;
   if (Rf_isMatrix(Sigma_sexp)) {
     Sigma = Rcpp::as<arma::mat>(Sigma_sexp);
   } else {
     // Should not happen for d > 1, but handle gracefully
     Sigma = arma::mat(d, d, arma::fill::eye);
     Sigma *= Rcpp::as<double>(Sigma_sexp);
   }
   
   // Cholesky decomposition
   arma::mat R = arma::chol(Sigma);
   
   // Return mu + R^T * z
   return mu + R.t() * z;
 }
 //' Bhattacharya coefficient
 // [[Rcpp::export]]
 arma::mat bc(SEXP mu1_sexp, SEXP mu2_sexp, SEXP sig1_sexp, SEXP sig2_sexp, bool diag_var) {
   // Convert mu1 to vector
   arma::vec mu1;
   if (Rf_isMatrix(mu1_sexp) || Rf_length(mu1_sexp) > 1) {
     mu1 = Rcpp::as<arma::vec>(mu1_sexp);
   } else {
     mu1 = arma::vec(1);
     mu1(0) = Rcpp::as<double>(mu1_sexp);
   }
   
   // Convert mu2 to vector
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
       // Full covariance case
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
       // Diagonal covariance case (multivariate)
       arma::mat sig1_mat = Rcpp::as<arma::mat>(sig1_sexp);
       arma::mat sig2_mat = Rcpp::as<arma::mat>(sig2_sexp);
       
       // Extract diagonal elements
       arma::vec sig1_diag = sig1_mat.diag();
       arma::vec sig2_diag = sig2_mat.diag();
       arma::vec sig = 0.5 * (sig1_diag + sig2_diag);
       
       arma::vec diff = mu1 - mu2;
       arma::vec diff_sq = arma::square(diff);
       
       // Element-wise division: (mu1 - mu2)^2 / sig
       double sum_term = arma::sum(diff_sq / sig);
       
       // Compute log terms
       double sum_log_sig = arma::sum(arma::log(sig));
       double sum_log_sig1_sig2 = arma::sum(arma::log(sig1_diag) + arma::log(sig2_diag));
       
       prod = std::exp(-sum_term / 8.0 - 0.5 * (sum_log_sig - 0.5 * sum_log_sig1_sig2));
     }
   } else {
     // Scalar case
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
   
   // Create result matrix
   arma::mat result(1, 2);
   result(0, 0) = prod;
   
   // Compute acos with bounds
   double bounded_prod = std::max(std::numeric_limits<double>::epsilon(),
                                  std::min(prod, 1.0 - 1e-16));
   result(0, 1) = std::acos(bounded_prod);
   
   return result;
 }