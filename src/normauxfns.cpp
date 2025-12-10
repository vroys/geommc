#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Simple Cholesky (upper-triangular), returns false if not PD
bool chol_upper(NumericMatrix &A) {
  int n = A.nrow();
  for (int j = 0; j < n; j++) {
    double sum = A(j, j);
    for (int k = 0; k < j; k++)
      sum -= A(k, j) * A(k, j);
    
    if (sum <= 0.0) return false;
    A(j, j) = std::sqrt(sum);
    
    for (int i = j + 1; i < n; i++) {
      double s = A(j, i);
      for (int k = 0; k < j; k++)
        s -= A(k, j) * A(k, i);
      A(j, i) = s / A(j, j);
    }
  }
  
  // zero lower triangle
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      A(i, j) = 0.0;
  
  return true;
}

// --- Backsolve R^T z = diff ---
NumericVector backsolve_transpose(const NumericMatrix &R, const NumericVector &diff) {
  int d = R.nrow();
  NumericVector z(d);
  for (int i = 0; i < d; i++) z[i] = diff[i];
  
  for (int j = 0; j < d; j++) {
    if (R(j, j) <= 0.0) stop("Cholesky diagonal non-positive");
    z[j] /= R(j, j);
    for (int i = j + 1; i < d; i++)
      z[i] -= z[j] * R(j, i);
  }
  return z;
}

// [[Rcpp::export]]
double ldens_mvnorm_cpp(SEXP y_sexp,
                        SEXP mean_sexp,
                        SEXP Sigma_sexp) {
  
  NumericVector y =
    (Rf_isMatrix(y_sexp) || Rf_length(y_sexp) > 1)
  ? NumericVector(y_sexp)
    : NumericVector::create(as<double>(y_sexp));
  
  NumericVector mu =
    (Rf_isMatrix(mean_sexp) || Rf_length(mean_sexp) > 1)
    ? NumericVector(mean_sexp)
      : NumericVector::create(as<double>(mean_sexp));
  
  int d = y.size();
  
  NumericMatrix Sigma;
  if (Rf_isMatrix(Sigma_sexp)) {
    Sigma = NumericMatrix(Sigma_sexp);
  } else {
    Sigma = NumericMatrix(1, 1);
    Sigma(0, 0) = as<double>(Sigma_sexp);
  }
  
  NumericMatrix R = clone(Sigma);
  if (!chol_upper(R))
    return R_NegInf;
  
  // z = y - mu
  NumericVector z(d);
  for (int i = 0; i < d; i++)
    z[i] = y[i] - mu[i];
  
  // Solve Rᵀ z = diff
  for (int j = 0; j < d; j++) {
    z[j] /= R(j, j);
    for (int i = j + 1; i < d; i++)
      z[i] -= z[j] * R(j, i);
  }
  
  double quad = 0.0, logdet = 0.0;
  for (int i = 0; i < d; i++) {
    quad += z[i] * z[i];
    logdet += std::log(R(i, i));
  }
  
  return -0.5 * d * std::log(2.0 * M_PI)
    - logdet
  - 0.5 * quad;
}

// [[Rcpp::export]]
double ldens_mvnormchol_cpp(SEXP y_sexp,
                            SEXP mean_sexp,
                            SEXP chol_sig_sexp) {
  
  // ---- y ----
  NumericVector y;
  if (Rf_isMatrix(y_sexp) || Rf_length(y_sexp) > 1) {
    y = NumericVector(y_sexp);
  } else {
    y = NumericVector(1);
    y[0] = as<double>(y_sexp);
  }
  
  // ---- mean ----
  NumericVector mean;
  if (Rf_isMatrix(mean_sexp) || Rf_length(mean_sexp) > 1) {
    mean = NumericVector(mean_sexp);
  } else {
    mean = NumericVector(1);
    mean[0] = as<double>(mean_sexp);
  }
  
  int d = y.size();
  
  // ---- chol.sig ----
  if (!Rf_isMatrix(chol_sig_sexp))
    return R_NegInf;
  
  NumericMatrix R(chol_sig_sexp);
  if (R.nrow() != R.ncol())
    return R_NegInf;
  
  // ---- z = y - mean ----
  NumericVector z(d);
  for (int i = 0; i < d; i++)
    z[i] = y[i] - mean[i];
  
  // ---- backsolve(R, z, transpose = TRUE) ----
  // Solves Rᵀ z = (y - mean)
  for (int j = 0; j < d; j++) {
    if (R(j, j) <= 0.0)
      return R_NegInf;  // protects log + division
    
    z[j] /= R(j, j);
    for (int i = j + 1; i < d; i++)
      z[i] -= z[j] * R(j, i);
  }
  
  // ---- logdet + quadratic form ----
  double logdet = 0.0;
  double quad   = 0.0;
  
  for (int i = 0; i < d; i++) {
    logdet += std::log(R(i, i));
    quad   += z[i] * z[i];
  }
  
  return -logdet
  - 0.5 * quad
  - 0.5 * d * std::log(2.0 * M_PI);
}

// [[Rcpp::export]]
NumericVector rmvnorm_cpp(SEXP mu_sexp,
                          SEXP Sigma_sexp) {
  
  // ---- mu ----
  NumericVector mu;
  if (Rf_isMatrix(mu_sexp) || Rf_length(mu_sexp) > 1) {
    mu = NumericVector(mu_sexp);
  } else {
    mu = NumericVector(1);
    mu[0] = as<double>(mu_sexp);
  }
  
  int d = mu.size();
  
  // ---- Sigma ----
  NumericMatrix Sigma;
  bool Sigma_is_scalar = false;
  if (Rf_isMatrix(Sigma_sexp)) {
    Sigma = NumericMatrix(Sigma_sexp);
    if (Sigma.nrow() != Sigma.ncol()) {
      stop("Sigma must be a square matrix");
    }
  } else {
    Sigma_is_scalar = true;
  }
  
  // ---- Generate standard normals ----
  NumericVector z(d);
  for (int i = 0; i < d; i++) {
    z[i] = R::rnorm(0.0, 1.0);
  }
  
  // ---- d = 1 case ----
  if (d == 1) {
    double s = Sigma_is_scalar ? as<double>(Sigma_sexp) : Sigma(0,0);
    NumericVector res(1);
    res[0] = mu[0] + std::sqrt(s) * z[0];
    return res;
  }
  
  // ---- d > 1 case ----
  // Cholesky decomposition (upper triangular)
  NumericMatrix R = Sigma_is_scalar
  ? NumericMatrix(d, d)  // Fill diagonal later
    : clone(Sigma);
  
  if (Sigma_is_scalar) {
    for (int i = 0; i < d; i++)
      R(i,i) = std::sqrt(as<double>(Sigma_sexp));
  } else {
    // Pure C++ Cholesky (upper triangular)
    for (int j = 0; j < d; j++) {
      double sum = R(j,j);
      for (int k = 0; k < j; k++)
        sum -= R(k,j) * R(k,j);
      if (sum <= 0.0)
        stop("Sigma is not positive definite");
      R(j,j) = std::sqrt(sum);
      for (int i = j+1; i < d; i++) {
        double s = R(j,i);
        for (int k = 0; k < j; k++)
          s -= R(k,j) * R(k,i);
        R(j,i) = s / R(j,j);
      }
    }
    // zero lower triangle
    for (int i = 0; i < d; i++)
      for (int j = 0; j < i; j++)
        R(i,j) = 0.0;
  }
  
  // ---- Multiply Rᵀ z (crossprod) ----
  NumericVector res(d, 0.0);
  for (int i = 0; i < d; i++) {
    double sum = 0.0;
    for (int j = 0; j < d; j++)
      sum += R(j,i) * z[j];  // crossprod(R, z)
    res[i] = mu[i] + sum;
  }
  
  return res;
}

// [[Rcpp::export]]
NumericMatrix bc_cpp(SEXP mu1_sexp,
                     SEXP mu2_sexp,
                     SEXP sig1_sexp,
                     SEXP sig2_sexp,
                     bool diag_var) {
  
  // --- mu1 and mu2 ---
  NumericVector mu1 = Rf_isMatrix(mu1_sexp) || Rf_length(mu1_sexp) > 1
  ? NumericVector(mu1_sexp)
    : NumericVector::create(as<double>(mu1_sexp));
  
  NumericVector mu2 = Rf_isMatrix(mu2_sexp) || Rf_length(mu2_sexp) > 1
  ? NumericVector(mu2_sexp)
    : NumericVector::create(as<double>(mu2_sexp));
  
  int d = mu1.size();
  
  double prod = 0.0;
  
  if (d > 1) {
    if (!diag_var) {
      // --- sig1 and sig2 must be matrices ---
      if (!Rf_isMatrix(sig1_sexp) || !Rf_isMatrix(sig2_sexp))
        stop("sig1 and sig2 must be matrices when diag_var = FALSE");
      
      NumericMatrix S1(sig1_sexp), S2(sig2_sexp);
      if (S1.nrow() != S1.ncol() || S2.nrow() != S2.ncol())
        stop("sig1 and sig2 must be square matrices");
      
      // sig = 0.5*(sig1 + sig2)
      NumericMatrix sig(d,d);
      for (int i=0;i<d;i++)
        for (int j=0;j<d;j++)
          sig(i,j) = 0.5*(S1(i,j)+S2(i,j));
      
      NumericMatrix R1 = clone(S1), R2 = clone(S2), R = clone(sig);
      if (!chol_upper(R1) || !chol_upper(R2) || !chol_upper(R))
        stop("One of the covariance matrices is not positive definite");
      
      NumericVector diff(d);
      for (int i=0;i<d;i++) diff[i] = mu1[i]-mu2[i];
      
      NumericVector z = backsolve_transpose(R, diff);
      
      double sum_z2 = 0.0, logdetR=0.0, logdetR1=0.0, logdetR2=0.0;
      for (int i=0;i<d;i++) {
        sum_z2 += z[i]*z[i];
        logdetR += std::log(R(i,i));
        logdetR1 += std::log(R1(i,i));
        logdetR2 += std::log(R2(i,i));
      }
      
      prod = std::exp(-sum_z2/8.0 - (logdetR - 0.5*(logdetR1+logdetR2)));
      
    } else {
      // diag.var = TRUE
      double sig1_val, sig2_val;
      if (Rf_isMatrix(sig1_sexp)) {
        NumericMatrix S1(sig1_sexp);
        sig1_val = S1(0,0);
      } else {
        sig1_val = as<double>(sig1_sexp);
      }
      if (Rf_isMatrix(sig2_sexp)) {
        NumericMatrix S2(sig2_sexp);
        sig2_val = S2(0,0);
      } else {
        sig2_val = as<double>(sig2_sexp);
      }
      
      double sig = 0.5*(sig1_val + sig2_val);
      
      double sum_diff2 = 0.0;
      for (int i=0; i<d; i++)
        sum_diff2 += (mu1[i]-mu2[i])*(mu1[i]-mu2[i]);
      
      prod = std::exp(-sum_diff2/(8.0*sig) 
                        - 0.5*d*(std::log(sig) - 0.5*(std::log(sig1_val) + std::log(sig2_val))));
          }
    
  } else {
    // d = 1
    double diff = mu1[0]-mu2[0];
    double s1 = as<double>(sig1_sexp);
    double s2 = as<double>(sig2_sexp);
    double sig = s1+s2;
    prod = std::exp(-diff*diff/(4.0*sig) - 0.5*(std::log(sig)-std::log(2.0*std::sqrt(s1*s2))));
  }
  
  // --- acos bounded ---
  double prod_clamped = std::min(std::max(prod, std::numeric_limits<double>::epsilon()), 1.0-1e-16);
  NumericMatrix out(1,2);
  out(0,0) = prod;
  out(0,1) = std::acos(prod_clamped);
  
  return out;
}