#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rwm(Function log_target,
               NumericVector initial,
               int n_iter,
               SEXP sig,
               bool return_sample = false)
{
  RNGScope scope;   // ensures RNG is initialized
  
  const int d = initial.size();
  NumericVector curr = clone(initial);
  NumericVector prop(d);
  
  //Build cholesky(sig)
  //If sig is a scalar, the covariance matrix is assumed to be sig * I_d.
  NumericMatrix chol_sig(d, d);
  
  if (Rf_isMatrix(sig)) {
    NumericMatrix Sig(sig);
    
    // Cholesky
    NumericMatrix A = clone(Sig);
    for (int j = 0; j < d; j++) {
      double s = 0.0;
      for (int k = 0; k < j; k++) s += A(j,k)*A(j,k);
      A(j,j) = std::sqrt(A(j,j) - s);
      for (int i = j+1; i < d; i++) {
        double s2 = 0.0;
        for (int k = 0; k < j; k++) s2 += A(i,k)*A(j,k);
        A(i,j) = (A(i,j) - s2) / A(j,j);
      }
    }
    chol_sig = A;
    
  } else if (Rf_isReal(sig)) {
    double s = REAL(sig)[0];
    double sd = std::sqrt(s);
    for (int i = 0; i < d; i++) chol_sig(i,i) = sd;
    
  } else {
    stop("sig must be scalar or matrix");
  }
  
  double log_tar_curr = as<double>(log_target(curr));
  int ctr_accep = 0;
  NumericVector mean_x(d);
  
  NumericMatrix sample_store;
  if (return_sample) {
    sample_store = NumericMatrix(n_iter,d);
  }
  
  for (int i = 0; i < n_iter; i++) {
    
    // prop = curr + chol_sig %*% rnorm(d)
    NumericVector z(d);
    for (int j = 0; j < d; j++) z[j] = R::rnorm(0,1);
    
    for (int j = 0; j < d; j++) {
      double sum = 0.0;
      for (int k = 0; k < d; k++)
        sum += chol_sig(j,k) * z[k];
      prop[j] = curr[j] + sum;
    }
    
    double log_tar_prop = as<double>(log_target(prop));
    double logr = log_tar_prop - log_tar_curr;
    
    if (logr >= 0 || std::log(R::runif(0,1)) < logr) {
      curr = clone(prop);
      log_tar_curr = log_tar_prop;
      ctr_accep++;
    }
    
    if (return_sample) {
      for (int j = 0; j < d; j++)
        sample_store(i,j) = curr[j];
    }
    
    double w = 1.0 / (i+1);
    for (int j = 0; j < d; j++)
      mean_x[j] += w * (curr[j] - mean_x[j]);
  }
  
  //SAFE retrieval of .Random.seed
  
  // Environment base_env = Environment::namespace_env("base");
  // 
  // // Ensure RNG seed exists (force RNG initialization)
  // if (!base_env.exists(".Random.seed")) {
  //   R::runif(0.0, 1.0);  // forces creation of .Random.seed
  // }
  // 
  // IntegerVector seed = base_env.exists(".Random.seed")
  //   ? as<IntegerVector>(base_env[".Random.seed"])
  //     : IntegerVector();   // fallback (should not happen)
  if (return_sample) {
    return List::create(
      _["mean"] = mean_x,
      _["samples"] = sample_store,
      _["acceptance_rate"] = (double)ctr_accep / n_iter,
      _["n_iter"] = n_iter
    );
  } else {
  return List::create(
    _["mean"] = mean_x,
    _["last_state"] = curr,
    _["acceptance_rate"] = (double)ctr_accep / n_iter,
    _["n_iter"] = n_iter
  //,_["rng_seed"] = seed
  );
}
}