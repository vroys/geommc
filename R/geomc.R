#' Markov chain Monte Carlo for discrete and continuous
#' distributions using geometric MH algorithms.
#' @rdname geomc
#' @importFrom stats rnorm runif optim
#' @importFrom cubature divonne
#' @importFrom magrittr %>%
#' @description geomc produces Markov chain samples from a target distribution.
#' The target distribution may be either a probability density function (pdf) or a probability mass function (pmf).
#' Users specify the target by providing an R function that evaluates the log un-normalized pdf or pmf.
#' The function geomc implements the geometric approach of Roy (2024) to move an initial (possibly uninformed)
#'  base density toward one or more approximate target densities, thereby constructing efficient proposal distributions
#'   for MCMC.
#'   The base density can be user-specified through its mean, covariance matrix, density function, and sampling function.
#'  When the base density is Gaussian, it may be specified using only its mean and covariance; however, providing the 
#'  density and sampling functions is recommended when available.
#'     If either or both of the mean and variance arguments and any of the density and sampling functions is
#' omitted, geomc automatically constructs a base density corresponding to a random-walk proposal with an appropriate scale.
#'  One or more approximate target densities
#' can be supplied. Just like the base density, each approximate target may be specified through its mean, covariance, 
#' density, and sampling functions.
#' Gaussian approximate targets may be specified using only their means and covariance matrices, although it is preferred to
#'  supply their densities and sampling functions as well. If either or both of the mean and variance
#' arguments and any of the density and sampling functions is missing for the approximate target density, then geomc 
#' automatically constructs a diffuse multivariate normal distribution as the approximate target.
#' If the argument gaus is set to FALSE, then both the base and approximate target can be specified through their
#' density functions and sampling functions. In this case, any user-supplied mean or covariance functions for these distributions are ignored.
#'  Conversely, if for either the base or the approximate target, the user specifies both a density function and a sampling 
#'  function but omits either the mean function or the covariance function, then gaus = FALSE is automatically enforced.
#' If the target distribution is a pmf (i.e., a discrete distribution), 
#' then gaus=FALSE and imp$enabled=TRUE (not the default values) need to be specified. This ensures that the algorithm uses 
#' importance sampling rather than numerical integration when computing \eqn{\langle \sqrt{f}, \sqrt{g} \rangle} for discrete targets.
#' @param log.target A function that evaluates the logarithm of the 
#'   un-normalized target density (pdf or pmf). Its first argument must be 
#'   the current state vector \code{x}, and it must return a numeric scalar. 
#'   If the target density requires additional parameters, the user-supplied 
#'   \code{log.target} must be written to accept them through \code{...}
#' @param initial is the initial state.
#' @param n.iter is the no. of samples needed.
#' @param eps is the value for epsilon perturbation. Default is 0.5.
#' @param ind is False if either the base density, \eqn{f} or the approximate target density, \eqn{g} depends on
#'  the current state \eqn{x}. Default is False.
#' @param gaus is True if both \eqn{f} and \eqn{g} are normal distributions. Default is True.
#' @param imp A list of three elements controlling optional importance sampling.
#'   This list has components:
#'   \describe{
#'     \item{\code{enabled}}{Logical. If \code{FALSE} (default),
#'       numerical integration is used to compute
#'       \eqn{\langle \sqrt{f}, \sqrt{g} \rangle}. If \code{TRUE},
#'       importance sampling is used instead.}
#'
#'     \item{\code{n.samp}}{A positive integer giving the number of Monte Carlo
#'       samples used when \code{enabled = TRUE}.}
#'
#'     \item{\code{samp.base}}{Logical. If \code{TRUE}, the samples in the
#'       importance sampler are drawn from the base density \eqn{f};
#'       otherwise they are drawn from the target density \eqn{g}.
#'       Default is \code{FALSE}.}
#'   }
#'   When \code{gaus = TRUE}, the \code{imp} argument is ignored.
#' @param a is the probability vector for the mixture proposal density. Default is the uniform distribution.
#' @param mean.base is the mean of the base density \eqn{f}, needs to be written as a function of the current state \eqn{x}.
#' @param var.base is the covariance matrix of the base density \eqn{f}, needs to be written as a function of the current state \eqn{x}.
#' @param dens.base is the density function of the base density \eqn{f}, needs to be written as a function \eqn{(y,x)} (in this order) of the proposed state \eqn{y} and the current state \eqn{x}, although it may not depend on \eqn{x}.
#' @param samp.base is the function to draw from the base density \eqn{f}, needs to be written as a function of the current state \eqn{x}.
#' @param mean.ap.tar is the vector of means of the densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function of the current state \eqn{x}. It must have the same dimension as \eqn{k} times the length of initial.
#' @param var.ap.tar is the matrix of covariance matrices of the densities \eqn{g_i(y|x), i=1,\dots,k} formed by column concatenation. It needs to be written as a function of the current value \eqn{x}. It must have the same dimension as the length of initial by \eqn{k} times the length of initial.
#' @param dens.ap.tar is the vector of densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function \eqn{(y,x)} (in this order) of the proposed state \eqn{y} and the current state \eqn{x}, although it may not depend on \eqn{x}.
#' @param samp.ap.tar is the function to draw from the densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function of (current state \eqn{x}, the indicator of mixing component \eqn{kk}). It must return a vector of the length of that of the initial.
#' @param ... additional arguments passed to \code{log.target}
#' @details
#' Using a geometric Metropolis-Hastings algorithm geom.mc produces Markov chains with the target as its stationary distribution. The details
#' of the method can be found in Roy (2024).
#' @return The function returns a list with the following components:
#' \item{\code{samples}}{A matrix containing the MCMC samples. Each column is one sample.}
#' \item{\code{acceptance.rate}}{The Metropolis-Hastings acceptance rate.}
#' \item{\code{gaus}}{The value of the logical gaus.}
#' \item{\code{ind}}{The value of the logical ind.}
#' \item{\code{model.case}}{An indicator specifying whether both, neither, or which of the functions \eqn{f} and \eqn{g} are missing.}
#' \item{\code{var.base}}{The variance used for the base density if not provided by the user.}
#' \item{\code{mean.ap.tar}}{The mean used for the approximate target density if not provided.}
#' \item{\code{var.ap.tar}}{The variance used for the approximate target density if not provided.}
#' @author Vivekananda Roy <vroy@iastate.edu>
#' @references Roy, V.(2024) A geometric approach to informative MCMC sampling https://arxiv.org/abs/2406.09010
#' @examples
#' log_target_mvnorm <- function(x, target.mean, target.Sigma) {d  <- length(x)
#' xc <- x - target.mean; Q  <- solve(target.Sigma); -0.5 * drop(t(xc) %*% Q %*% xc)}
#' result <- geomc(log.target=log_target_mvnorm,initial= c(0, 0),n.iter=500, target.mean=c(1, -2),
#'                target.Sigma=matrix(c(1.5, 0.7,0.7, 2.0), 2, 2))  
#'                # addidional arguments passed via ...
#' #target is multivariate normal, default choices
#' result$samples # the MCMC samples.
#' result$acceptance.rate # the acceptance rate.
#' #Additional returned values are 
#' result$var.base; result$mean.ap.tar; result$var.ap.tar; result$model.case; result$gaus; result$ind 
#' result<-geomc(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' initial=0,n.iter=500) #target is mixture of univariate normals, default choices
#' hist(result$samples)
#' result<-geomc(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' initial=0,n.iter=500, mean.base = function(x) x,var.base= function(x) 4,
#' dens.base=function(y,x) dnorm(y, mean=x,sd=2),samp.base=function(x) x+2*rnorm(1),
#' mean.ap.tar=function(x) c(0,10),var.ap.tar=function(x) c(1,1.4^2),
#' dens.ap.tar=function(y,x) c(dnorm(y),dnorm(y,mean=10,sd=1.4)),
#' samp.ap.tar=function(x,kk=1){if(kk==1){return(rnorm(1))} else{return(10+1.4*rnorm(1))}})
#' #target is mixture of univariate normals, random walk base density, an informed 
#' #choice for dens.ap.tar
#' hist(result$samples)
#' samp.ap.tar=function(x,kk=1){s.g=sample.int(2,1,prob=c(.5,.5))
#' if(s.g==1){return(rnorm(1))
#' }else{return(10+1.4*rnorm(1))}}
#' result<-geomc(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' initial=0,n.iter=500,gaus=FALSE,imp=list(enabled=TRUE,n.samp=100,samp.base=TRUE),
#' dens.base=function(y,x) dnorm(y, mean=x,sd=2),samp.base=function(x) x+2*rnorm(1),
#' dens.ap.tar=function(y,x) 0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4),
#' samp.ap.tar=samp.ap.tar)
#' #target is mixture of univariate normals, random walk base density, another 
#' #informed choice for dens.ap.tar
#' hist(result$samples)
#' size=5
#' result <- geomc(log.target = function(y) dbinom(y, size, 0.3, log = TRUE),
#' initial=0,n.iter=500,ind=TRUE,gaus=FALSE,imp=list(enabled=TRUE,n.samp=1000,samp.base=TRUE),
#' dens.base=function(y,x) 1/(size+1), samp.base= function(x) sample(seq(0,size,1),1),
#' dens.ap.tar=function(y,x) dbinom(y, size, 0.7),samp.ap.tar=function(x,kk=1) rbinom(1, size, 0.7))
#'  #target is binomial
#'  table(result$samples)
#' @export

geomc=function(log.target,initial,n.iter,eps=0.5,ind=FALSE,gaus=TRUE,imp=list(enabled=FALSE,n.samp=300,samp.base=TRUE),a=1,mean.base,var.base,dens.base,samp.base,mean.ap.tar,var.ap.tar,dens.ap.tar,samp.ap.tar,...){
  if(missing(log.target)) stop("log.target must be provided")

  logtarget.user <- log.target
  
  log.target <- function(x) logtarget.user(x, ...)
  
  if(missing(n.iter)) stop("n.iter must be provided")
  check_positive_integer(n.iter, "n.iter")

  if(log.target(initial)==-Inf || is.na(log.target(initial))) stop("the initial state must satisfy log.target(initial) > -Inf")

  check_fraction_01(eps, "eps")
  if (!is.logical(imp$enabled)) stop("imp$enabled must be logical")
  if (!is.logical(imp$samp.base)) stop("imp$samp.base must be logical")
  if (imp$enabled) check_positive_integer(imp$n.samp, "imp$n.samp")
  
  dd=length(initial)
  
  if(!missing(mean.base)){
  
    if(length(mean.base(initial))!=dd) stop("mean.base must return a vector of same length as initial")
    if(any(is.infinite(mean.base(initial)))||any(is.na(mean.base(initial))))stop("mean.base at initial is not finite")
  }
  if(!missing(var.base)){
    if (!all(dim(var.base(initial)) == c(dd, dd)))
      stop("var.base(x) must return a dd x dd covariance matrix.")
    
    if((dd==1 && var.base(initial)<0) ||(dd>1 &&!is_pd(var.base(initial)))){
      stop("var.base must be a positive definite matrix of dimension length(initial) times length(initial)")
    }
  }
  if(!missing(dens.base)){
    if(is.infinite(dens.base(initial,initial))||is.na(dens.base(initial,initial)))stop("dens.base at initial is not finite")
  }

  if(!missing(samp.base)){
    if(any(is.infinite(samp.base(initial)))||any(is.na(samp.base(initial))))stop("samp.base at initial does not return finite values")
  }

  if(!missing(dens.ap.tar)){
    if(any(is.infinite(dens.ap.tar(initial,initial)))||any(is.na(dens.ap.tar(initial,initial))))stop("dens.ap.tar at initial is not finite")
  }

  if(!missing(mean.ap.tar)){
    if(any(is.infinite(mean.ap.tar(initial)))||any(is.na(mean.ap.tar(initial))))stop("mean.ap.tar at initial is not finite")
    k=length(mean.ap.tar(initial))/dd
    if(!missing(var.ap.tar)){
      if(!checkFuncArgs(var.ap.tar, "x")) stop("var.ap.tar must be a function with only one argument \"x\"")

      if (!all(dim(var.ap.tar(initial)) == c(k * dd, k * dd)))
        stop("var.ap.tar(x) must return a (k*dd) x (k*dd) covariance matrix.")
  
      if((dd==1 && any(var.ap.tar(initial)<0)) ||(dd>1 && any(sapply(seq(1,k),function(ii) !is_pd(var.ap.tar(initial)[,((ii-1)*dd+1):(ii*dd)]))))){
        stop("var.ap.tar must be a set of positive definite matrices each of dimension length(initial) times length(initial) combined by columns")
      }
    }
    if(!missing(dens.ap.tar)){
      if(k!=length(dens.ap.tar(initial)))stop("dens.ap.tar and mean.ap.tar must return the means and densities of the same number of distributions")
    }
  }else if(!missing(dens.ap.tar)){
    k=length(dens.ap.tar(initial))
  }else{
    k=1
  }

  if(!missing(samp.ap.tar)){
    for(ii in 1:k){if(any(is.infinite(samp.ap.tar(initial,kk=ii)))||any(is.na(samp.ap.tar(initial,kk=ii))))stop("samp.ap.tar at initial does not return finite values")
    if(dd!=length(samp.ap.tar(initial,kk=ii))) stop("samp.ap.tar must return a sample of the same size as initial")
    }
  }

  if(imp$enabled){
    if(imp$samp.base){
      if(missing(samp.base)) stop("samp.base must be provided to implement importance sampling with samp.base=TRUE")
    }else{
      if(missing(samp.ap.tar)) stop("samp.ap.tar must be provided to implement importance sampling with samp.base=FALSE")
    }
  }
 
  if (((missing(mean.base) || missing(var.base)) &&
       !missing(dens.base) && !missing(samp.base)) || ((missing(mean.ap.tar) || missing(var.ap.tar)) &&
                                                       !missing(dens.ap.tar) && !missing(samp.ap.tar))) {
    gaus=FALSE
  }
  
  diag.v.ap=FALSE
  
  base_incomplete <- any(missing(dens.base), missing(samp.base)) &&
    any(missing(mean.base), missing(var.base))
  
  ap_tar_incomplete <- any(missing(dens.ap.tar), missing(samp.ap.tar)) &&
    any(missing(mean.ap.tar), missing(var.ap.tar))
  
  if (base_incomplete && ap_tar_incomplete) {
    model.case <- "both_missing"
  } else if (!base_incomplete && ap_tar_incomplete) {
    model.case <- "only_ap_tar_missing"
  } else if (base_incomplete && !ap_tar_incomplete) {
    model.case <- "only_base_missing"
  } else {
    model.case <- "none_missing"
  }
  
  if (model.case == "only_base_missing"){
    ffns=genfvar(log.target,initial,dd)
    s.base=ffns$sd.base
    v.base=ffns$var.base
    dens.base=function(y,x) exp(ldens_mvnorm_diag(y, x,v.base))
    samp.base=function(x) x+s.base*rnorm(dd)
    mean.base=function(x) x
    var.base=function(x) v.base*diag(dd)
    ind=FALSE
  }else if(ap_tar_incomplete){
    k=1
    fgfns=genfvargmeanvar(log.target,initial,dd)
    s.base=fgfns$sd.base
    v.base=fgfns$var.base
    m.ap.tar=fgfns$mean.ap.tar
    v.ap.tar=fgfns$var.ap.tar
    chol.ap.tar=chol(v.ap.tar)
    if(model.case == "both_missing"){
      ind=FALSE
      gaus=TRUE
      dens.base=function(y,x) exp(ldens_mvnorm_diag(y, x,v.base))
      samp.base=function(x) x+s.base*rnorm(dd)
      mean.base=function(x) x
      var.base=function(x) v.base*diag(dd)
      dens.ap.tar=function(y,x) exp(ldens_mvnormchol(y,m.ap.tar,chol.ap.tar))
      samp.ap.tar=function(x,kk=1) m.ap.tar+ crossprod(chol.ap.tar,rnorm(dd))
      mean.ap.tar=function(x) m.ap.tar
      var.ap.tar=function(x)  v.ap.tar
      diag.v.ap=fgfns$diag.v.ap
    }else if(model.case == "only_ap_tar_missing"){
      dens.ap.tar=function(y,x) exp(ldens_mvnormchol(y,m.ap.tar, chol.ap.tar))
      samp.ap.tar=function(x,kk=1) m.ap.tar+ crossprod(chol.ap.tar,rnorm(dd))
      mean.ap.tar=function(x) m.ap.tar
      var.ap.tar=function(x)  v.ap.tar
    } 
  }
  
  if (!missing(mean.base) && !missing(var.base) &&
      (missing(dens.base) || missing(samp.base))) {
    dens.base <- function(y, x) {
      exp(ldens_mvnorm(y, mean.base(x), var.base(x)))
    }
    samp.base <- function(x) {
      rmvnorm_chol(mean.base(x), var.base(x))
    }
  }
  
  if (!missing(mean.ap.tar) && !missing(var.ap.tar) &&
      (missing(dens.ap.tar) || missing(samp.ap.tar))) {
    dens.ap.tar <- function(y, x) {
      mu  <- mean.ap.tar(x)
      if (dd == 1) {
        Sig <- matrix(var.ap.tar(x),nrow=dd)
      } else {
        Sig <- var.ap.tar(x)
      }
      if (k == 1) {
        return(exp(ldens_mvnorm(y, mu, Sig)))
      } else {
        return(vapply(
          seq_len(k),
          function(ii) {
            idx <- ((ii - 1) * dd + 1):(ii * dd)
            exp(ldens_mvnorm(y, mu[idx], Sig[, idx, drop = FALSE]))
          },
          numeric(1)
        ))
      }
    }
    samp.ap.tar <- function(x, kk = 1) {
      mu  <- mean.ap.tar(x)
      if (dd == 1) {
        Sig <- matrix(var.ap.tar(x),nrow=dd)
      } else {
        Sig <- var.ap.tar(x)
      }
      if (k == 1) {
        return(rmvnorm_chol(mu, Sig))
      } else {
        idx <- ((kk - 1) * dd + 1):(kk * dd)
        return(rmvnorm_chol(mu[idx], Sig[, idx, drop = FALSE]))
      }
    }
  }
  
    
  if(length(a)!=k){a=rep.int(1/k,k)}
  ctr_accep =0
  samps = matrix(NA,nrow=dd,ncol=n.iter)
  curr = initial
  log.tar_curr=log.target(curr)
  eps = eps
  a=a/sum(a)

  calc  <- compute_prod_theta(curr, k, gaus, imp,
                              mean.base, var.base, dens.base, samp.base,
                              mean.ap.tar, var.ap.tar, dens.ap.tar, samp.ap.tar,
                              diag.v.ap, dd)
  
  if (ind) {
    prod  <- calc[, 1]
    theta <- calc[, 2]
  } else {
    prod_c  <- calc[, 1]
    theta_c <- calc[, 2]
  }
  
  for (i in 1:n.iter) {
    if (ind) {
      proposed <- as.vector(samp_phi(
        a, curr, prod, theta, eps,
        dens.base, dens.ap.tar, samp.base, samp.ap.tar
      ))
    } else {
      proposed <- as.vector(samp_phi(
        a, curr, prod_c, theta_c, eps,
        dens.base, dens.ap.tar, samp.base, samp.ap.tar
      ))
    }
    log.tar_prop <- log.target(proposed)
    if (!ind) {
      calc_p   <- compute_prod_theta(proposed, k, gaus, imp,
                                     mean.base, var.base, dens.base, samp.base,
                                     mean.ap.tar, var.ap.tar, dens.ap.tar, samp.ap.tar,
                                     diag.v.ap, dd)
      prod_p   <- calc_p[, 1]
      theta_p  <- calc_p[, 2]
    }
    logr <- if (ind) {
      log.tar_prop +
        log_phi(a, curr,proposed,prod,theta,eps,dens.base,dens.ap.tar) -
        log.tar_curr -
        log_phi(a, proposed,curr,prod,theta,eps,dens.base,dens.ap.tar)
    } else {
      log.tar_prop +
        log_phi(a,curr,proposed,prod_p,theta_p,eps,dens.base,dens.ap.tar) -
        log.tar_curr -
        log_phi(a,proposed,curr,prod_c,theta_c,eps,dens.base,dens.ap.tar)
    }
    
    if (is.na(logr)) logr <- -Inf
    if (logr >= 0 || log(runif(1)) < logr) {
      curr        <- proposed
      log.tar_curr <- log.tar_prop
      ctr_accep    <- ctr_accep + 1
      
      if (!ind) {
        prod_c  <- prod_p
        theta_c <- theta_p
      }
    }
    samps[, i] <- curr
  }
  
    if (model.case == "both_missing") {
    return(list(
      var.base = v.base,
      mean.ap.tar = m.ap.tar,
      var.ap.tar = v.ap.tar,
      samples = samps,
      acceptance.rate = ctr_accep / n.iter,
      model.case = model.case,gaus=gaus,ind=ind
    ))
  } else if (model.case == "only_ap_tar_missing") {
    return(list(
      mean.ap.tar = m.ap.tar,
      var.ap.tar = v.ap.tar,
      samples = samps,
      acceptance.rate = ctr_accep / n.iter,
      model.case = model.case,gaus=gaus,ind=ind
    ))
  }else if (model.case == "only_base_missing") {
    return(list(
      var.base = v.base,
      samples = samps,
      acceptance.rate = ctr_accep / n.iter,
      model.case = model.case,gaus=gaus,ind=ind
    ))
  } else {
    return(list(
      samples = samps,
      acceptance.rate = ctr_accep / n.iter,
      model.case = model.case,gaus=gaus,ind=ind
    ))
  }
}