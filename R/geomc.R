#' Markov chain Monte Carlo for discrete and continuous
#' distributions using geometric MH algorithms.
#' @rdname geomc
#' @importFrom stats rnorm runif optim
#' @importFrom cubature divonne
#' @description geomc produces Markov chain samples from a user-defined target, which
#' may be either a probability density function (pdf) or a probability mass function (pmf).
#' The target is provided as an R function returning the log of its unnormalized pdf or pmf.
#' @details
#' The function \code{geomc} implements the geometric approach of Roy (2024) to move an initial (possibly uninformed)
#'  base density \eqn{f} toward one or more approximate target densities \eqn{g_i, i=1,\dots,k}, thereby constructing 
#'  efficient proposal distributions for MCMC. The details of the method can be found in Roy (2024).
#'   
#'   The base density \eqn{f} can be user-specified through its mean, covariance matrix, density function, and sampling 
#'   function. One or more approximate target densities \eqn{g_i, i=1,\dots,k} can also be supplied. Just like the base 
#'   density, each approximate target may be specified through its mean, covariance, density, and sampling functions.
#' 
#'  A Gaussian (\eqn{f} or \eqn{g_i}) density must be specified in terms of its mean
#'  and covariance; optional density and sampling functions may also be supplied.
#' 
#'  Non-Gaussian densities, including discrete pmfs for either the base or approximate target, are specified solely 
#'  via their density and sampling functions. If for either the base or the approximate target, the user specifies 
#'  both a density function and a sampling function but omits either the mean function or the covariance function, 
#'  then \code{gaus = FALSE} is automatically enforced.
#'  
#' If either or both of the mean and variance arguments and any of the density and sampling functions is
#' omitted for the base density, \code{geomc} automatically constructs a base density corresponding to a 
#' random-walk proposal with an appropriate scale. If either or both of the mean and variance
#' arguments and any of the density and sampling functions is missing for the approximate target density, 
#' then \code{geomc} automatically constructs a diffuse multivariate normal distribution as the approximate target.
#' 
#' If the target distribution is a pmf (i.e., a discrete distribution) then the user must provide the base pmf \eqn{f} 
#' and one or more \eqn{g_i}'s. The package vignette \code{vignette("geommc")} provides several useful choices for \eqn{f}
#' and \eqn{g_i}. In addition, the user must set \code{gaus=FALSE} and either supply \code{bhat.coef} or set \code{imp$enabled=TRUE} 
#' (neither of which is the default for \code{gaus} or \code{imp}). This ensures that the algorithm uses either the user defined 
#' \code{bhat.coef} function or the importance sampling, rather than numerical integration or closed-form Gaussian expressions,
#' for computing the Bhattacharyya coefficient \eqn{\langle \sqrt{f}, \sqrt{g_i} \rangle} for discrete distributions.
#' 
#' For a discussion and several illustrative examples of the \code{geomc} function, see the package vignette
#' \code{vignette("geommc")}.
#' @param logp Either a function that evaluates the logarithm of the un-normalized target density (pdf or pmf) to sample from, or
#' a list containing at least one element named \code{log.target}. The list may optionally include any of the following elements
#' specifying additional functions for the base and approximate target densities: \code{mean.base}, \code{var.base}, \code{dens.base},
#' \code{samp.base}, \code{mean.ap.tar}, \code{var.ap.tar}, \code{dens.ap.tar}, \code{samp.ap.tar}, and \code{bhat.coef}. See below for details of these functions.
#' Any optional elements not provided are treated as missing, and default behavior is applied.
#' \itemize{
#' \item \code{log.target} A function that evaluates the logarithm of the
#'   un-normalized target density (pdf or pmf). Its first argument must be
#'   the current state vector \eqn{x}, and it must return a numeric scalar.
#'   If the target density requires additional parameters, the user-supplied
#'   \code{log.target} must be written to accept them through \code{...}
#'   \item \code{mean.base} is the mean of the base density \eqn{f}, needs to be written as a function of the current state \eqn{x}.
#' \item \code{var.base} is the covariance matrix of the base density \eqn{f}, needs to be written as a function of the current state \eqn{x}.
#'  \item \code{dens.base} is the density function of the base density \eqn{f}, needs to be written as a function \eqn{(y,x)} (in this order) of the proposed state \eqn{y} and the current state \eqn{x}, although it may not depend on \eqn{x}.
#'  \item \code{samp.base} is the function to draw from the base density \eqn{f}, needs to be written as a function of the current state \eqn{x}.
#'  \item \code{mean.ap.tar} is the vector of means of the densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function of the current state \eqn{x}. It must have the same dimension as \eqn{k} times the length of initial.
#'  \item \code{var.ap.tar} is the matrix of covariance matrices of the densities \eqn{g_i(y|x), i=1,\dots,k} formed by column concatenation. It needs to be written as a function of the current state \eqn{x}. It must have the same dimension as the length of initial by \eqn{k} times the length of initial.
#'  \item \code{dens.ap.tar} is the vector of densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function \eqn{(y,x)} (in this order) of the proposed state \eqn{y} and the current state \eqn{x}, although it may not depend on \eqn{x}.
#'  \item \code{samp.ap.tar} is the function to draw from the densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function of (current state \eqn{x}, the indicator of component \eqn{kk}). It must return a vector of the length of that of the initial.
#'  \item \code{bhat.coef} is the vector of Bhattacharyya coefficients \eqn{\langle \sqrt{f(\cdot|x)}, \sqrt{g_i(\cdot|x)} \rangle, i=1,\dots,k}. It needs to be written as a function of the current state \eqn{x}.
#' When \code{gaus= TRUE}, the \code{bhat.coef} argument is ignored and the built-in function is used.
#' }
#' @param initial is the initial state.
#' @param n.iter is the no. of samples needed.
#' @param eps is the value for epsilon perturbation. Default is 0.5.
#' @param ind is False if either the base density, \eqn{f} or the approximate target density, \eqn{g} depends on
#'  the current state \eqn{x}. Default is False.
#' @param gaus is True if both \eqn{f} and \eqn{g} are normal distributions. Default is True.
#' @param imp A list of three elements controlling optional importance sampling.
#'   This list has components:
#'   \itemize{
#'     \item \code{enabled} Logical. If \code{FALSE} (default),
#'       numerical integration is used to compute the Bhattacharyya coefficient
#'       \eqn{\langle \sqrt{f}, \sqrt{g_i} \rangle}. If \code{TRUE},
#'       importance sampling is used instead.
#'     \item \code{n.samp} A positive integer giving the number of Monte Carlo
#'       samples used when \code{enabled = TRUE}.
#'     \item \code{samp.base} Logical. If \code{TRUE}, the samples in the
#'       importance sampler are drawn from the base density \eqn{f};
#'       otherwise they are drawn from the approximate target density \eqn{g}.
#'       Default is \code{FALSE}.
#'   }
#'   When \code{gaus = TRUE}, the \code{imp} argument is ignored. Also, when \code{bhat.coef} is provided, 
#'   the \code{imp} argument is ignored.
#' @param a is the probability vector for the mixture proposal density. Default is the uniform distribution.
#' @param show.progress Logical. Whether to show progress during sampling. Default: TRUE.
#' @param ... additional arguments passed to \code{log.target}
#' @return The function returns a list with the following components:
#' \item{\code{samples}}{A matrix containing the MCMC samples. Each row is one sample.}
#' \item{\code{acceptance.rate}}{The Metropolis-Hastings acceptance rate.}
#' \item{\code{log.p}}{A vector with the logarithm of the (un-normalized) target density for each MCMC sample.}
#' \item{\code{gaus}}{The value of the logical gaus.}
#' \item{\code{ind}}{The value of the logical ind.}
#' \item{\code{model.case}}{Describes whether default settings are used for both functions 
#' \eqn{f} and \eqn{g}, for only one of them, or for neither.}
#' \item{\code{var.base}}{The default variance used for the random-walk base density if not provided by the user.}
#' \item{\code{mean.ap.tar}}{The default mean vector used for the approximate target density if not provided 
#' by the user.}
#' \item{\code{var.ap.tar}}{The default covariance matrix used for the approximate target density if not provided
#'  by the user.}
#' @author Vivekananda Roy <vroy@iastate.edu>
#' @references Roy, V.(2024) A geometric approach to informative MCMC sampling https://arxiv.org/abs/2406.09010
#' @examples
#' #target is multivariate normal, sampling using geomc with default choices
#' log_target_mvnorm <- function(x, target.mean, target.Sigma) {d  <- length(x)
#' xc <- x - target.mean; Q  <- solve(target.Sigma); -0.5 * drop(t(xc) %*% Q %*% xc)}
#' result <- geomc(logp=log_target_mvnorm,initial= c(0, 0),n.iter=500, target.mean=c(1, -2),
#'                target.Sigma=matrix(c(1.5, 0.7,0.7, 2.0), 2, 2))
#'                # additional arguments passed via ...
#' result$samples # the MCMC samples.
#' result$acceptance.rate # the acceptance rate.
#' result$log.p # the value of logp at the MCMC samples.
#' #Additional returned values are
#' result$var.base; result$mean.ap.tar; result$var.ap.tar; result$model.case; result$gaus; result$ind
#' #target is the posterior of (\eqn{\mu}, \eqn{\sigma^2}) with iid data from 
#' #N(\eqn{\mu}, \eqn{\sigma^2}) and N(mu0, \eqn{tau0^2}) prior on \eqn{\mu} 
#' #and an inverse–gamma (alpha0, beta0) prior on \eqn{\sigma^2}
#' log_target <- function(par, x, mu0, tau0, alpha0, beta0) {
#' mu  <- par[1];sigma2 <- par[2]
#' if (sigma2 <= 0) return(-Inf)
#' n  <- length(x); SSE <- sum((x - mu)^2)
#' val <- -(n/2) * log(sigma2) - SSE / (2 * sigma2) - (mu - mu0)^2 / (2 * tau0^2)
#' -(alpha0 + 1) * log(sigma2) -beta0 / sigma2
#' return(val)}
#' # sampling using geomc with default choices
#' result=geomc(logp=log_target,initial=c(0,1),n.iter=1000,x=1+rnorm(100),mu0=0,          
#' tau0=1,alpha0=2.01,beta0=1.01)
#' colMeans(result$samples)
#' #target is mixture of univariate normals, sampling using geomc with default choices
#' set.seed(3);result<-geomc(logp=list(log.target=function(y) 
#' log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4))),
#' initial=0,n.iter=1000) 
#' hist(result$samples)
#' #target is mixture of univariate normals, sampling using geomc with a random walk base density,
#' # and an informed choice for dens.ap.tar
#' result<-geomc(logp=list(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' mean.base = function(x) x,
#' var.base= function(x) 4, dens.base=function(y,x) dnorm(y, mean=x,sd=2),
#' samp.base=function(x) x+2*rnorm(1),
#' mean.ap.tar=function(x) c(0,10),var.ap.tar=function(x) c(1,1.4^2),
#' dens.ap.tar=function(y,x) c(dnorm(y),dnorm(y,mean=10,sd=1.4)),
#' samp.ap.tar=function(x,kk=1){if(kk==1){return(rnorm(1))} else{return(10+1.4*rnorm(1))}}),
#' initial=0,n.iter=500)
#' hist(result$samples)
#' #target is mixture of univariate normals, sampling using geomc with a random walk base density, 
#' # and another informed choice for dens.ap.tar
#' samp.ap.tar=function(x,kk=1){s.g=sample.int(2,1,prob=c(.5,.5))
#' if(s.g==1){return(rnorm(1))
#' }else{return(10+1.4*rnorm(1))}}
#' result<-geomc(logp=list(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' dens.base=function(y,x) dnorm(y, mean=x,sd=2),samp.base=function(x) x+2*rnorm(1),
#' dens.ap.tar=function(y,x) 0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4), samp.ap.tar=samp.ap.tar),
#' initial=0,n.iter=500,gaus=FALSE,imp=list(enabled=TRUE,n.samp=100,samp.base=TRUE))
#' hist(result$samples)
#' #target is mixture of bivariate normals, sampling using geomc with random walk base density,
#' # and an informed choice for dens.ap.tar
#' log_target_mvnorm_mix <- function(x, mean1, Sigma1, mean2, Sigma2) {
#' return(log(0.5*exp(geommc:::ldens_mvnorm(x, mean1,Sigma1))+
#' 0.5*exp(geommc:::ldens_mvnorm(x,mean2,Sigma2))))}
#' result <- geomc(logp=list(log.target=log_target_mvnorm_mix, mean.base = function(x) x, 
#' var.base= function(x) 2*diag(2), mean.ap.tar=function(x) c(0,0,10,10),
#' var.ap.tar=function(x) cbind(diag(2),2*diag(2))),initial= c(5, 5),n.iter=500, mean1=c(0, 0), 
#' Sigma1=diag(2), mean2=c(10, 10), Sigma2=2*diag(2))
#' colMeans(result$samples)
#' #While the geomc with random walk base successfully moves back and forth between the two modes,
#' # the random walk Metropolis is trapped in a mode, as run below 
#' result<-geommc:::rwm(log_target= function(x) log_target_mvnorm_mix(x,c(0, 0), diag(2), 
#' c(10, 10), 2*diag(2)), initial= c(5, 5), n_iter = 500, sig = 2,return_sample = TRUE)
#' colMeans(result$samples)
#'  #target is binomial, sampling using geomc
#' size=5
#' result <- geomc(logp=list(log.target = function(y) dbinom(y, size, 0.3, log = TRUE),
#' dens.base=function(y,x) 1/(size+1), samp.base= function(x) sample(seq(0,size,1),1),
#' dens.ap.tar=function(y,x) dbinom(y, size, 0.7),samp.ap.tar=function(x,kk=1) rbinom(1, size, 0.7)),
#' initial=0,n.iter=500,ind=TRUE,gaus=FALSE,imp=list(enabled=TRUE,n.samp=1000,samp.base=TRUE))
#'  table(result$samples)
#' @export

geomc=function(logp,initial,n.iter,eps=0.5,ind=FALSE,gaus=TRUE,imp=list(enabled=FALSE,n.samp=100,samp.base=TRUE),a=1,show.progress=TRUE,...){
  if (missing(logp)) stop("logp must be provided")
  
  if (is.function(logp)) {
    logtarget.user <- logp
    mean.base <- var.base <- dens.base <- samp.base <- NULL
    mean.ap.tar <- var.ap.tar <- dens.ap.tar <- samp.ap.tar <- NULL
    bhat.coef <- NULL
  } else if (is.list(logp)) {

    if (!("log.target" %in% names(logp)))
      stop("If logp is a list, it must contain an element named 'log.target'.")
    
    if (!is.function(logp$log.target))
      stop("logp$log.target must be a function")

    logtarget.user <- logp$log.target

    mean.base   <- get_opt(logp, "mean.base")
    var.base    <- get_opt(logp, "var.base")
    dens.base   <- get_opt(logp, "dens.base")
    samp.base   <- get_opt(logp, "samp.base")

    mean.ap.tar <- get_opt(logp, "mean.ap.tar")
    var.ap.tar  <- get_opt(logp, "var.ap.tar")
    dens.ap.tar <- get_opt(logp, "dens.ap.tar")
    samp.ap.tar <- get_opt(logp, "samp.ap.tar")
    
    bhat.coef <- get_opt(logp, "bhat.coef")
    
    optional_funcs <- list(
      mean.base = mean.base, var.base = var.base, 
      dens.base = dens.base, samp.base = samp.base,
      mean.ap.tar = mean.ap.tar, var.ap.tar = var.ap.tar,
      dens.ap.tar = dens.ap.tar, samp.ap.tar = samp.ap.tar,
      bhat.coef = bhat.coef
    )
    
    for (fname in names(optional_funcs)) {
      if (!is.null(optional_funcs[[fname]]) && !is.function(optional_funcs[[fname]])) {
        stop(sprintf("logp$%s must be a function if provided", fname))
      }
    }
    } else {
    stop("logp must be either a function or a list.")
  }

  log.target <- function(x) logtarget.user(x, ...)

    if (missing(initial)) {
      stop("The 'initial' argument must be provided.")
    }

  if (!is.numeric(initial)) {
    stop("The 'initial' argument must be a numeric vector.")
  }

  if (length(initial) < 1) {
    stop("The 'initial' vector must have length at least 1.")
  }
  
  if (any(!is.finite(initial))) {
    stop("All elements of 'initial' must be finite (no NA, NaN, or Inf).")
  }
  
  test_eval <- try(log.target(initial), silent = TRUE)
  
  if (inherits(test_eval, "try-error")) {
    stop("Could not evaluate log.target at 'initial'. Please check that log-density function is compatible with the initial values provided.")
  }
  
  if (!is.numeric(test_eval) || length(test_eval) != 1) {
    stop("log.target must return a single numeric value")
  }

  if(missing(n.iter)) stop("n.iter must be provided")
  check_positive_integer(n.iter, "n.iter")

  logp_initial <- log.target(initial)
  if(is.na(logp_initial)) {
    stop("log.target(initial) returns NA")
  }
  if(logp_initial == -Inf) {
    stop("log.target(initial) is -Inf. The initial state must have positive (non-zero) density.")
  }
  if(logp_initial == Inf) {
    stop("log.target(initial) is Inf. Please check your log-density function.")
  }
  
  check_fraction_01(eps, "eps")
  
  if (!is.logical(ind) || length(ind) != 1 || is.na(ind))
    stop("ind must be a single logical value (TRUE or FALSE)")
  
  if (!is.logical(gaus) || length(gaus) != 1 || is.na(gaus))
    stop("gaus must be a single logical value (TRUE or FALSE)")
  
  if (!is.list(imp))
    stop("imp must be a list")
  
  if (!("enabled" %in% names(imp)))
    stop("imp must contain element 'enabled'")
  if (!("n.samp" %in% names(imp)))
    stop("imp must contain element 'n.samp'")
  if (!("samp.base" %in% names(imp)))
    stop("imp must contain element 'samp.base'")
  
  if (!is.logical(imp$enabled) || length(imp$enabled) != 1 || is.na(imp$enabled))
    stop("imp$enabled must be a single logical value")
  if (!is.logical(imp$samp.base) || length(imp$samp.base) != 1 || is.na(imp$samp.base))
    stop("imp$samp.base must be a single logical value")
  if (imp$enabled) check_positive_integer(imp$n.samp, "imp$n.samp")
  
  if (!is.numeric(a) || any(is.na(a)) || any(!is.finite(a)) || any(a < 0))
    stop("a must be a numeric vector of non-negative finite values")
  if (sum(a) == 0)
    stop("a must have at least one positive element")

  a <- a / sum(a)
  
  if (!is.logical(show.progress) || length(show.progress) != 1 || is.na(show.progress))
    stop("show.progress must be a single logical value (TRUE or FALSE)")
  
  dd=length(initial)

  if(!is.null(mean.base)){
    mb <- try(mean.base(initial), silent = TRUE)
    if (inherits(mb, "try-error"))
      stop("mean.base(initial) produced an error")
    if(!is.numeric(mb))
      stop("mean.base must return a numeric vector")
    if(length(mb)!=dd) 
      stop(sprintf("mean.base must return a vector of length %d (same as initial)", dd))
    if(any(is.infinite(mb))||any(is.na(mb)))
      stop("mean.base(initial) is not finite")
  }
  
  if(!is.null(var.base)){
    vb <- try(var.base(initial), silent = TRUE)
    if (inherits(vb, "try-error"))
      stop("var.base(initial) produced an error")
    
    if (dd == 1) {
      if (!is.numeric(vb) || length(vb) != 1 || vb <= 0) {
        stop("var.base(x) must return a positive scalar variance when initial has length 1.")
      }
    } else {
      if (!is.matrix(vb) || any(dim(vb) != c(dd, dd))) {
        stop(sprintf("var.base(x) must return a %d x %d covariance matrix.", dd, dd))
      }
      if (!isSymmetric(vb)) {
        stop("var.base(x) must return a symmetric matrix")
      }
      if (!is_pd(vb)) {
        stop("var.base(x) must be positive definite.")
      }
    }
  }
  
  if(!is.null(dens.base)){
    db <- try(dens.base(initial, initial), silent = TRUE)
    if (inherits(db, "try-error"))
      stop("dens.base(initial, initial) produced an error")
    if(!is.numeric(db) || length(db) != 1)
      stop("dens.base must return a single numeric value")
    if(is.infinite(db)||is.na(db))
      stop("dens.base(initial, initial) is not finite")
    if(db < 0)
      stop("dens.base must return non-negative values (it's a density)")
  }
  
  if(!is.null(samp.base)){
    sb <- try(samp.base(initial), silent = TRUE)
    if (inherits(sb, "try-error"))
      stop("samp.base(initial) produced an error")
    if(!is.numeric(sb))
      stop("samp.base must return a numeric vector")
    if(length(sb) != dd)
      stop(sprintf("samp.base must return a vector of length %d (same as initial)", dd))
    if(any(is.infinite(sb))||any(is.na(sb)))
      stop("samp.base(initial) does not return finite values")
  }
  
  if(!is.null(dens.ap.tar)){
    dap <- try(dens.ap.tar(initial, initial), silent = TRUE)
    if (inherits(dap, "try-error"))
      stop("dens.ap.tar(initial, initial) produced an error")
    if(!is.numeric(dap))
      stop("dens.ap.tar must return a numeric vector")
    if(any(is.infinite(dap))||any(is.na(dap)))
      stop("dens.ap.tar(initial, initial) is not finite")
    if(any(dap < 0))
      stop("dens.ap.tar must return non-negative values (it's a density)")
  }
  
  if(!is.null(mean.ap.tar)){
    mapt <- try(mean.ap.tar(initial), silent = TRUE)
    if (inherits(mapt, "try-error"))
      stop("mean.ap.tar(initial) produced an error")
    if(!is.numeric(mapt))
      stop("mean.ap.tar must return a numeric vector")
    if(any(is.infinite(mapt))||any(is.na(mapt)))
      stop("mean.ap.tar(initial) is not finite")
    
    k=length(mapt)/dd
    if(k != round(k) || k < 1)
      stop(sprintf("length of mean.ap.tar(initial) must be a positive multiple of length(initial) = %d", dd))
    k <- as.integer(k)
    
    if(!is.null(var.ap.tar)){
      vat <- try(var.ap.tar(initial), silent = TRUE)
      if (inherits(vat, "try-error"))
        stop("var.ap.tar(initial) produced an error")
      
      if (dd == 1) {
        if (!is.numeric(vat) || length(vat) != k || any(vat <= 0)) {
          stop(sprintf(
            "var.ap.tar(x) must return a numeric vector of length %d with positive values when initial has length 1.",
            k
          ))
        }
      } else {
        if (!is.matrix(vat) || any(dim(vat) != c(dd, k * dd))) {
          stop(sprintf("var.ap.tar(x) must return a %d x %d covariance matrix.", dd, k * dd))
        }
        non_pd <- sapply(seq_len(k), function(ii) {
          block <- vat[, ((ii - 1) * dd + 1):(ii * dd)]
          !isSymmetric(block) || !is_pd(block)
        })
        if (any(non_pd)) {
          stop("var.ap.tar must be a set of positive definite symmetric matrices, each of dimension length(initial) x length(initial), combined by columns.")
        }
      }
    }
    if(!is.null(dens.ap.tar)){
      dap_len <- length(dens.ap.tar(initial, initial))
      if(k != dap_len)
        stop(sprintf("Inconsistency: dens.ap.tar returns %d values but mean.ap.tar implies k=%d components. These must match.", dap_len, k))
    }
  } else if(!is.null(dens.ap.tar)){
    dap_len <- length(dens.ap.tar(initial, initial))
    k <- dap_len
    if(k < 1)
      stop("dens.ap.tar must return at least one density value")
  } else{
    k=1
  }
  
  
  if(!is.null(samp.ap.tar)){
    for(ii in 1:k){
      sat <- try(samp.ap.tar(initial, kk=ii), silent = TRUE)
      if (inherits(sat, "try-error"))
        stop(sprintf("samp.ap.tar(initial, kk=%d) produced an error", ii))
      if(!is.numeric(sat))
        stop("samp.ap.tar must return a numeric vector")
      if(any(is.infinite(sat))||any(is.na(sat)))
        stop(sprintf("samp.ap.tar(initial, kk=%d) does not return finite values", ii))
      if(length(sat)!=dd) 
        stop(sprintf("samp.ap.tar must return a vector of length %d (same as initial)", dd))
    }
  }

  if(!is.null(bhat.coef)){
    bhval <- try(bhat.coef(initial), silent = TRUE)
    if (inherits(bhval, "try-error"))
      stop("bhat.coef(initial) produced an error")
    if(!is.numeric(bhval))
      stop("bhat.coef must return a numeric vector")
    if(any(is.infinite(bhval))||any(is.na(bhval)))
      stop("bhat.coef(initial) is not finite")
    if(any(bhval < 0))
      stop("bhat.coef must return non-negative values")
    if(any(bhval >1))
      stop("bhat.coef must return values less than or equal to one")
    bhval_len <- length(bhat.coef(initial))
    if(k != bhval_len)
      stop(sprintf("Inconsistency: bhat.coef returns %d values but mean.ap.tar implies k=%d components. These must match.", bhval_len, k))
    if(bhval_len < 1)
      stop("bhat.coef must return at least one value")
  }
  
  
  if(imp$enabled){
    if(imp$samp.base){
      if(is.null(samp.base)) 
        stop("samp.base must be provided to implement importance sampling with samp.base=TRUE")
    }else{
      if(is.null(samp.ap.tar)) 
        stop("samp.ap.tar must be provided to implement importance sampling with samp.base=FALSE")
    }
  }

  if (((is.null(mean.base) || is.null(var.base)) &&
       !is.null(dens.base) && !is.null(samp.base)) || ((is.null(mean.ap.tar) || is.null(var.ap.tar)) &&
                                                       !is.null(dens.ap.tar) && !is.null(samp.ap.tar))) {
    gaus=FALSE
  }

  diag.v.ap=FALSE

  base_incomplete <- any(is.null(dens.base), is.null(samp.base)) &&
    any(is.null(mean.base), is.null(var.base))

  ap_tar_incomplete <- any(is.null(dens.ap.tar), is.null(samp.ap.tar)) &&
    any(is.null(mean.ap.tar), is.null(var.ap.tar))

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

  if (!is.null(mean.base) && !is.null(var.base) &&
      (is.null(dens.base) || is.null(samp.base))) {
    dens.base <- function(y, x) {
      exp(ldens_mvnorm(y, mean.base(x), var.base(x)))
    }
    samp.base <- function(x) {
      rmvnorm(mean.base(x), var.base(x))
    }
  }

  if (!is.null(mean.ap.tar) && !is.null(var.ap.tar) &&
      (is.null(dens.ap.tar) || is.null(samp.ap.tar))) {
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
        return(rmvnorm(mu, Sig))
      } else {
        idx <- ((kk - 1) * dd + 1):(kk * dd)
        return(rmvnorm(mu[idx], Sig[, idx, drop = FALSE]))
      }
    }
  }
  
  if(gaus){
    bhat.coef <- NULL
  }

  if (isTRUE(show.progress)) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      warning("Install the `progress` package to show progress bars.")
      show.progress <- FALSE
    } else {
      pb <- progress::progress_bar$new(
        total = n.iter,
        format = "  [:bar] :percent (:current/:total)",
        clear = FALSE,
        width = 60
      )
    }
  }


  if(length(a)!=k){a=rep.int(1/k,k)}
  ctr_accep <-0
  samps <- matrix(NA,nrow=n.iter,ncol=dd)
  log.tar.store<-numeric(n.iter)
  curr <- initial
  log.tar_curr<-log.target(curr)
  eps <- eps
  
  calc  <- compute_prod_theta(curr, k, gaus, imp,
                              mean.base, var.base, dens.base, samp.base,
                              mean.ap.tar, var.ap.tar, dens.ap.tar, samp.ap.tar,
                              bhat.coef, diag.v.ap, dd)

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
                                     bhat.coef, diag.v.ap, dd)
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
    if (logr >= 0 || log(runif(1L)) < logr) {
      curr        <- proposed
      log.tar_curr <- log.tar_prop
      ctr_accep    <- ctr_accep + 1

      if (!ind) {
        prod_c  <- prod_p
        theta_c <- theta_p
      }
    }
    samps[i,] <- curr
    log.tar.store[i]<-log.tar_curr

    if (isTRUE(show.progress)) pb$tick()
  }

    if (model.case == "both_missing") {
    return(list(
      var.base = v.base,
      mean.ap.tar = m.ap.tar,
      var.ap.tar = v.ap.tar,
      samples = samps,
      log.p = log.tar.store,
      acceptance.rate = ctr_accep / n.iter,
      model.case = "No user-specified values detected; using default settings for base and ap.tar.",gaus=gaus,ind=ind
    ))
  } else if (model.case == "only_ap_tar_missing") {
    return(list(
      mean.ap.tar = m.ap.tar,
      var.ap.tar = v.ap.tar,
      samples = samps,
      log.p = log.tar.store,
      acceptance.rate = ctr_accep / n.iter,
      model.case = "No user-specified value detected; using the default setting for ap.tar.",gaus=gaus,ind=ind
    ))
  }else if (model.case == "only_base_missing") {
    return(list(
      var.base = v.base,
      samples = samps,
      log.p = log.tar.store,
      acceptance.rate = ctr_accep / n.iter,
      model.case = "No user-specified value detected; using the default setting for base.",gaus=gaus,ind=ind
    ))
  } else {
    return(list(
      samples = samps,
      log.p = log.tar.store,
      acceptance.rate = ctr_accep / n.iter,
      model.case = "User-specified settings for both base and ap.tar are being used.",gaus=gaus,ind=ind
    ))
  }
}
