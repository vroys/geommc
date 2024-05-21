#' Markov chain Monte Carlo for discrete and continuous
#' distributions using geometric MH algorithms.
#' @rdname geomc
#' @importFrom stats rnorm runif
#' @importFrom matrixcalc is.positive.definite
#' @importFrom cubature divonne
#' @importFrom magrittr %>%
#' @description geomc produces Markov chain samples from a target distribution.
#' The target can be a pdf or pmf. Users specify the target distribution by an R function that evaluates
#' the log un-normalized pdf or pmf. geomc uses the geometric approach of Roy (2024)  to move an uninformed
#' base density (e.g. a random walk proposal) towards different global/local approximations of the
#' target density. The base density can be specified along with its mean, covariance matrix, and a function
#' for sampling from it. Gaussian densities can be specified by mean and variance only, although it is preferred to supply its density
#' and sampling functions as well. If either or both of the mean and variance arguments and any of the density and sampling functions is
#' missing, then a base density corresponding to a random walk with an appropriate scale parameter is used. One or more approximate target densities
#' can be specified along with their means, covariance matrices, and a function for sampling from the densities.
#' Gaussian densities can be specified by mean and variance only, although it is preferred to supply their densities and sampling
#' functions as well. If either or both of the mean and variance
#' arguments and any of the density and sampling functions is missing for the approximate target density, then a normal distribution with mean computed from
#' a pilot run of a random walk Markov chain and a diagonal covariance matrix with a large variance is used.
#' If the Argument gaus is set as FALSE then both the base and the approximate target can be specified by their
#' densities and functions for sampling from it. That is, if gaus=FALSE, the functions specifying the means and variances of
#' both the base and the approximate target densities are not used.
#' If the target is a pmf (discrete distribution), then gaus=FALSE and imp \eqn{[1]}=TRUE (not the default values) need to be specified.
#' @param log.target is the logarithm of the (un-normalized) target density function, needs to be written as a function of the current value \eqn{x}.
#' @param initial is the initial values.
#' @param n.iter is the no. of samples needed.
#' @param eps is the value for epsilon perturbation. Default is 0.5.
#' @param ind is False if either the base density, \eqn{f} or the approximate target density, \eqn{g} depends on \eqn{x}. Default is False.
#' @param gaus is True if both \eqn{f} and \eqn{g} are normal distributions. Default is True.
#' @param imp is a vector of three elements. If gaus is TRUE, then the imp argument is not used.
#' imp \eqn{[1]} is False  if numerical integration is used, otherwise, importance sampling is used. Default is False.
#' imp \eqn{[2]} (n.samp) is no of samples in importance sampling.
#' imp \eqn{[3]} (samp.base) is True if samples from \eqn{f} is used, otherwise samples from \eqn{g} is used. Default is False.
#' @param a is the probability vector for the mixture proposal density. Default is the uniform distribution.
#' @param mean.base is the mean of the base density \eqn{f}, needs to be written as a function of current value \eqn{x}.
#' @param var.base is the covariance matrix of the base density \eqn{f}, needs to be written as a function of current value \eqn{x}.
#' @param dens.base is the density function of the base density \eqn{f}, needs to be written as a function \eqn{(y,x)} (in this order) of the proposed value \eqn{y} and the current value \eqn{x}, although it may not depend on \eqn{x}.
#' @param samp.base is the function to draw from the base density \eqn{f}, needs to be written as a function of current value \eqn{x}.
#' @param mean.ap.tar is the vector of means of the densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function of current value \eqn{x}. It must have the same dimension as \eqn{k} times the length of initial.
#' @param var.ap.tar is the matrix of covariance matrices of the densities \eqn{g_i(y|x), i=1,\dots,k} formed by column concatenation. It needs to be written as a function of current value \eqn{x}. It must have the same dimension as the length of initial by \eqn{k} times the length of initial.
#' @param dens.ap.tar is the vector of densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function \eqn{(y,x)} (in this order) of the proposed value \eqn{y} and the current value \eqn{x}, although it may not depend on \eqn{x}.
#' @param samp.ap.tar is the function to draw from the densities \eqn{g_i(y|x), i=1,\dots,k}. It needs to be written as a function of (current value \eqn{x}, the indicator of mixing component \eqn{kk}). It must return a vector of the length of that of the initial.
#' @details
#' Using a geometric Metropolis-Hastings algorithm geom.mc produces Markov chains with the target as its stationary distribution. The details
#' of the method can be found in Roy (2024).
#' @return The function returns a list with the following elements:
#'  \item{\code{samples}}{A matrix containing the MCMC samples. Each column is one sample.}
#' \item{\code{acceptance.rate}}{The acceptance rate.}
#' @author Vivekananda Roy <vroy@iastate.edu>
#' @references Roy, V.(2024) A geometric approach to informative MCMC sampling
#' @examples
#' result <- geomc(log.target=function(y) dnorm(y,log=TRUE),initial=0,n.iter=500) #univariate normal
#' result$samples # the MCMC samples.
#' result$acceptance.rate # the acceptance rate.
#' result<-geomc(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' initial=0,n.iter=500) #mixture of univariate normals, default choices
#' hist(result$samples)
#' result<-geomc(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' initial=0,n.iter=500, mean.base = function(x) x,var.base= function(x) 4,
#' dens.base=function(y,x) dnorm(y, mean=x,sd=2),samp.base=function(x) x+2*rnorm(1),
#' mean.ap.tar=function(x) c(0,10),var.ap.tar=function(x) c(1,1.4^2),
#' dens.ap.tar=function(y,x) c(dnorm(y),dnorm(y,mean=10,sd=1.4)),
#' samp.ap.tar=function(x,kk=1){if(kk==1){return(rnorm(1))} else{return(10+1.4*rnorm(1))}})
#' #mixture of univariate normals, an informed choice for dens.ap.tar
#' hist(result$samples)
#' samp.ap.tar=function(x,kk=1){s.g=sample.int(2,1,prob=c(.5,.5))
#' if(s.g==1){return(rnorm(1))
#' }else{return(10+1.4*rnorm(1))}}
#' result<-geomc(log.target=function(y) log(0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4)),
#' initial=0,n.iter=500,gaus=FALSE,imp=c(TRUE,n.samp=100,samp.base=TRUE),
#' dens.base=function(y,x) dnorm(y, mean=x,sd=2),samp.base=function(x) x+2*rnorm(1),
#' dens.ap.tar=function(y,x) 0.5*dnorm(y)+0.5*dnorm(y,mean=10,sd=1.4),
#' samp.ap.tar=samp.ap.tar)
#' #mixture of univariate normals, another informed choice for dens.ap.tar
#' hist(result$samples)
#' result <- geomc(log.target=function(y) -0.5*crossprod(y),initial=rep(0,4),
#' n.iter=500) #multivariate normal
#' size=5
#' result <- geomc(log.target = function(y) dbinom(y, size, 0.3, log = TRUE),
#' initial=0,n.iter=500,ind=TRUE,gaus=FALSE,imp=c(TRUE,n.samp=1000,samp.base=TRUE),
#' dens.base=function(y,x) 1/(size+1), samp.base= function(x) sample(seq(0,size,1),1),
#' dens.ap.tar=function(y,x) dbinom(y, size, 0.7),samp.ap.tar=function(x,kk=1) rbinom(1, size, 0.7))
#'  #binomial
#'  table(result$samples)
#' @export

geomc=function(log.target,initial,n.iter,eps=0.5,ind=FALSE,gaus=TRUE,imp=c(FALSE,n.samp=1000,samp.base=FALSE),a=1,mean.base,var.base,dens.base,samp.base,mean.ap.tar,var.ap.tar,dens.ap.tar,samp.ap.tar){
  if(missing(log.target)) stop("log.target must be provided")

  if(!checkFuncArgs(log.target, "x")) stop("log.target must be a function with only one argument \"x\"")

  if(missing(n.iter)) stop("n.iter must be provided")

  if(n.iter < 1) stop("n.iter must be larger than or equal to 1")


  dd=length(initial)
  if(eps<0 || eps>1) stop("eps must be a proper fraction")

  if(!missing(mean.base)){
    if(!checkFuncArgs(mean.base, "x")) stop("mean.base must be a function with only one argument \"x\"")

    if(length(mean.base(initial))!=dd) stop("mean.base must return a vector of same length as initial")
  }
  if(!missing(var.base)){
    if(!checkFuncArgs(var.base, "x"))stop("var.base must be a function with only one argument \"x\"")

    if(length(var.base(initial))!=dd^2)stop("x and var.base have non-conforming size")

    if((dd==1 && var.base(initial)<0) ||(dd>1 &&!is.positive.definite(var.base(initial), tol = sqrt(.Machine$double.eps)))){
      stop("var.base must be a positive definite matrix of dimension length(initial) times length(initial)")
    }
  }
  if(!missing(dens.base)){
    if(!checkFuncArgs(dens.base, c("y", "x")))stop("dens.base must be a function with only two arguments \"y\" and \"x\"")
  }

  if(!missing(samp.base)){
    if(!checkFuncArgs(samp.base, "x")) stop("samp.base must be a function with only one argument \"x\"")
  }

  if(!missing(dens.ap.tar)){
    if(!checkFuncArgs(dens.ap.tar, c("y", "x")))stop("dens.ap.tar must be a function with only two arguments \"y\" and \"x\"")
  }


  if(!missing(mean.ap.tar)){
    if(!checkFuncArgs(mean.ap.tar, "x")) stop("mean.ap.tar must be a function with only one argument \"x\"")
    k=length(mean.ap.tar(initial))/dd
    if(!missing(var.ap.tar)){
      if(!checkFuncArgs(var.ap.tar, "x")) stop("var.ap.tar must be a function with only one argument \"x\"")

      if(length(var.ap.tar(initial))!=k*dd^2) stop("x and var.ap.tar have non-conforming size")

      if((dd==1 && any(var.ap.tar(initial)<0)) ||(dd>1 && any(sapply(seq(1,k),function(ii) !is.positive.definite(var.ap.tar(initial)[,((ii-1)*dd+1):(ii*dd)], tol = sqrt(.Machine$double.eps)))))){
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
    if(dd!=length(samp.ap.tar(initial,kk=1))) stop("samp.ap.tar must return a sample of the same size as initial")
  }

  if(imp[1]){
    if(imp[3]){
      if(missing(samp.base)) stop("samp.base must be provided to implement importance sampling with samp.base=TRUE")
    }else{
      if(missing(samp.ap.tar)) stop("samp.ap.tar must be provided to implement importance sampling with samp.base=FALSE")
    }
  }

  if(all(!missing(mean.base),!missing(var.base),any(missing(dens.base),missing(samp.base)))){
    dens.base=function(y,x) exp(ldens_mvnorm(y, mean.base(x),var.base(x)))
    samp.base=function(x) if(dd==1){return(mean.base(x)+sqrt(var.base(x))*rnorm(1))}else{return(mean.base(x)+chol(var.base(x))%*%rnorm(dd))}
  }
  if(all(!missing(mean.ap.tar),!missing(var.ap.tar),any(missing(dens.ap.tar),missing(samp.ap.tar)))){
    dens.ap.tar=function(y,x) {
      if(k==1){return(exp(ldens_mvnorm(y,mean.ap.tar(x),var.ap.tar(x))))}else{
        mean.ap=mean.ap.tar(x)
        var.ap= var.ap.tar(x)
        if(dd==1){return(sapply(seq(1,k),function(ii) exp(ldens_mvnorm(y,mean.ap[ii],var.ap[ii]))))}else{
          return(sapply(seq(1,k),function(ii) exp(ldens_mvnorm(y,mean.ap[((ii-1)*dd+1):(ii*dd)],var.ap[,((ii-1)*dd+1):(ii*dd)]))))
        }
      }
    }
    samp.ap.tar=function(x,kk=1) {
      if(k==1){if(dd==1){return(mean.ap.tar(x)+sqrt(var.ap.tar(x))*rnorm(1))
      }else{return(mean.ap.tar(x)+chol(var.ap.tar(x))%*%rnorm(dd))}}else{
        mean.ap=mean.ap.tar(x)
        var.ap= var.ap.tar(x)
        if(dd==1){return(mean.ap[kk]+sqrt(var.ap[kk])*rnorm(1))
        }else{return(mean.ap[((kk-1)*dd+1):(kk*dd)]+chol(var.ap[,((kk-1)*dd+1):(kk*dd)])%*%rnorm(dd))}
      }
    }
  }
  TargetOnly=FALSE
  if(any(all(any(missing(dens.base),missing(samp.base)),any(missing(mean.base),missing(var.base))),all(any(missing(dens.ap.tar),missing(samp.ap.tar)),any(missing(mean.ap.tar),missing(var.ap.tar))))){

    fgfns=genfvargmean(log.target,initial,dd)
    if(all(all(any(missing(dens.base),missing(samp.base)),any(missing(mean.base),missing(var.base))),all(any(missing(dens.ap.tar),missing(samp.ap.tar)),any(missing(mean.ap.tar),missing(var.ap.tar))))){
      TargetOnly=TRUE
      dens.base=function(y,x) exp(ldens_mvnorm_diag(y, x,fgfns$var.base))
      samp.base=function(x) x+fgfns$sd.base*rnorm(dd)
      mean.base=function(x) x
      var.base=function(x) fgfns$var.base*diag(dd)
      k=1
      dens.ap.tar=function(y,x) exp(ldens_mvnorm_diag(y,fgfns$mean.ap.tar, 30^2))
      samp.ap.tar=function(x,kk=1) x+30*rnorm(dd)
      mean.ap.tar=function(x) fgfns$mean.ap.tar
      var.ap.tar=function(x)  30^2*diag(dd)
    }else if(all(any(all(!missing(dens.base),!missing(samp.base)),all(!missing(mean.base),!missing(var.base))),all(any(missing(dens.ap.tar),missing(samp.ap.tar)),any(missing(mean.ap.tar),missing(var.ap.tar))))){
      k=1
      dens.ap.tar=function(y,x) exp(ldens_mvnorm_diag(y,fgfns$mean.ap.tar, 30^2))
      samp.ap.tar=function(x,kk=1) x+30*rnorm(dd)
      mean.ap.tar=function(x) fgfns$mean.ap.tar
      var.ap.tar=function(x)  30^2*diag(dd)
    }else{
      dens.base=function(y,x) exp(ldens_mvnorm_diag(y, x,fgfns$var.base))
      samp.base=function(x) x+fgfns$sd.base*rnorm(dd)
      mean.base=function(x) x
      var.base=function(x) fgfns$var.base*diag(dd)
    }
  }

  if(length(a)!=k){a=rep.int(1/k,k)}
  ctr_accep =0
  samps = matrix(NA,nrow=dd,ncol=n.iter)
  curr = initial
  log.tar_curr=log.target(curr)
  eps = eps
  a=a/sum(a)
  if(!gaus){
    calc=thet(k,curr,imp,dens.base,dens.ap.tar,samp.base,samp.ap.tar)
    if(!ind){
      prod_c=calc[,1]
      theta_c=calc[,2]
      for (i in 1:n.iter) {
        proposed = as.vector(samp_phi(a,curr,prod_c,theta_c,eps,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
        calc=thet(k,proposed,imp,dens.base,dens.ap.tar,samp.base,samp.ap.tar)
        prod_p=calc[,1]
        theta_p=calc[,2]
        log.tar_prop=log.target(proposed)
        logr = (log.tar_prop+log_phi(a,curr,proposed,prod_p,theta_p,eps,dens.base,dens.ap.tar)-log.tar_curr-log_phi(a,proposed,curr,prod_c,theta_c,eps,dens.base,dens.ap.tar))%>%replace(is.na(.),-Inf)
        if (logr >=0 || log(runif(1)) < logr) {
          curr = proposed
          prod_c=prod_p
          theta_c=theta_p
          log.tar_curr=log.tar_prop
          ctr_accep = ctr_accep+1
        }
        samps[,i] = curr
      }
      return(list(samples=samps,acceptance.rate=ctr_accep/n.iter))
    }else{
      prod=calc[,1]
      theta=calc[,2]
      for (i in 1:n.iter) {
        proposed = as.vector(samp_phi(a,curr,prod,theta,eps,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
        log.tar_prop=log.target(proposed)
        logr = (log.tar_prop+log_phi(a,curr,proposed,prod,theta,eps,dens.base,dens.ap.tar)-log.tar_curr-log_phi(a,proposed,curr,prod,theta,eps,dens.base,dens.ap.tar))%>%replace(is.na(.),-Inf)
        if (logr >=0 || log(runif(1)) < logr) {
          curr = proposed
          log.tar_curr=log.tar_prop
          ctr_accep = ctr_accep+1
        }
        samps[,i] = curr
      }
      return(list(samples=samps,acceptance.rate=ctr_accep/n.iter))
    }
  }else{
    if(!ind){
        if(dd==1){
          calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(curr),mean.ap.tar(curr)[seq((ii-1)*dd+1,ii*dd)],var.base(curr),matrix(var.ap.tar(curr),nrow=1)[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
      }else{
      calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(curr),mean.ap.tar(curr)[seq((ii-1)*dd+1,ii*dd)],var.base(curr),var.ap.tar(curr)[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
      }
      #calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(curr),mean.ap.tar(curr)[seq((ii-1)*dd+1,ii*dd)],var.base(curr),as.matrix(var.ap.tar(curr))[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
      prod_c=calc[,1]
      theta_c=calc[,2]
      for (i in 1:n.iter) {
        proposed = as.vector(samp_phi(a,curr,prod_c,theta_c,eps,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
        log.tar_prop=log.target(proposed)
        if(dd==1){
          calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(proposed),mean.ap.tar(proposed)[seq((ii-1)*dd+1,ii*dd)],var.base(proposed),matrix(var.ap.tar(proposed),nrow=1)[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
        }else{
          calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(proposed),mean.ap.tar(proposed)[seq((ii-1)*dd+1,ii*dd)],var.base(proposed),var.ap.tar(proposed)[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
        }
        #calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(proposed),mean.ap.tar(proposed)[seq((ii-1)*dd+1,ii*dd)],var.base(proposed),as.matrix(var.ap.tar(proposed))[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
        prod_p=calc[,1]
        theta_p=calc[,2]
        logr = (log.tar_prop+log_phi(a,curr,proposed,prod_p,theta_p,eps,dens.base,dens.ap.tar)-log.tar_curr-log_phi(a,proposed,curr,prod_c,theta_c,eps,dens.base,dens.ap.tar))%>%replace(is.na(.),-Inf)
        if (logr >=0 || log(runif(1)) < logr) {
          curr = proposed
          prod_c=prod_p
          theta_c=theta_p
          log.tar_curr=log.tar_prop
          ctr_accep = ctr_accep+1
        }
        samps[,i] = curr
      }
      return(list(samples=samps,acceptance.rate=ctr_accep/n.iter))
    }else{
      if(dd==1){
        calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(curr),mean.ap.tar(curr)[seq((ii-1)*dd+1,ii*dd)],var.base(curr),matrix(var.ap.tar(curr),nrow=1)[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
      }else{
        calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(curr),mean.ap.tar(curr)[seq((ii-1)*dd+1,ii*dd)],var.base(curr),var.ap.tar(curr)[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
      }
      #calc=matrix(unlist(lapply(seq(1,k),function(ii) bc(mean.base(curr),mean.ap.tar(curr)[seq((ii-1)*dd+1,ii*dd)],var.base(curr),as.matrix(var.ap.tar(curr))[,seq((ii-1)*dd+1,ii*dd),drop=FALSE],TargetOnly))), ncol = 2, byrow = TRUE)
      prod=calc[,1]
      theta=calc[,2]
      for (i in 1:n.iter) {
        proposed = as.vector(samp_phi(a,curr,prod,theta,eps,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
        log.tar_prop=log.target(proposed)
        logr = (log.tar_prop+log_phi(a,curr,proposed,prod,theta,eps,dens.base,dens.ap.tar)-log.tar_curr-log_phi(a,proposed,curr,prod,theta,eps,dens.base,dens.ap.tar))%>%replace(is.na(.),-Inf)
        if (logr >=0 || log(runif(1)) < logr) {
          curr = proposed
          log.tar_curr=log.tar_prop
          ctr_accep = ctr_accep+1
        }
        samps[,i] = curr
      }
      return(list(samples=samps,acceptance.rate=ctr_accep/n.iter))
    }
  }
}
