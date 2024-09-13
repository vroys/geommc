#' The log-unnormalized posterior probability of a model for Bayesian
#' variable selection.
#' @rdname logp.vs
#' @description  Calculates the log-unnormalized posterior probability of a model.
#' @param model The indices of active variables.
#' @param X The \eqn{n\times p} covariate matrix without intercept.
#' @param y The response vector of length \eqn{n}.
#' @param lam0 The precision parameter for \eqn{\beta_0}. Default: 0 (corresponding to improper uniform prior).
#' @param a0 The shape parameter for prior on \eqn{\sigma^2}. Default: 0.
#' @param b0 The scale parameter for prior on \eqn{\sigma^2}. Default: 0.
#' @param lam The slab precision parameter.
#' @param w The prior inclusion probability of each variable.
#' @return The log-unnormalized posterior probability of the model.
#' @author Vivekananda Roy
#' @references Roy, V.(2024) A geometric approach to informative MCMC sampling https://arxiv.org/abs/2406.09010
#' @examples
#' n=50; p=100; nonzero = 3
#' trueidx <- 1:3
#' nonzero.value <- 4
#' TrueBeta <- numeric(p)
#' TrueBeta[trueidx] <- nonzero.value
#' rho <- 0.5
#' xone <- matrix(rnorm(n*p), n, p)
#' X <- sqrt(1-rho)*xone + sqrt(rho)*rnorm(n)
#' y <- 0.5 + X %*% TrueBeta + rnorm(n)
#' result <- geomc.vs(X=X, y=y)
#' logp.vs(result$median.model,X,y,lam = nrow(X)/ncol(X)^2,w = sqrt(nrow(X))/ncol(X))
#' @export

logp.vs <- function(model,X,y,lam0=0, a0=0, b0=0,lam,w)
{
  if(lam0<0) stop("lam0 must be a non-negative number")
  if(a0<0) stop("a0 must be a non-negative number")
  if(b0<0) stop("b0 must be a non-negative number")
  if(lam<=0) stop("lam must be a positive number")
  if(w<=0 || w>=1) stop("w must be a proper fraction")
  if(is.logical(model)) model = which(model);

  model = as.integer(model)
  if(any(any(model>ncol(X)),any(model<1))) stop("The model must be a subset of column numbers of X")
  
  nn = nrow(X)
  if(nn!=length(y)) stop("The number of rows of X must match with the length of y")
  logw = log(w/(1-w))
  if(lam0==0){
    mult.c=0.5*(nn-1)+a0
    add.c=2*b0
    if(b0==0){
      ys = scale(y)
    }else{
      ys = y-mean(y)
    }
    }else{
      ys = y-mean(y)
      mult.c=0.5*nn+a0
      add.c=2*b0+(nn*lam0*mean(y)^2)/(nn+lam0)
      }
  p0 = length(model)
  if(p0 == 0)
      return(-mult.c*log(add.c+sum(ys^2)))
x.g = scale(X[,model,drop=FALSE])
  xtx = crossprod(x.g) + diag(x = lam,nrow = p0)

  R = chol(xtx)

  z = backsolve(R,crossprod(x.g,ys),transpose = T)

  logp = 0.5*p0*log(lam) - sum(log(diag(R))) - mult.c*log(sum(ys^2) - sum(z^2)+add.c) + p0*logw
  
  return(logp)
}
