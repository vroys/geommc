#' The log-unnormalized posterior probability of a model for Bayesian
#' variable selection.
#' @rdname logp.vs
#' @description  Calculates the log-unnormalized posterior probability of a model.
#' @param model The indices of active variables.
#' @param X The \eqn{n\times p} covariate matrix without intercept.
#' @param y The response vector of length \eqn{n}.
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

logp.vs <- function(model,X,y,lam,w)
{
  if(lam<=0) stop("lam must be a positive number")
  if(w<=0 || w>=1) stop("w must be a proper fraction")
  if(is.logical(model)) model = which(model);

  model = as.integer(model)
  if(any(any(model>ncol(X)),any(model<1))) stop("The model must be a subset of column numbers of X")
  

  n = nrow(X)
  logw = log(w/(1-w))

  p0 = length(model)
  if(p0 == 0)
    return(-0.5*(n-1)*log(n-1))

 ys = scale(y)
x.g = scale(X[,model,drop=FALSE])
  xtx = crossprod(x.g) + diag(x = lam,nrow = p0)

  R = chol(xtx)

  z = backsolve(R,crossprod(x.g,ys),transpose = T)

  logp = 0.5*p0*log(lam) - sum(log(diag(R))) - 0.5*(n-1)*log(sum(ys^2) - sum(z^2)) + p0*logw

  return(logp)
}
