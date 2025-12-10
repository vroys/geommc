#' Markov chain Monte Carlo for Bayesian variable selection
#' using a geometric MH algorithm.
#' @rdname geomc.vs
#' @importFrom stats sd
#' @importFrom Matrix Diagonal sparseMatrix tcrossprod crossprod diag colMeans colSums t rowMeans rowSums solve
#' @description geomc.vs uses a geometric approach to MCMC for performing Bayesian variable selection.
#' It produces MCMC samples from the posterior density of a Bayesian hierarchical feature selection model.
#' @param X The \eqn{n\times p} covariate matrix without intercept. The following classes are supported:
#' \code{matrix} and \code{dgCMatrix}.  No need to center or scale this matrix manually. Scaling is performed implicitly and
#' regression coefficients are returned on the original scale.
#' @param y The response vector of length \eqn{n}. No need to center or scale.
#' @param initial is the initial model (the set of active variables). Default: Null model.
#' @param n.iter is the no. of samples needed. Default: 50.
#' @param burnin is the value of burnin used to compute the median probability model. Default: 1.
#' @param eps is the value for epsilon perturbation. Default: 0.5.
#' @param symm indicates if the base density is of symmetric RW-MH. Default: True.
#' @param move.prob is the vector of ('addition', 'deletion', 'swap') move probabilities. Default: (0.4,0.4,0.2).
#' move.prob is used only when symm is set to False.
#' @param lam0 The precision parameter for \eqn{\beta_0}. Default: 0 (corresponding to improper uniform prior).
#' @param a0 The shape parameter for prior on \eqn{\sigma^2}. Default: 0.
#' @param b0 The scale parameter for prior on \eqn{\sigma^2}. Default: 0.
#' @param lam The slab precision parameter. Default: \eqn{n/p^2}
#' as suggested by the theoretical results of Li, Dutta, Roy (2023).
#' @param w The prior inclusion probability of each variable. Default: \eqn{\sqrt{n}/p}.
#' @param model.summary If true, additional summaries are returned. Default: FALSE.
#' @param model.threshold The threshold probability to select the covariates for
#' the median model (median.model) and the weighted average model (wam).
#' A covariate will be included in median.model (wam) if its marginal inclusion
#' probability (weighted marginal inclusion probability) is greater than the threshold. Default: 0.5.
#' @param show.progress Logical. Whether to show progress during sampling. Default: TRUE.
#' @details
#' geomc.vs provides MCMC samples using the geometric MH algorithm of Roy (2024)
#' for variable selection based on a hierarchical Gaussian linear model with priors placed
#' on the regression coefficients as well as on the model space as follows:
#' \deqn{y | X, \beta_0,\beta,\gamma,\sigma^2,w,\lambda \sim N(\beta_01 + X_\gamma\beta_\gamma,\sigma^2I_n)}
#' \deqn{\beta_i|\beta_0,\gamma,\sigma^2,w,\lambda \stackrel{indep.}{\sim} N(0, \gamma_i\sigma^2/\lambda),~i=1,\ldots,p,}
#' \deqn{\beta_0|\gamma,\sigma^2,w,\lambda \sim N(0, \sigma^2/\lambda_0)}
#' \deqn{\sigma^2|\gamma,w,\lambda \sim Inv-Gamma (a_0, b_0)}
#' \deqn{\gamma_i|w,\lambda \stackrel{iid}{\sim} Bernoulli(w)}
#' where \eqn{X_\gamma} is the \eqn{n \times |\gamma|} submatrix of \eqn{X} consisting of
#' those columns of \eqn{X} for which \eqn{\gamma_i=1} and similarly, \eqn{\beta_\gamma} is the
#' \eqn{|\gamma|} subvector of \eqn{\beta} corresponding to \eqn{\gamma}. The density \eqn{\pi(\sigma^2)} of 
#' \eqn{\sigma^2 \sim Inv-Gamma (a_0, b_0)} has the form \eqn{\pi(\sigma^2) \propto (\sigma^2)^{-a_0-1} \exp(-b_0/\sigma^2)}.
#'  The functions in the package also allow the non-informative prior \eqn{(\beta_{0}, \sigma^2)|\gamma, w \sim 1 / \sigma^{2}} 
#'  which is obtained by setting \eqn{\lambda_0=a_0=b_0=0}.
#' geomc.vs provides the empirical MH acceptance rate and MCMC samples from the posterior pmf of the models \eqn{P(\gamma|y)}, which is available
#' up to a normalizing constant.

#' If \eqn{\code{model.summary}} is set TRUE, geomc.vs also returns other model summaries. In particular, it returns the 
#' marginal inclusion probabilities (mip) computed by the Monte Carlo average as well as the weighted marginal 
#' inclusion probabilities (wmip) computed with weights \deqn{w_i =
#' P(\gamma^{(i)}|y)/\sum_{k=1}^K P(\gamma^{(k)}|y), i=1,2,...,K} where \eqn{\gamma^{(k)}, k=1,2,...,K} are the distinct
#' models sampled. Thus, if \eqn{N_k} is the no. of times the \eqn{k}th distinct model \eqn{\gamma^{(k)}} is repeated in the MCMC samples,
#' the mip for the \eqn{j}th variable is \deqn{mip_j =
#' \sum_{k=1}^{K} N_k I(\gamma^{(k)}_j = 1)/n.iter} and
#' wmip for the \eqn{j}th variable is \deqn{wmip_j =
#' \sum_{k=1}^K w_k I(\gamma^{(k)}_j = 1).}
#' The median.model is the model containing variables \eqn{j} with \eqn{mip_j >
#' \code{model.threshold}} and the wam is the model containing variables \eqn{j} with \eqn{wmip_j >
#' \code{model.threshold}}. Note that \eqn{E(\beta|\gamma, y)}, the conditional posterior mean of \eqn{\beta} given a model \eqn{\gamma} is
#' available in closed form (see Li, Dutta, Roy (2023) for details). geomc.vs returns two estimates (beta.mean, beta.wam) of the posterior mean 
#' of \eqn{\beta} computed as
#' \deqn{ beta.mean = \sum_{k=1}^{K} N_k E(\beta|\gamma^{(k)},y)/n.iter} and
#' \deqn{beta.wam = \sum_{k=1}^K w_k E(\beta|\gamma^{(k)},y),} respectively.
#' @return A list with components
#' \item{samples}{MCMC samples from  \eqn{P(\gamma|y)} returned as a n.iter\eqn{\times p} sparse \code{lgCMatrix}.}
#' \item{\code{acceptance.rate}}{The acceptance rate based on all samples.}
#' \item{\code{mip}}{The \eqn{p} vector of marginal inclusion probabilities of all variables based on post burnin samples.}
#' \item{\code{median.model}}{The  median probability model based on post burnin samples.}
#' \item{\code{beta.mean}}{The Monte Carlo estimate of posterior  mean of \eqn{\beta} (the \eqn{p+1} vector c(intercept, regression
#'   coefficients)) based on post burnin samples.}
#'   \item{\code{wmip}}{The \eqn{p} vector of weighted marginal inclusion probabilities of all variables based on post burnin samples.}
#' \item{\code{wam}}{The weighted average model based on post burnin samples.}
#' \item{\code{beta.wam}}{The model probability weighted estimate of posterior mean of \eqn{\beta} (the \eqn{p+1} vector c(intercept, regression
#'   coefficients)) based on post burnin samples.}
#' \item{\code{log.post}}{The n.iter vector of log of the unnormalized marginal posterior pmf \eqn{P(\gamma|y)} evaluated
#'   at the samples.}
#' @author Vivekananda Roy
#' @references Roy, V.(2024) A geometric approach to informative MCMC
#'   sampling https://arxiv.org/abs/2406.09010
#' @references Li, D., Dutta, S., Roy, V.(2023) Model Based Screening Embedded Bayesian 
#' Variable Selection for Ultra-high Dimensional Settings, Journal of Computational and 
#' Graphical Statistics, 32, 61-73
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
#' result <- geomc.vs(X=X, y=y,model.summary = TRUE)
#' result$samples # the MCMC samples
#' result$acceptance.rate #the acceptance.rate
#' result$mip #marginal inclusion probabilities
#' result$wmip #weighted marginal inclusion probabilities
#' result$median.model #the median.model
#' result$wam #the weighted average model
#' result$beta.mean #the posterior mean of regression coefficients
#' result$beta.wam #another estimate of the posterior mean of regression coefficients
#' result$log.post #the log (unnormalized) posterior probabilities of the MCMC samples.
#' @export
geomc.vs=function(X,y,initial=NULL,n.iter=50,burnin=1,eps=0.5,symm=TRUE, move.prob=c(.4,.4,.2), lam0=0, a0=0, b0=0, lam = nrow(X)/ncol(X)^2, w = sqrt(nrow(X))/ncol(X), model.summary=FALSE, model.threshold = 0.5,show.progress=TRUE){
  check_positive_integer(n.iter, "n.iter")
  check_positive_integer(burnin, "burnin")
  if(burnin > n.iter) stop("burnin must be  smaller than or equal to n.iter")
  check_fraction_01(eps, "eps")
  
  if (!is.numeric(lam0) || length(lam0) != 1 || is.na(lam0) || lam0 < 0)
    stop("lam0 must be a single non-negative number")
  
  if (!is.numeric(a0) || length(a0) != 1 || is.na(a0) || a0 < 0)
    stop("a0 must be a single non-negative number")
  
  if (!is.numeric(b0) || length(b0) != 1 || is.na(b0) || b0 < 0)
    stop("b0 must be a single non-negative number")
  
  if (!is.numeric(lam) || length(lam) != 1 || is.na(lam) || lam <= 0)
    stop("lam must be a single positive number")
  
  if (!is.numeric(w) || length(w) != 1 || is.na(w) || w <= 0 || w >= 1)
    stop("w must be a single proper fraction in (0,1)")
  
  if (!is.numeric(model.threshold) || length(model.threshold) != 1 ||
      is.na(model.threshold) || model.threshold <= 0 || model.threshold >= 1)
    stop("model.threshold must be a single proper fraction in (0,1)")
  if (!symm) {
    if (!is.numeric(move.prob) || length(move.prob) != 3 ||
        any(is.na(move.prob)) || any(move.prob < 0) ||
        sum(move.prob) == 0)
      stop("move.prob must be a numeric probability vector of length 3 with non-negative entries")
    
    move.prob <- move.prob / sum(move.prob)
  }
   ctr_accep =0
   ncovar <- ncol(X)
   nn <- length(y)
   
   if (nn != nrow(X))
     stop("The number of rows of X must match the length of y")
   
   if (any(is.na(X)) || any(is.na(y)))
     stop("X and y must not contain missing values")

   if (is.null(initial)) {
     initial <- integer(0)
   } else {
     if (!is.numeric(initial) || any(initial %% 1 != 0))
       stop("initial must be NULL or a vector of integer column indices")
     if (any(initial < 1) || any(initial > ncovar))
       stop("initial must be a subset of column numbers of X")
     initial <- as.integer(initial)
   }
   
  xbar = colMeans(X)
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
  if(class(X)[1] == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X,xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }
  Xty = D*as.numeric(crossprod(X,ys))
  yty=sum(ys^2)
  logw = log(w/(1-w))
curr = sort.int(initial) 
eps = eps
log.post<-numeric(n.iter)
size <- integer(n.iter)
indices <- integer(n.iter*20)
ctr.ind=0

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

logp.add.cur=addvar_vs(curr,n=nn,p=ncovar, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, logw, D=D, xbar=xbar)$logp
logp.del.cur=delvar_vs(curr, p=ncovar,x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, logw, D=D, xbar=xbar)$logp
logp.swap.cur=swapvar_vs(curr, n=nn,p=ncovar, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, logw, D=D, xbar=xbar)$logp
calc=thet_vs(curr,logp.add.cur,logp.del.cur,logp.swap.cur,ncovar,symm,move.prob)
prod_c=calc[,1]
theta_c=calc[,2]
logp_curr=logp.vs.in(curr,X,yty,Xty,mult.c,add.c,lam,logw)
        for (i in 1:n.iter) {
          proposed = samp_phi_vs(curr,prod_c,theta_c,ncovar,logp.add.cur,logp.del.cur,logp.swap.cur,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw,eps)
          logp.add.p=addvar_vs(proposed, n=nn,p=ncovar,x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, logw, D=D, xbar=xbar)$logp
          logp.del.p=delvar_vs(proposed, p=ncovar,x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, logw, D=D, xbar=xbar)$logp
          logp.swap.p=swapvar_vs(proposed, n=nn,p=ncovar,x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, logw, D=D, xbar=xbar)$logp
          calc=thet_vs(proposed,logp.add.p,logp.del.p,logp.swap.p,ncovar,symm,move.prob)
          prod_p=calc[,1]
          theta_p=calc[,2]
          logp_prop=logp.vs.in(proposed,X,yty,Xty,mult.c,add.c,lam,logw)
          logr = logp_prop+log_phi_vs(curr,proposed,prod_p,theta_p,ncovar,logp.add.p,logp.del.p,logp.swap.p,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw,eps)-logp_curr-log_phi_vs(proposed,curr,prod_c,theta_c,ncovar,logp.add.cur,logp.del.cur,logp.swap.cur,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw,eps)
            if (logr >=0 || log(runif(1)) < logr){
                curr = proposed
                prod_c=prod_p
                theta_c=theta_p
                logp.add.cur=logp.add.p
                logp.del.cur=logp.del.p
                logp.swap.cur=logp.swap.p
                logp_curr=logp_prop
                ctr_accep = ctr_accep+1
            }
          log.post[i] <-logp_curr
           if(length(curr)>0){
          size[i] <- length(curr)
          indices[(ctr.ind+1):(ctr.ind+size[i])] <- curr
          ctr.ind <- ctr.ind+size[i]
           }
          if (isTRUE(show.progress)) pb$tick()
        }
indices <- indices[indices>0]
cumsize <- cumsum(size)
samps <- sparseMatrix(j=indices,p = c(0,cumsize),index1 = T,dims = c(n.iter,ncovar), x = T)
if(!model.summary){
  return(list(samples=samps,acceptance.rate=ctr_accep/n.iter))
}
samps.postburn<-samps[burnin:n.iter,,drop=FALSE]
log.postburn<-log.post[burnin:n.iter]
logpost.uniq.ind <-  !duplicated(log.postburn)
log.postburn.uniq<-log.postburn[logpost.uniq.ind]
weight.unnorm <- exp(log.postburn.uniq - max(log.postburn.uniq))
weight <- weight.unnorm/sum(weight.unnorm)
n.rep<-sapply(log.postburn.uniq, FUN=function(x){sum(x==log.postburn)})
samps.uniq<-samps.postburn[logpost.uniq.ind,,drop=FALSE]
modelsize.uniq<-rowSums(samps.uniq)
no.model.uniq<-sum(logpost.uniq.ind)
MIP<-colMeans(samps.postburn)
WMIP <- as.vector(weight %*% samps.uniq)
med.model<-which(MIP>=model.threshold)
model.WAM<-which(WMIP>=model.threshold)
beta.est <- matrix(0, (ncovar+1),no.model.uniq)

for(i in 1:no.model.uniq){
  if(modelsize.uniq[i]==0){
    beta.est[1, i] <- (nn/(nn+lam0))*mean(y)
  }else{
    model_i <- samps.uniq[i, ] == 1
    x.est <- cbind(rep(1, nn), scale(X[, model_i, drop = FALSE], center = xbar[model_i], scale = 1/D[model_i]))
    beta <- solve(crossprod(x.est) + diag(c(lam0, lam*rep(1, modelsize.uniq[i]))), crossprod(x.est, y))
    beta.est[c(T, model_i), i] <- c(beta[1]-sum(beta[-1]*xbar[model_i]*D[model_i]),beta[-1] * D[model_i])
  }
}
beta.m<-rowSums(beta.est%*%diag(n.rep))/(n.iter-burnin+1)
beta.wam <- rowSums(beta.est%*%diag(weight))
return(list(samples=samps,acceptance.rate=ctr_accep/n.iter,mip=MIP,median.model=med.model,beta.mean=beta.m,wmip=WMIP,wam=model.WAM,beta.wam=beta.wam,log.post=log.post))
}
