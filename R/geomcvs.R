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
#' @param model.threshold The threshold probability to select the covariates for
#' the median model (median.model) and the weighted average model (wam).
#' A covariate will be included in median.model (wam) if its marginal inclusion
#' probability (weighted marginal inclusion probability) is greater than the threshold. Default: 0.5.
#' @param lam0 The precision parameter for \eqn{\beta_0}. Default: 0 (corresponding to improper uniform prior).
#' @param a0 The shape parameter for prior on \eqn{\sigma^2}. Default: 0.
#' @param b0 The scale parameter for prior on \eqn{\sigma^2}. Default: 0.
#' @param lam The slab precision parameter. Default: \eqn{n/p^2}
#' as suggested by the theoretical results of Li, Dutta, Roy (2023).
#' @param w The prior inclusion probability of each variable. Default: \eqn{\sqrt{n}/p}.
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
#' geomc.vs provides MCMC samples from the posterior pmf of the models \eqn{P(\gamma|y)}, which is available
#' up to a normalizing constant.

#' geomc.vs also returns the marginal inclusion probabilities (mip)
#' computed by the Monte Carlo average as well as the weighted marginal inclusion
#' probabilities (wmip) computed with weights \deqn{w_i =
#' P(\gamma^{(i)}|y)/\sum_{k=1}^K P(\gamma^{(k)}|y), i=1,2,...,K} where \eqn{K} is the number of distinct
#' models sampled. Thus, based on the samples \eqn{\gamma^{(k)}, k=1,2,...,n.iter} mip for the \eqn{j}th variable is \deqn{mip_j =
#' \sum_{k=1}^{n.iter} I(\gamma^{(k)}_j = 1)/n.iter} and
#' wmip is as \deqn{wmip_j =
#' \sum_{k=1}^K w_k I(\gamma^{(k)}_j = 1).}
#' The median.model is the model containing variables \eqn{j} with \eqn{mip_j >
#' \code{model.threshold}} and the wam is the model containing variables \eqn{j} with \eqn{wmip_j >
#' \code{model.threshold}}. The conditional posterior mean of \eqn{\beta} given a model is
#' available in closed form. geomc.vs returns the posterior means of \eqn{\beta} conditional on the median.model and
#' the wam.
#' @return A list with components
#' \item{samples}{MCMC samples from  \eqn{P(\gamma|y)} returned as \eqn{p \times}n.iter sparse \code{dgCMatrix}.}
#' \item{\code{acceptance.rate}}{The acceptance rate based on all samples.}
#' \item{\code{mip}}{The \eqn{p} vector of marginal inclusion probabilities of all variables based on post burnin samples.}
#' \item{\code{median.model}}{The  median probability model based on post burnin samples.}
#' \item{\code{beta.med}}{The posterior  mean of \eqn{\beta} (the \eqn{p+1} vector c(intercept, regression
#'   coefficients)) conditional on the median.model based on post burnin samples.}
#'   \item{\code{wmip}}{The \eqn{p} vector of weighted marginal inclusion probabilities of all variables based on post burnin samples.}
#' \item{\code{wam}}{The weighted average model based on post burnin samples.}
#' \item{\code{beta.wam}}{The posterior  mean of \eqn{\beta} (the \eqn{p+1} vector c(intercept, regression
#'   coefficients)) conditional on the wam based on post burnin samples.}
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
#' result <- geomc.vs(X=X, y=y)
#' result$samples # the MCMC samples
#' result$acceptance.rate #the acceptance.rate
#' result$mip #marginal inclusion probabilities
#' result$wmip #weighted marginal inclusion probabilities
#' result$median.model #the median.model
#' result$wam #the weighted average model
#' result$beta.med #the posterior mean of regression coefficients for the median.model
#' result$beta.wam #the posterior mean of regression coefficients for the wam
#' result$log.post #the log (unnormalized) posterior probabilities of the MCMC samples.
#' @export
geomc.vs=function(X,y,initial=NULL,n.iter=50,burnin=1,eps=0.5,symm=TRUE, move.prob=c(.4,.4,.2), model.threshold = 0.5, lam0=0, a0=0, b0=0, lam = nrow(X)/ncol(X)^2, w = sqrt(nrow(X))/ncol(X)){
   if(n.iter < 1) stop("n.iter must be larger than or equal to 1")
   if(burnin < 1) stop("burnin must be larger than or equal to 1")
   if(burnin > n.iter) stop("burnin must be  smaller than or equal to n.iter")
   if(eps<=0 || eps>=1) stop("eps must be a proper fraction")
   if(lam0<0) stop("lam0 must be a non-negative number")
   if(a0<0) stop("a0 must be a non-negative number")
   if(b0<0) stop("b0 must be a non-negative number")
   if(lam<=0) stop("lam must be a positive number")
   if(w<=0 || w>=1) stop("w must be a proper fraction")
   if(model.threshold<=0 || model.threshold>=1) stop("model.threshold must be a proper fraction")
   if(!symm) if(any(move.prob <0) || length(move.prob)!=3) stop("move.prob must be a probability vector of length three")
   move.prob=move.prob/sum(move.prob)
   ctr_accep =0
   ncovar <- ncol(X)
   nn<-length(y)
   if(nn!=nrow(X)) stop("The number of rows of X must match with the length of y")
   if(any(any(initial>ncovar),any(initial<1))) stop("The initial model must be a subset of column numbers of X")

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
curr = sort.int(initial) # Initial value
eps = eps
log.post<-numeric(n.iter)
size <- integer(n.iter)
indices <- integer(n.iter*20)
ctr.ind=0

logp.add.cur=addvar_vs(curr, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, w=w, D=D, xbar=xbar)$logp
logp.del.cur=delvar_vs(curr, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, w=w, D=D, xbar=xbar)$logp
logp.swap.cur=swapvar_vs(curr, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, w=w, D=D, xbar=xbar)$logp
calc=thet_vs(curr,logp.add.cur,logp.del.cur,logp.swap.cur,ncovar,symm,move.prob)
prod_c=calc[,1]
theta_c=calc[,2]
logp_curr=logp.vs(curr,X,ys,lam0,a0,b0,lam,w)
        for (i in 1:n.iter) {
          proposed = samp_phi_vs(curr,prod_c,theta_c,ncovar,logp.add.cur,logp.del.cur,logp.swap.cur,symm,move.prob,X,ys,lam0,a0,b0,lam,w,eps)
          logp.add.p=addvar_vs(proposed, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, w=w, D=D, xbar=xbar)$logp
          logp.del.p=delvar_vs(proposed, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, w=w, D=D, xbar=xbar)$logp
          logp.swap.p=swapvar_vs(proposed, x=X, yty=yty, xty=Xty, mult.c, add.c, lam=lam, w=w, D=D, xbar=xbar)$logp
          calc=thet_vs(proposed,logp.add.p,logp.del.p,logp.swap.p,ncovar,symm,move.prob)
          prod_p=calc[,1]
          theta_p=calc[,2]
          logp_prop=logp.vs(proposed,X,ys,lam0,a0,b0,lam,w)
          logr = logp_prop+log_phi_vs(curr,proposed,prod_p,theta_p,ncovar,logp.add.p,logp.del.p,logp.swap.p,symm,move.prob,X,ys,lam0,a0,b0,lam,w,eps)-logp_curr-log_phi_vs(proposed,curr,prod_c,theta_c,ncovar,logp.add.cur,logp.del.cur,logp.swap.cur,symm,move.prob,X,ys,lam0,a0,b0,lam,w,eps)
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
        }

   indices <- indices[indices>0]
  cumsize <- cumsum(size)
   samps <- sparseMatrix(i=indices,p = c(0,cumsize),index1 = T,dims = c(ncovar,n.iter), x = 1)
   samps.postburn<-samps[,burnin:n.iter,drop=FALSE]
   log.postburn<-log.post[burnin:n.iter]
   logpost.uniq.ind <-  !duplicated(log.postburn)
   weight.unnorm = exp(log.postburn[logpost.uniq.ind] - max(log.postburn))
   weight = weight.unnorm/sum(weight.unnorm)

   MIP=rowMeans(samps.postburn)
   WMIP=as.vector(samps.postburn[,logpost.uniq.ind,drop=FALSE]%*%weight)
   med.model<-which(MIP>=model.threshold)
   model.WAM<-which(WMIP>=model.threshold)
   if (length(med.model)==0){
     beta.med<-numeric((ncovar+1))
     beta.med[1] <- (nn/(nn+lam0))*mean(y)
   }else{
     x.est <- cbind(rep(1, nn), scale(X[, med.model], center = xbar[med.model], scale = 1/D[med.model]))
     beta <- solve(crossprod(x.est) + diag(c(lam0, lam*rep(1, length(med.model)))), crossprod(x.est, y))
     beta.med<-numeric(ncovar)
     beta.med[med.model]<-beta[-1] * D[med.model]
     beta.med <- c(beta[1]-sum(beta[-1]*xbar[med.model]*D[med.model]),beta.med)
   }
   if (length(model.WAM)==0){
     beta.wam<-numeric((ncovar+1))
     beta.wam[1] <- (nn/(nn+lam0))*mean(y)
   }else{
     x.est1 <- cbind(rep(1, nn), scale(X[, model.WAM], center = xbar[model.WAM], scale = 1/D[model.WAM]))
     beta1 <- solve(crossprod(x.est1) + diag(c(lam0, lam*rep(1, length(model.WAM)))), crossprod(x.est1, y))
     beta.wam<-numeric(ncovar)
     beta.wam[model.WAM]<-beta1[-1] * D[model.WAM]
     beta.wam <- c(beta1[1]-sum(beta1[-1]*xbar[model.WAM]*D[model.WAM]),beta.wam)
     }
        return(list(samples=samps,acceptance.rate=ctr_accep/n.iter,mip=MIP,median.model=med.model,beta.med=beta.med,wmip=WMIP,wam=model.WAM,beta.wam=beta.wam,log.post=log.post))
        }
