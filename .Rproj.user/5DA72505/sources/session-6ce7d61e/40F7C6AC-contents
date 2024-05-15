#' @noRd
#' @keywords internal
thet_vs=function(model,logp.add,logp.del,logp.swap,ncovar,symm,move.prob){
  dim=length(model)
  logp.best<-max(logp.add,logp.del,logp.swap)
  cc=sum(exp(logp.add-logp.best))+sum(exp(logp.del-logp.best))+sum(exp(logp.swap-logp.best))
  if(symm){
    if(dim ==0){
      logp.add1 <- logp.add-logp.best
      prod=sum(sqrt(exp(logp.add1)/cc))/sqrt(2*ncovar)
    }else if(dim ==ncovar){
      logp.del1 <- logp.del-logp.best
      prod=sum(sqrt(exp(logp.del1)/cc))/sqrt(2*ncovar)
    }else{
      logp.add1 <- logp.add-logp.best
      logp.del1 <- logp.del-logp.best
      logp.swap1 <- logp.swap-logp.best
      prod=sum(sqrt(c(exp(logp.add1),exp(logp.del1))/cc))/sqrt(2*ncovar)+sum(sqrt(exp(logp.swap1)/cc))/sqrt(2*dim*(ncovar-dim))
    }
  }else{
    if(dim ==0){
      logp.add1 <- logp.add-logp.best
      prod=sum(sqrt(exp(logp.add1)/cc))*sqrt(move.prob[1]/ncovar)
    }else if(dim ==ncovar){
      logp.del1 <- logp.del-logp.best
      prod=sum(sqrt(exp(logp.del1)/cc))*sqrt(move.prob[2]/ncovar)
    }else{
      logp.add1 <- logp.add-logp.best
      logp.del1 <- logp.del-logp.best
      logp.swap1 <- logp.swap-logp.best
      prod=sum(sqrt(exp(logp.add1)/cc))*sqrt(move.prob[1]/(ncovar-dim))+sum(sqrt(exp(logp.del1)/cc))*sqrt(move.prob[2]/dim)+sum(sqrt(exp(logp.swap1)/cc))*sqrt(move.prob[3]/(dim*(ncovar-dim)))
    }
  }
  return(unname(cbind(prod,acos(pmin(pmax(prod,0.0),1.0)))))
}
#model is the indices of \gamma==1
#' @noRd
samp_f_vs=function(model,ncovar,symm,move.prob){
  dimgam <- length(model)
  p<-ncovar
  if(dimgam>p){
    stop("model size must be smaller than the no. of covariates")
  }
    if (dimgam == 0|| dimgam==p){
      index <- sample.int(p,size=1)
      if (dimgam == 0){
        model<-index
      }else{
        model<-setdiff(seq(1,p),index)
      }
    }else{
      if(symm){
        move.prob=c(0.5*c(p-dimgam,dimgam)/p,0.5) #add, del, swap prob vector
        }
      move=sample.int(3,1,prob=move.prob)
      if(move==1){
        index <- sample(setdiff(seq(1,p),model),size=1)
        model<-sort.int(union(model,index))
    }else if(move==2){
      index <- sample(model,size=1)
      model<-setdiff(model,index)
    } else{
      swapin <- sample(model,size=1) # position of randomly chosen included variable
      swapout <- sample(setdiff(seq(1,p),model),size=1)# position of randomly chosen excluded variable
model=sort.int(setdiff(union(model,swapout),swapin))
    }
    }
  return(model=model)
}

dens_f_vs<- function(proposed,curr,ncovar,symm,move.prob){
  p1<-length(proposed)
  p0<-length(curr)
  if(abs(p1-p0)>1){ #we can drop this as the proposed values are only from nbd(curr)
    return(0)
  }else{
    if(symm){
              move.prob=c(0.5*c(ncovar-p0,p0)/ncovar,0.5) #add, del, swap prob vector
    }
      if(length(setdiff(union(proposed,curr),intersect(proposed,curr)))==2){
        return(move.prob[3]/(p0*(ncovar-p0)))
      }else if(length(setdiff(proposed,curr))==1){
          return(move.prob[1]/(ncovar-p0))
      }else if(length(setdiff(curr,proposed))==1){
           return(move.prob[2]/p0)
      }else{
        return(0)
      }
      }
}

dens_g_vs<- function(proposed,logp.add,logp.del,logp.swap,X,ys,lam,w){
  logp.best<-max(logp.add,logp.del,logp.swap)
cc<-sum(exp(logp.add-logp.best))+sum(exp(logp.del-logp.best))+sum(exp(logp.swap-logp.best))
  return(exp(logp.vs(proposed,X,ys,lam,w)-logp.best)/cc)
}

samp_g_vs=function(model,ncovar,logp.add,logp.del,logp.swap){
  p0 <- length(model)
  p<-ncovar
  if(p0>p){
    stop("model size must be smaller than the no. of covariates")
  }
  if(p0==0){
    logp.best<-max(logp.add)
    logp.add1<-logp.add-logp.best
    model <- sample.int(p,size=1,prob=exp(logp.add1))
  }else if(p0==p){
    logp.best<-max(logp.del)
    logp.del1<-logp.del-logp.best
       index<-sample.int(p,size=1,prob=exp(logp.del1))
        model<-setdiff(seq(1,p),index)
    }else{
      logp.best<-max(logp.add,logp.del,logp.swap)
      logp.add1 <- logp.add-logp.best
      logp.del1 <- logp.del-logp.best
      logp.swap1 <- logp.swap-logp.best
      move=sample.int(3,1,prob=c(sum(exp(logp.add1)),sum(exp(logp.del1)),sum(exp(logp.swap1))))
      if(move==1){
        #index <- sample(setdiff(seq(1,p),model),size=1,prob=exp(logp.add[logp.add>-Inf]))
        logp.best<-max(logp.add)
        logp.add1<-logp.add-logp.best
        index <- sample.int(p,size=1,prob=exp(logp.add1))
        model<-sort.int(union(model,index))
    }else if(move==2){
      #index <- sample(model,size=1,prob=exp(logp.del[logp.del>-Inf]))
      logp.best<-max(logp.del)
      logp.del1<-logp.del-logp.best
      index <- sample.int(p,size=1,prob=exp(logp.del1))
      model<-setdiff(model,index)
    } else{
      logp.best<-max(logp.swap)
      logp.swap1<-logp.swap-logp.best
      index<-sample.int(p0*p,1,prob=exp(logp.swap1))# position of randomly chosen element of swap matrix
      swapin_index <-ceiling(index/p)
      swapin <- model[swapin_index] # position of randomly chosen included variable
      swapout<- index-(swapin_index-1)*p # position of randomly chosen excluded variable
model=sort.int(setdiff(union(model,swapout),swapin))
    }
  }
  return(model=model)
}

#h(x|curr) with prod_curr #this returns a vector of length k
dens_h_vs=function(prop,curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w){
  kk=length(prod)
  val=numeric(kk)
  ind.j=which(prod<1)
  val[ind.j]=(sqrt(dens_g_vs(prop,logp.add,logp.del,logp.swap,X,ys,lam,w)[ind.j]) - prod[ind.j]*sqrt(dens_f_vs(prop,curr,ncovar,symm,move.prob)))^2/(1-prod[ind.j]^2)
    return(val)
}
#sample from u(.|curr)
samp_u_vs = function(curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,kk=1){
    if(runif(1)<=1/(1+prod[kk]^2)) return(samp_g_vs(curr,ncovar,logp.add,logp.del,logp.swap))
    else return(samp_f_vs(curr,ncovar,symm,move.prob))
}
#sample from h(.|curr)#is called only if prod[kk]<1
samp_h_vs = function(curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w,kk=1){
    success <- FALSE
    while (!success) {
          #M= (1+prod^2)/(1-prod^2)
            prop=samp_u_vs(curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,kk)
            success <- log(runif(1)) < log(dens_h_vs(prop,curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w)[kk])+log(1-prod[kk]^2)-log(dens_g_vs(prop,logp.add,logp.del,logp.swap,X,ys,lam,w)[kk]+ prod[kk]^2*dens_f_vs(prop,curr,ncovar,symm,move.prob))
  }
  return(prop)
}


log_phi_vs = function(prop,curr,prod,theta,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w,eps){
#  kk=length(a)
 # val=numeric(kk)
  kk=1
  val=0
  if(eps==1){
    val=log(prod^2*dens_f_vs(prop,curr,ncovar,symm,move.prob)+(1-prod^2)*dens_h_vs(prop,curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w))
    } else{
      coef=(cos(eps*theta))^2
      val=log(coef*dens_f_vs(prop,curr,ncovar,symm,move.prob)+(1-coef)*dens_h_vs(prop,curr,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w))
    }
  return(sum(val))
#  return(sum(a*val))
}
#sample from phi(|curr)
samp_phi_vs = function(model,prod,theta,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w,eps){
  #kk=length(a)
  kk=1
  s.a=1
  #if(kk>1){
  #  s.a=sample.int(kk,1,prob=a)
   # }
    if(prod[s.a]==1){
        return(samp_f_vs(model,ncovar,symm,move.prob))
            }else{
    if(eps==1){
        u = runif(1)
        if(u<=prod[s.a]^2){
            return(samp_f_vs(model,ncovar,symm,move.prob))
        }else
        {return(samp_h_vs(model,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w,s.a))
            }
    }else{
        coef=(cos(eps*theta[s.a]))^2
        u = runif(1)
        if(u<=coef){
          return(samp_f_vs(model,ncovar,symm,move.prob))
        }else{
          return(samp_h_vs(model,prod,ncovar,logp.add,logp.del,logp.swap,symm,move.prob,X,ys,lam,w,s.a))
          }
    }
    }
}
### A function to compute the log (unormalized) posterior probabilities for all models
### in the "added" set for the current model.
### codes of addvar_vs, delvar_vs and swapvar_vs are modified from the codes in the R package bravo
addvar_vs <- function (model, x, ys, xty, lam, w, R0 = NULL, v0 = NULL, D,xbar){
  n <- nrow(x)
  p <- ncol(x)
  p0 = length(model)


  yty <- n - 1
  xtx <- n - 1
  logw <- log(w/(1-w))

  if (p0 == 0) { # There is no variable in the model
    R1 <- sqrt(xtx+lam)
    RSS <- yty - (xty/R1)^2
    logp <- 0.5*log(lam)-0.5*log(xtx+lam)-0.5*(n-1)*log(RSS) + logw
    j = which.max(logp)
    return(list(logp=logp,R=sparseMatrix(i=1,j=1,x=R1,triangular=T),
                v=xty[j]/R1,which.max=j, RSS=RSS))
  }

  # else there's at least one variable

  D0 <- D[model]
  x0bar <- xbar[model]
  x0t <- D0*{t(x[,model,drop=FALSE]) - x0bar}

  if(is.null(R0)) {
    R0 = chol(tcrossprod(x0t) + diag(x = lam,nrow = p0));
  }
  if(is.null(v0)) {
    v0 = backsolve(R0,xty[model],transpose = T)
  }


  S1 <- backsolve(R0,x0t,transpose = T);
  if(p0 == 1) S1 = matrix(S1,nrow = 1)

  S1 <-  S1 - rowMeans(S1); # For numerical stability.

  S  <- S1 %*% x;

  S <- S %*% Diagonal(p,x=D)

  if(class(S)[1] == "dgeMatrix") {
    sts <- colSumSq_dge(S@x,S@Dim)
  } else {
    sts <- colSumSq_matrix(S)
  }

  sts[model] <- 0;
  s0 <- sqrt({xtx+lam} - sts)

  u <- (xty-crossprod(S, v0))/s0
  u[model] = 0;
  logdetR1 <- sum(log(diag(R0))) + log(s0)
  RSS <- {yty - sum(v0^2)} - u^2
  RSS[model] = 1 # whatever, they are going to be set to -Inf
  logp <- 0.5*(p0+1)*log(lam) - logdetR1 - 0.5*(n-1)*log(RSS) + (p0+1)*logw
  logp = as.numeric(logp)
  logp[model] <- -Inf
  j = which.max(logp)

  # now update R1tinv and v1 and return as a list
  sj = S[,j]
  s0j = s0[j]
  R1 = rbind(cbind(R0,sj),c(rep(0,p0),s0j))
  R1 =as.matrix(R1)                       #methods::as(R1,"dtCMatrix") depricated
  v1 = c(v0,u[j])

  RSS[model] <- -Inf
  return(list(logp=logp, R = R1, v=v1, which.max=j, RSS=RSS))
}
### A function to compute the log (unormalized) posterior probabilities for all models
### in the "deleted" set for the current model.
delvar_vs <- function(model, x, xty, lam, w, D, xbar) {
  n <- nrow(x)
  p <- ncol(x)
  p0 <- length(model)
  yty <- n-1
  logw <- log(w/(1-w))
  logp <- numeric(p)
  if (p0 == 0) {  ##added it
    logp <- rep(-Inf,p) ##added it
    RSS.del <- NULL ##added it
  }else{
    if (p0 == 1) {
    logp.del <- -(n-1)/2*log(yty)
    RSS.del <- NULL ##added it
  } else {
    logp.del <- numeric(p0)
    RSS.del <- numeric(p0)
    x0 <- scale(x[, model, drop=F])
    xgx <- crossprod(x0) + lam*diag(p0)
    for (j in 1:p0) {
      # delete one variable in the current model
      model.temp <- model[-j]
      R0 <- chol(xgx[-j, -j])
      logdetR0 <- sum(log(diag(R0)))
      if(is.nan(logdetR0)) logdetR0 = Inf
      RSS0 <- yty - sum(backsolve(R0, xty[model.temp], transpose = T)^2)
      if(RSS0 <= 0) RSS0 = .Machine$double.eps
      logp.del[j] <- 0.5*(p0-1)*log(lam) - logdetR0 - 0.5*(n-1)*log(RSS0) + (p0-1)*logw
      RSS.del[j] <- RSS0
    }
  }
  logp[model] <- logp.del
  logp[-model] <- -Inf
  RSS.del[model] <- RSS.del
    RSS.del[-model] <- -Inf
  }
  return(list(logp=logp, RSS=RSS.del))
}


### A function to compute the log (unormalized) posterior probabilities for all models
### in the "swapped" set for the current model.
swapvar_vs <- function(model, x, ys, xty, lam, w, D, xbar, swapOnly=F) {
  n <- nrow(x)
  p <- ncol(x)
  p0 <- length(model)
  xtx <- n - 1
  yty <- n-1
  logw <- log(w/(1-w))
    if (p0 == 0) {  ##added it
    logp <- rep(-Inf,p) ##added it
    RSS.mat <- NULL ##added it
  }else{
  if (swapOnly == T) {
    logp <- matrix(0, nrow=p, ncol=p0)
    RSS.mat <- matrix(0, nrow=p, ncol=p0)
    if (p0 == 1) {
      r1 <- addvar_vs(model=NULL, x=x, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
      logp[, 1] <- r1$logp
      logp[model, 1] <- -Inf
      RSS.mat[, 1] <- r1$RSS
      RSS.mat[model, 1] <- -Inf
    } else {
      x0 <- scale(x[, model, drop=F])
      xgx <- crossprod(x0) + lam*diag(p0)
      x0x1 <- crossprod(x0, x)
      x0x <- x0x1 %*% Diagonal(p,x=D)
      for (j in 1:p0) {
        # delete one variable in the current model
        model.temp <- model[-j]
        R0 <- chol(xgx[-j, -j])
        logdetR0 <- sum(log(diag(R0)))
        if(is.nan(logdetR0)) logdetR0 = Inf
        v0 <- backsolve(R0, xty[model.temp], transpose = T)

        # add back another variable from the remaning variables
        S <- backsolve(R0, x0x[-j, , drop=F], transpose = T)
        if(length(model.temp) == 1) S = matrix(S, nrow = 1)
        if(class(S)[1] == "dgeMatrix") {
          sts <- colSumSq_dge(S@x,S@Dim)
        } else {
          sts <- colSumSq_matrix(S)
        }
        sts[model] <- 0;
        s0 <- sqrt({xtx+lam} - sts)

        u <- (xty-crossprod(S, v0))/s0
        u[model] = 0;
        logdetR1 <- sum(log(diag(R0))) + log(s0)
        RSS <- {yty - sum(v0^2)} - u^2
        RSS[model] = 1 # whatever, they are going to be set to -Inf
        logp1 <- 0.5*(p0)*log(lam) - logdetR1 - 0.5*(n-1)*log(RSS) + p0*logw
        logp[, j] = as.numeric(logp1)
        logp[model, j] <- -Inf
        RSS.mat[, j] <- RSS
        RSS.mat[model, j] <- -Inf
      }
    }
  } else {
    logp <- matrix(0, nrow=p, ncol=p0+1)
    RSS.mat <- matrix(0, nrow=p, ncol=p0+1)
    if (p0 == 1) {
      logp.del <- -(n-1)/2*log(yty)
      RSS.del <- 0
      r1 <- addvar_vs(model=NULL, x=x, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
      logp[, 2] <- r1$logp
      logp[model, 2] <- -Inf
      RSS.mat[, 2] <- r1$RSS
      RSS.mat[model, 2] <- -Inf
    } else {
      logp.del <- numeric(p0)
      RSS.del <- numeric(p0)
      x0 <- scale(x[, model, drop=F])
      xgx <- crossprod(x0) + lam*diag(p0)
      x0x1 <- crossprod(x0, x)
      x0x <- x0x1 %*% Diagonal(p,x=D)
      for (j in 1:p0) {
        # delete one variable in the current model
        model.temp <- model[-j]
        R0 <- chol(xgx[-j, -j])
        logdetR0 <- sum(log(diag(R0)))
        if(is.nan(logdetR0)) logdetR0 = Inf
        v0 <- backsolve(R0, xty[model.temp], transpose = T)
        RSS0 <- yty - sum(backsolve(R0, xty[model.temp], transpose = T)^2)
        if(RSS0 <= 0) RSS0 = .Machine$double.eps
        logp.del[j] <- 0.5*(p0-1)*log(lam) - logdetR0 - 0.5*(n-1)*log(RSS0) + (p0-1)*logw
        RSS.del[j] <- RSS0

        # add back another variable from the remaning variables
        S <- backsolve(R0, x0x[-j, , drop=F], transpose = T)
        if(length(model.temp) == 1) S = matrix(S, nrow = 1)
        if(class(S)[1] == "dgeMatrix") {
          sts <- colSumSq_dge(S@x,S@Dim)
        } else {
          sts <- colSumSq_matrix(S)
        }
        sts[model] <- 0;
        s0 <- sqrt({xtx+lam} - sts)

        u <- (xty-crossprod(S, v0))/s0
        u[model] = 0;
        logdetR1 <- sum(log(diag(R0))) + log(s0)
        RSS <- {yty - sum(v0^2)} - u^2
        RSS[model] = 1 # whatever, they are going to be set to -Inf
        logp1 <- 0.5*(p0)*log(lam) - logdetR1 - 0.5*(n-1)*log(RSS) + p0*logw
        logp[, j+1] = as.numeric(logp1)
        logp[model, j+1] <- -Inf
        RSS.mat[, j+1] <- RSS
        RSS.mat[model, j+1] <- -Inf
      }
    }
    logp[model, 1] <- logp.del
    logp[-model, 1] <- -Inf
    RSS.mat[model, 1] <- RSS.del
    RSS.mat[-model, 1] <- -Inf
  }
  }
  return(list(logp=logp, RSS=RSS.mat))
}
