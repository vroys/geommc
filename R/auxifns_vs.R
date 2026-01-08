#' @noRd
#' @keywords internal
thet_vs <- function(model, logp.add, logp.del, logp.swap, logp.best, log_cc, ncovar, symm, move.prob) {
  dim <- length(model)
  half_logp.add <- 0.5 * (logp.add - logp.best - log_cc)
  half_logp.del <- 0.5 * (logp.del - logp.best - log_cc)
  half_logp.swap <- 0.5 * (logp.swap - logp.best - log_cc)
  
  if (symm) {
    inv_sqrt_2ncovar <- 1 / sqrt(2 * ncovar)
    
    if (dim == 0) {
      prod <- sum(exp(half_logp.add)) * inv_sqrt_2ncovar
    } else if (dim == ncovar) {
      prod <- sum(exp(half_logp.del)) * inv_sqrt_2ncovar
    } else {
      prod <- (sum(exp(half_logp.add)) + sum(exp(half_logp.del))) * inv_sqrt_2ncovar +
        sum(exp(half_logp.swap)) / sqrt(2 * dim * (ncovar - dim))
    }
  } else {
    if (dim == 0) {
      prod <- sum(exp(half_logp.add)) * sqrt(move.prob[1] / ncovar)
    } else if (dim == ncovar) {
      prod <- sum(exp(half_logp.del)) * sqrt(move.prob[2] / ncovar)
    } else {
      prod <- sum(exp(half_logp.add)) * sqrt(move.prob[1] / (ncovar - dim)) +
        sum(exp(half_logp.del)) * sqrt(move.prob[2] / dim) +
        sum(exp(half_logp.swap)) * sqrt(move.prob[3] / (dim * (ncovar - dim)))
    }
  }
  
  prod_clamped <- pmax(pmin(prod, 1 - 1e-16), .Machine$double.eps)
  
  return(unname(cbind(prod, acos(prod_clamped))))
}
#model is the indices of \gamma==1
#' @noRd
samp_f_vs <- function(model, ncovar, symm, move.prob) {
  dimgam <- length(model)
  p <- ncovar
  
  if (dimgam > p)
    stop("model size must be smaller than the no. of covariates")
  
  # Boundary cases
  if (dimgam == 0L) {
    return(sample.int(p, 1L))
  }
  
  if (dimgam == p) {
    return(model[-sample.int(p, 1L)])
  }
  
  # Interior case
  if (symm) {
    move.prob <- c(0.5 * (p - dimgam) / p,
                   0.5 * dimgam / p,
                   0.5)
  }
  
  move <- sample.int(3L, 1L, prob = move.prob)
  
  if (move == 1L) {            # add
    out <- sample.int(p, 1L)
    while (out %in% model)
      out <- sample.int(p, 1L)
    model <- sort.int(c(model, out))
    
  } else if (move == 2L) {     # delete
    model <- model[-sample.int(dimgam, 1L)]
    
  } else {                     # swap
    in_pos <- sample.int(dimgam, 1L)
    out <- sample.int(p, 1L)
    while (out %in% model)
      out <- sample.int(p, 1L)
    model[in_pos] <- out
    model <- sort.int(model)
  }
  
  return(model=model)
}
dens_f_vs <- function(proposed, curr, ncovar, symm, move.prob) {
  p0 <- length(curr)
  p1 <- length(proposed)
  
  # Fast rejection
  if (abs(p1 - p0) > 1L)
    return(0)
  
  # Determine move type by size difference
  if (p1 == p0) {
    # could be swap or identical
    # check symmetric difference size == 2
    if (p0 == 0L || p0 == ncovar)
      return(0)
    
    # count mismatches
    n_diff <- sum(!curr %in% proposed) + sum(!proposed %in% curr)
    if (n_diff != 2L)
      return(0)
    
    if (symm) {
      move.prob <- c(0.5 * (ncovar - p0) / ncovar,
                     0.5 * p0 / ncovar,
                     0.5)
    }
    return(move.prob[3L] / (p0 * (ncovar - p0)))
  }
  
  if (p1 == p0 + 1L) {         # add
    if (symm) {
      move.prob <- c(0.5 * (ncovar - p0) / ncovar,
                     0.5 * p0 / ncovar,
                     0.5)
    }
    return(move.prob[1L] / (ncovar - p0))
  }
  
  if (p1 == p0 - 1L) {         # delete
    if (symm) {
      move.prob <- c(0.5 * (ncovar - p0) / ncovar,
                     0.5 * p0 / ncovar,
                     0.5)
    }
    return(move.prob[2L] / p0)
  }
  
  return(0)
}

dens_g_vs<- function(proposed,logp.best,log_cc,X,yty,Xty,mult.c,add.c,lam,logw){
  return(exp(logp_vs_in(proposed,X,yty,Xty,mult.c,add.c,lam,logw)-logp.best-log_cc))
}

samp_g_vs <- function(model, ncovar, logp.add, logp.del, logp.swap, logp.best)
{
  p0 <- length(model)
  p  <- ncovar
  
  if (p0 == 0L) {
    lp  <- logp.add - max(logp.add)
    idx <- sample.int(p, 1L, prob = exp(lp))
    return(idx)
  }
  
  if (p0 == p) {
    lp  <- logp.del - max(logp.del)
    pos <- sample.int(p0, 1L, prob = exp(lp))
    return(model[-pos])
  }
  
  
  # Weights of add, del, swap
  w_add  <- sum(exp(logp.add-logp.best))
  w_del  <- sum(exp(logp.del-logp.best))
  w_swap <- sum(exp(logp.swap-logp.best))
  
  move <- sample.int(3L, 1L, prob = c(w_add, w_del, w_swap))
  
  # ---------- ADD ----------
  if (move == 1L) {
    lp  <- logp.add - max(logp.add)
    idx <- sample.int(p, 1L, prob = exp(lp))
    return(sort.int(c(model, idx)))
  }
  
  # ---------- DELETE ----------
  if (move == 2L) {
    lp  <- logp.del - max(logp.del)
    idx <- sample.int(p0, 1L, prob = exp(lp[model]))
    return(model[-idx])
  }
  
  # ---------- SWAP ----------
  lp <- logp.swap - max(logp.swap)
  idx    <- sample.int(p * p0, 1L, prob = exp(lp))
  
  swapin      <- ((idx - 1L) %% p) + 1L
  swapout_pos <- ((idx - 1L) %/% p) + 1L
  
  model[swapout_pos] <- swapin
  
  return(sort.int(model))
}

#h(x|curr) with prod_curr
dens_h_vs=function(prop,curr,prod,ncovar,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw){
  if (prod ==1) return(numeric(1))
  
  dens_g_vals <- dens_g_vs(prop,logp.best,log.cc,X,yty,Xty,mult.c,add.c,lam,logw)
  dens_f_val  <- sqrt(dens_f_vs(prop,curr,ncovar,symm,move.prob))
  
  tmp <- sqrt(dens_g_vals) - prod * dens_f_val
  return((tmp^2) / (1 - prod^2))
}
#sample from u(.|curr)
samp_u_vs = function(curr,prod,ncovar,logp.add,logp.del,logp.swap,logp.best,symm,move.prob){
  threshold <- 1 / (1 + prod * prod)
    if(runif(1L) <= threshold) return(samp_g_vs(curr,ncovar,logp.add,logp.del,logp.swap,logp.best))
    else return(samp_f_vs(curr,ncovar,symm,move.prob))
}
#sample from h(.|curr)#is called only if prod<1
samp_h_vs = function(curr,prod,ncovar,logp.add,logp.del,logp.swap,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw){
  p2     <- prod*prod
  one_m_p2 <- 1 - p2
  log_one_m_p2 <- log(one_m_p2)
  
  max.try <- max(200000, 2000 * (1 + p2) / one_m_p2)
    success <- FALSE
    attempts <- 0
    while (!success) {
      attempts <- attempts + 1
      if (attempts > max.try) {
        stop(print(
          "Rejection sampler for h is taking too long to generate a valid proposal. 
         Try reducing eps or adjusting other tuning parameters (lam, w, move.prob), 
         or use a different starting model."))
      }
      
            prop=samp_u_vs(curr,prod,ncovar,logp.add,logp.del,logp.swap,logp.best,symm,move.prob)
            dg <- dens_g_vs(prop,logp.best,log.cc,X,yty,Xty,mult.c,add.c,lam,logw)
            df <- dens_f_vs(prop,curr,ncovar,symm,move.prob)
            num <- dens_h_vs(prop,curr,prod,ncovar,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw)
            log.accept <-log(num)+log_one_m_p2-log(dg+ p2*df)
            if (log(runif(1L)) < log.accept) {
              success <- TRUE
            }
  }
  return(prop)
}


log_phi_vs = function(prop,curr,prod,theta,ncovar,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw,eps){
  if(eps==1){
    p2<- prod^2 
    return(log(p2*dens_f_vs(prop,curr,ncovar,symm,move.prob)+(1-p2)*dens_h_vs(prop,curr,prod,ncovar,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw)))
    } else{
      coef=(cos(eps*theta))^2
      return(log(coef*dens_f_vs(prop,curr,ncovar,symm,move.prob)+(1-coef)*dens_h_vs(prop,curr,prod,ncovar,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw)))
    }
}
#sample from phi(|curr)
samp_phi_vs = function(model,prod,theta,ncovar,logp.add,logp.del,logp.swap,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw,eps){
  p2 <- prod^2
    if(p2==1){
        return(samp_f_vs(model,ncovar,symm,move.prob))
            }else{
              u = runif(1L)
    if(eps==1){
        if(u<=p2){
            return(samp_f_vs(model,ncovar,symm,move.prob))
        }else
        {return(samp_h_vs(model,prod,ncovar,logp.add,logp.del,logp.swap,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw))
            }
    }else{
        coef=cos(eps*theta)^2
        if(u<=coef){
          return(samp_f_vs(model,ncovar,symm,move.prob))
        }else{
          return(samp_h_vs(model,prod,ncovar,logp.add,logp.del,logp.swap,logp.best,log.cc,symm,move.prob,X,yty,Xty,mult.c,add.c,lam,logw))
          }
    }
    }
}
### A function to compute the log (unormalized) posterior probabilities for all models
### in the "added" set for the current model.
### codes of addvar_vs, delvar_vs and swapvar_vs are modified from the codes in the R package bravo
addvar_vs <- function (model,n,p, x, yty, xty, mult.c, add.c, lam, logw, R0 = NULL, v0 = NULL, D,xbar){
  p0 = length(model)

  xtx <- n - 1

  if (p0 == 0) { # There is no variable in the model
    R1 <- sqrt(xtx+lam)
    RSS <- yty - (xty/R1)^2+add.c
    logp <- 0.5*log(lam)-0.5*log(xtx+lam)-mult.c*log(RSS) + logw
    return(list(logp=logp))
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
  RSS <- {yty - sum(v0^2)} - u^2+add.c
  RSS[model] = 1 # whatever, they are going to be set to -Inf
  logp <- 0.5*(p0+1)*log(lam) - logdetR1 - mult.c*log(RSS) + (p0+1)*logw
  logp = as.numeric(logp)
  logp[model] <- -Inf
  return(list(logp=logp))
}
### A function to compute the log (unormalized) posterior probabilities for all models
### in the "deleted" set for the current model.
delvar_vs <- function(model,p, x, yty, xty, mult.c, add.c, lam, logw, D, xbar) {
  p0 <- length(model)
  logp <- numeric(p)
  if (p0 == 0) {  ##added it
    logp <- rep(-Inf,p) ##added it
  }else{
    if (p0 == 1) {
    logp.del <- -mult.c*log(add.c+yty)
    #RSS.del <- NULL ##added it
  } else {
    logp.del <- numeric(p0)
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
      logp.del[j] <- 0.5*(p0-1)*log(lam) - logdetR0 - mult.c*log(RSS0+add.c) + (p0-1)*logw
    }
  }
  logp[model] <- logp.del
  logp[-model] <- -Inf
  }
  return(list(logp=logp))
}


### A function to compute the log (unormalized) posterior probabilities for all models
### in the "swapped" set for the current model.
swapvar_vs <- function(model, n, p, x, yty, xty, mult.c, add.c, lam, logw, D, xbar) {
  p0 <- length(model)
  xtx <- n - 1
    if (p0 == 0) {  ##added it
    logp <- rep(-Inf,p) ##added it
  }else{
  #if (swapOnly == T) {
    logp <- matrix(0, nrow=p, ncol=p0)
    if (p0 == 1) {
      r1 <- addvar_vs(model=NULL, n, p, x=x, yty=yty, xty=xty, mult.c, add.c, lam=lam, logw=logw, D=D, xbar=xbar)
      logp[, 1] <- r1$logp
      logp[model, 1] <- -Inf
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
        RSS <- {yty - sum(v0^2)} - u^2+add.c
        RSS[model] = 1 # whatever, they are going to be set to -Inf
        logp1 <- 0.5*(p0)*log(lam) - logdetR1 - mult.c*log(RSS) + p0*logw
        logp[, j] = as.numeric(logp1)
        logp[model, j] <- -Inf
      }
    }
  }
  return(list(logp=logp))
}
