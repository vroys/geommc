#' @keywords internal
#' @noRd
check_positive_integer <- function(x, name) {
  if (!is.numeric(x) ||
      length(x) != 1 ||
      is.na(x) ||
      is.infinite(x) ||
      x <= 0 ||
      x != as.integer(x)) {
    stop(sprintf("%s must be a positive finite integer", name))
  }
}
check_fraction_01 <- function(x, name) {
  if (!is.numeric(x) ||
      length(x) != 1 ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0 ||
      x > 1) {
    stop(sprintf("%s must be a finite numeric scalar in [0, 1]", name))
  }
}
###log of multivariate normal density
ldens_mvnorm=function(y,mean,Sigma){
  y    <- as.numeric(y)
  mean <- as.numeric(mean)
  d    <- length(y)
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  z <- backsolve(R,y-mean,transpose=TRUE)
  logdet <- sum(log(diag(R)))
  return(-logdet-0.5*sum(z^2)-0.5*d*log(2 *pi))
}
###log of multivariate normal density for constant variance sig
ldens_mvnorm_diag=function(y,mean,sig){
  y    <- as.numeric(y)
  mean <- as.numeric(mean)
  d    <- length(y)
  if (!is.finite(sig) || sig <= 0) return(-Inf)
  sqdist <- sum((y - mean)^2)
  return(-0.5*sqdist/sig-0.5*d*log(2*pi*sig))
}
###log of multivariate normal density when cholesky factor of sig is given
ldens_mvnormchol=function(y,mean,chol.sig){
  y    <- as.numeric(y)
  mean <- as.numeric(mean)
  d    <- length(y)
  if (!is.matrix(chol.sig) || nrow(chol.sig) != ncol(chol.sig)) return(-Inf)
  z <- backsolve(chol.sig,y-mean,transpose=TRUE)
  logdet <-sum(log(diag(chol.sig)))
  return(-logdet-0.5*sum(z^2)-0.5*d*log(2 *pi))
}
rmvnorm_chol <- function(mu, Sigma) {
  d <- length(mu)
  z <- rnorm(d)
  if (d == 1) {
    return(mu + sqrt(Sigma) * z)
  }
  return(mu + crossprod(chol(Sigma), z))
}
is_pd <- function(mat, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(mat)) return(FALSE)
  ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  all(ev > tol)
}
####calculate msjd for MC output## mat is npara times niter matrix
msejd<-function(mat){
  if (!is.matrix(mat)) {
    mat <- matrix(mat, nrow = 1)
  }
  niter <- ncol(mat)
  if (niter < 2) return(0) 
  diffs <- mat[, 2:niter, drop = FALSE] - mat[, 1:(niter - 1), drop = FALSE]
  return(sum(diffs * diffs) / (niter - 1))
}
bc = function(mu1, mu2, sig1, sig2,diag.var){
  ddd=length(mu1)
  if(ddd>1){
    if(!diag.var){
      sig = 0.5*(sig1+sig2)
      R1 = chol(sig1)
      R2 = chol(sig2)
      R = chol(sig)
      z = backsolve(R,mu1-mu2,transpose=T)
      prod=exp(-sum(z^2)/8-(sum(log(diag(R)))-0.5*(sum(log(diag(R1)))+sum(log(diag(R2))))))
      }else{
        sig = 0.5*(sig1[1,1]+sig2[1,1])
        prod=exp(-sum((mu1-mu2)^2)/8/sig-0.5*ddd*(log(sig)-0.5*log(sig1[1,1]*sig2[1,1])))
        
      }
    }else{
    prod=exp(-(mu1-mu2)^2/(4*(sig1+sig2))-0.5*(log(sig1+sig2)-(log(2*sqrt(sig1*sig2)))))
  }
  return(unname(cbind(prod,acos(pmin(pmax(prod,.Machine$double.eps),1-1e-16)))))
}
#cmpute prod_curr= <f(y|curr), g(y|curr)> and \theta_curr
thet=function(curr,k,imp,dens.base,dens.ap.tar,samp.base,samp.ap.tar){
  prod=numeric(k)
  if (!imp$enabled) {
    prod <- sapply(seq_len(k), function(j) {
      out <- divonne(
        function(y) {
          f.ap  <- dens.ap.tar(y, curr)[j]
          f.bas <- dens.base(y, curr)
          sqrt(f.ap * f.bas)
        },
        lowerLimit = rep(-Inf, length(curr)),
        upperLimit = rep(Inf,  length(curr)),
        absTol = 1e-2
      )$integral
    })
  } else {
    if (imp$samp.base) {
      samples_list <- lapply(seq_len(imp$n.samp), function(i) {
        out <- samp.base(curr)
        matrix(out, ncol = 1)
      })
      samp.store <- do.call(cbind, samples_list)
      imp.store=apply(samp.store,2,function (y) sqrt(dens.ap.tar(y,curr)/dens.base(y,curr)))
      imp.store[!is.finite(imp.store)] <- 0
      prod <- if (k > 1) rowMeans(imp.store) else mean(imp.store)
    } else {
      for (j in seq_len(k)) {
        samples_list <- lapply(seq_len(imp$n.samp), function(i) {
          out <- samp.ap.tar(curr, kk = j)
          matrix(out, ncol = 1)
        })
        samp.store <- do.call(cbind, samples_list)
        imp.store=apply(samp.store,2,function (y) sqrt(dens.base(y,curr)/dens.ap.tar(y,curr)[j]))
        imp.store[!is.finite(imp.store)] <- 0
        prod[j] <- mean(imp.store)
      }
    }
  }
  angle <- acos(pmin(pmax(prod, .Machine$double.eps), 1 - 1e-16))
  return(unname(cbind(prod, angle)))
}
#cmpute prod_curr= <f(y|curr), g(y|curr)> and \theta_curr either by bc or thet
compute_prod_theta <- function(x, k, gaus, imp,
                               mean.base, var.base, dens.base, samp.base,
                               mean.ap.tar, var.ap.tar, dens.ap.tar, samp.ap.tar,
                               diag.v.ap, dd) {
  if (!gaus) {
    return(thet(x,k,imp,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
  }
  mb  <- mean.base(x)
  map <- mean.ap.tar(x)
  vb  <- var.base(x)
  vap <- var.ap.tar(x)
  
  out <- matrix(0, nrow = k, ncol = 2)
  
  for (ii in seq_len(k)) {
    idx <- ((ii-1)*dd + 1):(ii*dd)
    mean_ap_i <- map[idx]
    var_ap_i  <- if (dd == 1)
      matrix(vap, nrow = 1)[, idx, drop = FALSE]
    else
      vap[, idx, drop = FALSE]
    out[ii, ] <- bc(
      mb,
      mean_ap_i,
      vb,
      var_ap_i,
      diag.v.ap
    )
  }
  return(out)
}
#h(x|curr) with prod_curr #this returns a vector of length k
dens_h=function(x,curr,prod,dens.base,dens.ap.tar){
  kk=length(prod)
  val=numeric(kk)
  ind.j=which(prod<1)
  val[ind.j]=(sqrt(dens.ap.tar(x,curr)[ind.j]) - prod[ind.j]*sqrt(dens.base(x,curr)))^2/(1-prod[ind.j]^2)
  return(val)
}
#sample from u(.|curr)
samp_u = function(curr,prod,kk=1,samp.base,samp.ap.tar){
  if(runif(1)<=1/(1+prod[kk]^2)) return(samp.ap.tar(curr,kk))
  else return(samp.base(curr))
}
#sample from h(.|curr)#is called only if prod[kk]<1
samp_h = function(curr,prod,kk=1,dens.base,dens.ap.tar,samp.base,samp.ap.tar){
  success <- FALSE
  while (!success) {
    x<-samp_u(curr,prod,kk,samp.base,samp.ap.tar)
    rhs<- log(dens_h(x,curr,prod,dens.base,dens.ap.tar)[kk])+log1p(-prod[kk]^2)-log(dens.ap.tar(x,curr)[kk]+ prod[kk]^2*dens.base(x,curr))
    rhs[!is.finite(rhs)] <- -Inf
    success <- log(runif(1)) < rhs
  }
  return(x)
}

log_phi = function(a,x,curr,prod,theta,eps,dens.base,dens.ap.tar){
  kk=length(a)
  val=numeric(kk)
  if (eps == 1) {
    coef <- prod^2
  } else {
    coef <- cos(eps * theta)^2
  }
  coef <- min(max(coef, 0), 1)
  f1 <- dens.base(x, curr)                   
  f2 <- dens_h(x, curr, prod, dens.base, dens.ap.tar)
  log_f1 <- log(f1)
  log_f2 <- log(f2)
  log_c1 <- if (coef > 0) log(coef) else -Inf
  log_c2 <- if (coef < 1) log(1 - coef) else -Inf
  m <- pmax(log_c1 + log_f1, log_c2 + log_f2)
  val <- m + log(exp(log_c1 + log_f1 - m) + exp(log_c2 + log_f2 - m))
  return(sum(a*val))
}
#sample from phi(|curr)
samp_phi = function(a,curr,prod,theta,eps,dens.base,dens.ap.tar,samp.base,samp.ap.tar){
  kk=length(a)
  s.a <- if (kk > 1) sample.int(kk, 1, prob = a) else 1
  if(prod[s.a]==1){
    return(samp.base(curr))
  }
  if (eps == 1) {
    weight <- prod[s.a]^2
  } else {
    weight <- (cos(eps * theta[s.a]))^2
  }
  if (runif(1) <= weight) {
    return(samp.base(curr))
  } else {
    return(samp_h(curr, prod, s.a,
                  dens.base, dens.ap.tar,
                  samp.base, samp.ap.tar))
  }
}
genfvar <- function(log.target, initial, dd) {
  out <- rw_mc_cpp(log.target, initial, n_iter = 500, sig= 2.38^2 /dd)
  if (out$acceptance_rate >= 0.25 && out$acceptance_rate <= 0.7) {
    sig.rw <- 2.38^2 /dd
  } else if (out$acceptance_rate > 0.7) {
    for (i in 1:10) {
      sig.rw <- 2^i * 2.38^2/dd
      out <- rw_mc_cpp(log.target, initial, n_iter = 500, sig=sig.rw)
      if (out$acceptance_rate <= 0.7) break
    }
  } else {
    for (i in 1:10) {
      sig.rw <- 2.38^2/dd/(2^i)
      out <- rw_mc_cpp(log.target, initial, n_iter = 500, sig=sig.rw)
      if (out$acceptance_rate >= 0.25) break
    }
  }
  
  var.base <- sig.rw
  sd.base  <- sqrt(sig.rw)
  return(list(
    sd.base = sd.base,
    var.base = var.base
  ))
}
genfvargmeanvar <- function(log.target, initial, dd) {
  out <- rw_mc_cpp(log.target, initial, n_iter = 500, sig= 2.38^2 /dd)
  if (out$acceptance_rate >= 0.25 && out$acceptance_rate <= 0.7) {
    sig.rw <- 2.38^2 /dd
  } else if (out$acceptance_rate > 0.7) {
    for (i in 1:10) {
      sig.rw <- 2^i * 2.38^2/dd
      out <- rw_mc_cpp(log.target, initial, n_iter = 500, sig=sig.rw)
      if (out$acceptance_rate <= 0.7) break
    }
  } else {
    for (i in 1:10) {
      sig.rw <- 2.38^2/dd/(2^i)
      out <- rw_mc_cpp(log.target, initial, n_iter = 500, sig=sig.rw)
      if (out$acceptance_rate >= 0.25) break
    }
  }
  
  var.base <- sig.rw
  sd.base  <- sqrt(sig.rw)
  
  diag.v.ap<- FALSE
  
  opt <- try(
    optim(initial, fn = log.target, method = "L-BFGS-B",
          control = list(fnscale = -1), hessian = TRUE),
    silent = TRUE
  )
  
  optim_success <- !(inherits(opt, "try-error") || opt$convergence != 0)
  
  if (optim_success) {
    mean.ap.tar <- opt$par
    H <- opt$hessian
    hess_ok <- !any(is.na(H)) && all(is.finite(H))
    
    if (!hess_ok) {
      H <- try(numDeriv::hessian(log.target, mean.ap.tar), silent = TRUE)
    }

    hess_ok <- !(inherits(H, "try-error") || any(is.na(H)) || any(!is.finite(H)))
    
    if (hess_ok) {
      var.ap.tar <- tryCatch({
        V <- -solve(H)
        V <- 0.5 * (V + t(V))
        diag(V) <- pmax(diag(V), 100)
        Vpd <- as.matrix(Matrix::nearPD(V)$mat)
        Vpd
      }, error = function(e) {
        diag.v.ap <<- TRUE
        900 * diag(dd)
      })
    } else {
      var.ap.tar <- 900 * diag(dd)
      diag.v.ap=TRUE
    }
    return(list(
      sd.base = sd.base,
      var.base = var.base,
      mean.ap.tar = mean.ap.tar,
      var.ap.tar = var.ap.tar, 
      diag.v.ap=diag.v.ap
    ))
  }
  
  out1 <- rw_mc_cpp(log.target, out$last_state, n_iter = 4500, sig=sig.rw)
  
  
  mean.ap.tar <- (out$mean * out$n_iter + out1$mean * out1$n_iter) / (out$n_iter + out1$n_iter)
  
  var.ap.tar <- 900 * diag(dd)
  diag.v.ap=TRUE
  
  return(list(
    sd.base = sd.base,
    var.base = var.base,
    mean.ap.tar = mean.ap.tar,
    var.ap.tar = var.ap.tar,
    diag.v.ap=diag.v.ap
  ))
}

checkFuncArgs <- function(func, args) {
  if(!inherits(func, 'function')) {
    return(FALSE)
  }
  func.args <- formals(func)
  if(length(func.args) != length(args)) {
    return(FALSE)
  }
  #  if(!all(names(func.args) == args)) {
  #   return(FALSE)
  # }
  return(TRUE)
}
