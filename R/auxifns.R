#' @keywords internal
#' @noRd
get_opt <- function(lst, name) {
if (name %in% names(lst)) lst[[name]] else NULL
}
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
###log of multivariate normal density for constant variance sig
ldens_mvnorm_diag=function(y,mean,sig){
  y    <- as.numeric(y)
  mean <- as.numeric(mean)
  d    <- length(y)
  if (!is.finite(sig) || sig <= 0) return(-Inf)
  sqdist <- sum((y - mean)^2)
  return(-0.5*sqdist/sig-0.5*d*log(2*pi*sig))
}
is_pd <- function(mat) {
  tryCatch({chol(mat); TRUE}, error = function(e) FALSE)
}
####calculate msjd for MC output## mat is npara times niter matrix
msejd<-function(mat){
  if (!is.matrix(mat)) {
    mat <- matrix(mat, ncol = 1)
  }
  niter <- nrow(mat)
  if (niter < 2) return(0) 
  diffs <- mat[2:niter, ,drop = FALSE] - mat[1:(niter - 1), ,drop = FALSE]
  return(sum(diffs * diffs) / (niter - 1))
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
  ind.j=which(prod<1)
  if (length(ind.j) == 0) return(numeric(kk))
  val=numeric(kk)
  sqrt_tar <- sqrt(dens.ap.tar(x,curr))
  sqrt_base <- sqrt(dens.base(x, curr))
  val[ind.j] <- (sqrt_tar[ind.j] - prod[ind.j] * sqrt_base)^2 / (1 - prod[ind.j]^2)
  return(val)
}
#sample from u(.|curr)
samp_u = function(curr,prod,kk=1,samp.base,samp.ap.tar){
  if(runif(1L)<=1/(1+prod[kk]^2)) return(samp.ap.tar(curr,kk))
  else return(samp.base(curr))
}
#sample from h(.|curr)#is called only if prod[kk]<1
samp_h = function(curr,prod,kk=1,dens.base,dens.ap.tar,samp.base,samp.ap.tar){
  max.try = max(200000, 2000*(1+prod^2)/(1-prod^2))
  success <- FALSE
  attempts <- 0
  
  while (!success) {
    attempts <- attempts + 1
    if (attempts > max.try) {
      stop(print(
        "Rejection sampler for h is taking too long to generate a proposal. 
        Try reducing eps, and adjusting other input functions, 
         or use a different starting vector."
      ))
    }
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
  s.a <- if (kk > 1) sample.int(kk, 1L, prob = a) else 1
  if(prod[s.a]==1){
    return(samp.base(curr))
  }
  if (eps == 1) {
    weight <- prod[s.a]^2
  } else {
    weight <- (cos(eps * theta[s.a]))^2
  }
  if (runif(1L) <= weight) {
    return(samp.base(curr))
  } else {
    return(samp_h(curr, prod, s.a,
                  dens.base, dens.ap.tar,
                  samp.base, samp.ap.tar))
  }
}
genfvar <- function(log.target, initial, dd) {
  n_pilot <- pmin(500, pmax(150, ceiling(5000 / dd)))
  
  sig.rw <- 2.38^2 / dd
  
  out <- rwm(log.target, initial, n_iter = n_pilot, sig = sig.rw)
  acc_rate <- out$acceptance_rate
  
  if (acc_rate >= 0.45 && acc_rate <= 0.7) {
    return(list(sd.base = sqrt(sig.rw), var.base = sig.rw))
  }
  
  base_sig <- 2.38^2 / dd
  
  if (acc_rate > 0.7) {
    mult_low <- 1
    mult_high <- 1024  # 2^10
    
    for (iter in 1:10) {
      mult <- exp((log(mult_low) + log(mult_high)) / 2)
      sig.rw <- mult * base_sig
      
      out <- rwm(log.target, initial, n_iter = n_pilot, sig = sig.rw)
      acc_rate <- out$acceptance_rate
      
      if (acc_rate >= 0.45 && acc_rate <= 0.7) break
      
      if (acc_rate > 0.7) {
        mult_low <- mult
      } else {
        mult_high <- mult
      }
    }
    
  } else {
    mult_low <- 1 / 1024  # 2^(-10)
    mult_high <- 1
    
    for (iter in 1:10) {
      mult <- exp((log(mult_low) + log(mult_high)) / 2)
      sig.rw <- mult * base_sig
      
      out <- rwm(log.target, initial, n_iter = n_pilot, sig = sig.rw)
      acc_rate <- out$acceptance_rate
      
      if (acc_rate >= 0.45 && acc_rate <= 0.7) break
      
      if (acc_rate < 0.45) {
        mult_high <- mult
      } else {
        mult_low <- mult
      }
    }
  }
  
  return(list(
    sd.base = sqrt(sig.rw), 
    var.base = sig.rw
  ))
}
genfvargmeanvar <- function(log.target, initial, dd) {
  n_pilot <- pmin(500, pmax(100, ceiling(5000 / dd)))
  
  sig.rw <- 2.38^2 / dd
  
  out <- rwm(log.target, initial, n_iter = n_pilot, sig = sig.rw)
  acc_rate <- out$acceptance_rate
  
  if (acc_rate < 0.45 || acc_rate > 0.7) {
    base_sig <- 2.38^2 / dd
    
    if (acc_rate > 0.7) {
      mult_low <- 1
      mult_high <- 1024
      
      for (iter in 1:10) {
        mult <- exp((log(mult_low) + log(mult_high)) / 2)
        sig.rw <- mult * base_sig
        
        out <- rwm(log.target, initial, n_iter = n_pilot, sig = sig.rw)
        acc_rate <- out$acceptance_rate
        
        if (acc_rate >= 0.45 && acc_rate <= 0.7) break
        
        if (acc_rate > 0.7) {
          mult_low <- mult
        } else {
          mult_high <- mult
        }
      }
    } else {
      mult_low <- 1 / 1024
      mult_high <- 1
      
      for (iter in 1:10) {
        mult <- exp((log(mult_low) + log(mult_high)) / 2)
        sig.rw <- mult * base_sig
        
        out <- rwm(log.target, initial, n_iter = n_pilot, sig = sig.rw)
        acc_rate <- out$acceptance_rate
        
        if (acc_rate >= 0.45 && acc_rate <= 0.7) break
        
        if (acc_rate < 0.45) {
          mult_high <- mult
        } else {
          mult_low <- mult
        }
      }
    }
  }
  
  var.base <- sig.rw
  sd.base <- sqrt(sig.rw)
  
  
  diag.v.ap <- FALSE
  
  start_point <- if (!is.null(out$last_state)) out$last_state else initial
  
  skip_optim <- dd > 100  
  
  if (!skip_optim) {
    opt <- try(
      optim(start_point, fn = log.target, method = "L-BFGS-B",
            control = list(fnscale = -1, maxit = 100 * dd), 
            hessian = TRUE),
      silent = TRUE
    )
    
    optim_success <- !(inherits(opt, "try-error") || opt$convergence != 0)
  } else {
    optim_success <- FALSE
  }
  
  if (optim_success) {
    mean.ap.tar <- opt$par
    H <- opt$hessian
    hess_ok <- !any(is.na(H)) && all(is.finite(H))
    
    if (!hess_ok && dd <= 50) {
      H <- try(numDeriv::hessian(log.target, mean.ap.tar), silent = TRUE)
      hess_ok <- !(inherits(H, "try-error") || any(is.na(H)) || any(!is.finite(H)))
    } else if (!hess_ok) {
      hess_ok <- FALSE  
    }
    
    if (hess_ok) {
      var.ap.tar <- tryCatch({
        V <- -solve(H)
        V <- 0.5 * (V + t(V))
        diag(V) <- pmax(diag(V), 1e-4)
        
        if (dd > 20) {
          max_off_diag <- max(abs(V[row(V) != col(V)]))
          min_diag <- min(abs(diag(V)))
          if (max_off_diag < 0.1 * min_diag) {
            diag.v.ap <<- TRUE
            Vpd <- diag(diag(V))
          } else {
            is.pd <- tryCatch({chol(V); TRUE}, error = function(e) FALSE)
            
            if (is.pd) {
              Vpd <- V
            } else if (dd <= 50) {
              Vpd <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
            } else {
              diag.v.ap <<- TRUE
              Vpd <- diag(pmax(diag(V), 0.01 * mean(abs(diag(V)))))
            }
          }
        } else {
          is.pd <- tryCatch({chol(V); TRUE}, error = function(e) FALSE)
          
          if (is.pd) {
            Vpd <- V
          } else {
            Vpd <- as.matrix(Matrix::nearPD(V, corr = FALSE)$mat)
          }
        }
        Vpd
      }, error = function(e) {
        diag.v.ap <<- TRUE
        diag(rep(400, dd))
      })
    }else {
      diag.v.ap <<- TRUE
      var.ap.tar <- diag(rep(400, dd))
    }
      
      if (!diag.v.ap && all(var.ap.tar == diag(diag(var.ap.tar)))) {
        diag.v.ap <- TRUE
      }
    
      return(list(
        sd.base = sd.base,
        var.base = var.base,
        mean.ap.tar = mean.ap.tar,
        var.ap.tar = var.ap.tar,
        diag.v.ap = diag.v.ap
      ))
    }
  
  total_needed <- pmax(1000, pmin(5000, ceiling(50000 / dd)))
  remaining_iter <- total_needed - n_pilot
  
  if (remaining_iter > 0) {
    out1 <- rwm(log.target, out$last_state, 
                      n_iter = remaining_iter, sig = sig.rw)
    
    mean.ap.tar <- (out$mean * n_pilot + out1$mean * remaining_iter) / total_needed
  } else {
    mean.ap.tar <- out$mean
  }
  
  var.ap.tar <- diag(rep(400, dd))
  diag.v.ap <- TRUE
  
  return(list(
    sd.base = sd.base,
    var.base = var.base,
    mean.ap.tar = mean.ap.tar,
    var.ap.tar = var.ap.tar,
    diag.v.ap = diag.v.ap
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
