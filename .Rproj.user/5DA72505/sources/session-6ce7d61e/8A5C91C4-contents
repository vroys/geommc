#' @keywords internal
#' @noRd
###log of multivariate normal density
ldens_mvnorm=function(y,mean,sig){
  pp=length(y)
  R = chol(sig)
  z = backsolve(R,y-mean,transpose=TRUE)
  return(-sum(log(diag(R)))-0.5*sum(z^2)-0.5*pp*log(2 *pi))
}
###log of multivariate normal density for constant variance sig
ldens_mvnorm_diag=function(y,mean,sig){
  pp=length(y)
  return(-0.5*sum((y-mean)^2)/sig-0.5*pp*log(2*pi*sig))
}
###log of multivariate normal density when cholesky factor of sig is given
ldens_mvnormchol=function(y,mean,chol.sig){
  pp=length(y)
  z = backsolve(chol.sig,y-mean,transpose=TRUE)
  return(-sum(log(diag(chol.sig)))-0.5*sum(z^2)-0.5*pp*log(2 *pi))
}
####calculate msjd for MC output## mat is npara times niter matrix
msejd<-function(mat){
  if(is.matrix(mat)==FALSE){
    niter<-length(mat)
    mat.i.1<-mat[2:niter]
    mat.i<-mat[1:(niter-1)]
    dist<-sum((mat.i.1-mat.i)^2)/(niter-1)
  }
  else {
    niter=ncol(mat)
    mat.i.1<-mat[,2:niter]
    mat.i<-mat[,1:(niter-1)]
    dist<-sum((mat.i.1-mat.i)^2)/(niter-1)
  }
  return(dist)
}
bc = function(mu1, mu2, sig1, sig2,TargetOnly){
  ddd=length(mu1)
  if(ddd>1){
    if(TargetOnly){
      sig = 0.5*(sig1[1,1]+sig2[1,1])
       prod=exp(-sum((mu1-mu2)^2)/8/sig-0.5*ddd*(log(sig)-0.5*log(sig1[1,1]*sig2[1,1])))
    }else{
    sig = 0.5*(sig1+sig2)
    R1 = chol(sig1)
    R2 = chol(sig2)
    R = chol(sig)
    z = backsolve(R,mu1-mu2,transpose=T)
    prod=exp(-sum(z^2)/8-(sum(log(diag(R)))-0.5*(sum(log(diag(R1)))+sum(log(diag(R2))))))
      }
    }else{
    prod=exp(-(mu1-mu2)^2/(4*(sig1+sig2))-0.5*(log(sig1+sig2)-(log(2*sqrt(sig1*sig2)))))
  }
  return(unname(cbind(prod,acos(pmin(pmax(prod,0),1.0)))))
}
#cmpute prod_curr= <f(y|curr), g(y|curr)> and \theta_curr
thet=function(k,curr,imp,dens.base,dens.ap.tar,samp.base,samp.ap.tar){
  y =NULL
  curr=curr
  prod=numeric(k)
  if(!imp[1]){
    prod=sapply(seq(1,k),function(i) divonne(function(y) sqrt(dens.ap.tar(y,curr)[i]*dens.base(y,curr)), lowerLimit=rep(-Inf,length(curr)), upperLimit =rep(Inf,length(curr)),absTol = 1e-2)$integral)
  }else{
    if(imp[3]){
      samp.store= matrix(unlist(lapply(seq(1,imp[2]),function(ii) samp.base(curr))), ncol = imp[2])
      imp.store=apply(samp.store,2,function (y) sqrt(dens.ap.tar(y,curr)/dens.base(y,curr)))
      imp.store=imp.store%>%replace(is.infinite(.),0)%>%replace(is.na(.),0)
      if(k>1){
        prod=rowMeans(imp.store)
      }else{
        prod=mean(imp.store)
      }
    }else{
      samp.store=lapply(seq(1,k),function(jj) lapply(seq(1,imp[2]),function(ii) samp.ap.tar(curr,kk=jj)))
      dens.ap=dens.ap.tar(y,curr)
      for(j in 1:k){
        imp.store=apply(matrix(unlist(samp.store[[j]]), ncol = imp[2]),2,function (y) sqrt(dens.base(y,curr)/dens.ap[j]))
        imp.store=imp.store%>%replace(is.infinite(.),0)%>%replace(is.na(.),0)
        prod[j]=mean(imp.store)
      }
    }
  }
  return(unname(cbind(prod,acos(pmin(pmax(prod,0),1.0)))))
}

#h(x|curr) with prod_curr #this returns a vector of length k
dens_h=function(x,curr,prod,dens.base,dens.ap.tar){
  kk=length(prod)
  val=numeric(kk)
  # val[prod==1]=log(dens.base(x,curr))
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
    x=samp_u(curr,prod,kk,samp.base,samp.ap.tar)
    success <- log(runif(1)) < (log(dens_h(x,curr,prod,dens.base,dens.ap.tar)[kk])+log(1-prod[kk]^2)-log(dens.ap.tar(x,curr)[kk]+ prod[kk]^2*dens.base(x,curr)))%>%replace(is.na(.),-Inf)
  }
  return(x)
}

log_phi = function(a,x,curr,prod,theta,eps,dens.base,dens.ap.tar){
  kk=length(a)
  val=numeric(kk)
  if(eps==1){
    val=log(prod^2*dens.base(x,curr)+(1-prod^2)*dens_h(x,curr,prod,dens.base,dens.ap.tar))
  } else{
    coef=(cos(eps*theta))^2
    val=log(coef*dens.base(x,curr)+(1-coef)*dens_h(x,curr,prod,dens.base,dens.ap.tar))
  }
  return(sum(a*val))
}
#sample from phi(|curr)
samp_phi = function(a,curr,prod,theta,eps,dens.base,dens.ap.tar,samp.base,samp.ap.tar){
  kk=length(a)
  s.a=1
  if(kk>1){
    s.a=sample.int(kk,1,prob=a)
  }
  if(prod[s.a]==1){
    return(samp.base(curr))
  }else{
    if(eps==1){
      u = runif(1)
      if(u<=prod[s.a]^2){
        return(samp.base(curr))
      }else
      {return(samp_h(curr,prod,s.a,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
      }
    }else{
      coef=(cos(eps*theta[s.a]))^2
      u = runif(1)
      if(u<=coef){
        return(samp.base(curr))
      }else{
        return(samp_h(curr,prod,s.a,dens.base,dens.ap.tar,samp.base,samp.ap.tar))
      }
    }
  }
}
genfvargmean=function(log.target,initial,dd){
  out<-mcmc::metrop(log.target,initial,nbatch=1,blen=500,scale=2.38/sqrt(dd))
  if(out$accept >= 0.2 && out$accept <= 0.7){
    out1<-mcmc::metrop(out,nbatch=1,blen=4500)
    return(list(sd.base=2.38/sqrt(dd),var.base=2.38^2/dd,mean.ap.tar=c(out$blen*out$batch+out1$blen*out1$batch)/(out$blen+out1$blen)))
  }else if(out$accept>0.7){
    for(i in 1:10){
      sc=2.38/sqrt(dd)/2^i
      out<-mcmc::metrop(log.target,initial,nbatch=1,blen=500,scale=sc)
      if(out$accept <= 0.7){
        out1<-mcmc::metrop(out,nbatch=1,blen=4500,scale=sc)
        return(list(sd.base=sc,var.base=sc^2,mean.ap.tar=c(out$blen*out$batch+out1$blen*out1$batch)/(out$blen+out1$blen)))
      }
    }
    out1<-mcmc::metrop(out,nbatch=1,blen=4500,scale=sc)
    return(list(sd.base=sc,var.base=sc^2,mean.ap.tar=c(out$blen*out$batch+out1$blen*out1$batch)/(out$blen+out1$blen)))
  }else{
    for(i in 1:10){
      sc=2^i*2.38/sqrt(dd)
      out<-mcmc::metrop(log.target,initial,nbatch=1,blen=500,scale=sc)
      if(out$accept >= 0.2){
        out1<-mcmc::metrop(out,nbatch=1,blen=4500,scale=sc)
        return(list(sd.base=sc,var.base=sc^2,mean.ap.tar=(out$blen*out$batch+out1$blen*out1$batch)/(out$blen+out1$blen)))
      }
    }
    out1<-mcmc::metrop(out,nbatch=1,blen=4500,scale=sc)
    return(list(sd.base=sc,var.base=sc^2,mean.ap.tar=c(out$blen*out$batch+out1$blen*out1$batch)/(out$blen+out1$blen)))
  }
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
