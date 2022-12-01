FIGP <- function(G, d, y, nu, nug, 
                 kernel=c("linear", "nonlinear")[1],
                 theta.init=ifelse(kernel=="linear", 0.01, 1),
                 theta.lower=ifelse(kernel=="linear", 1e-6, 1e-2),
                 theta.upper=ifelse(kernel=="linear", 0.1, 100)){
  
  n <- length(y)
  nlsep <- function(par, G, d, Y, rnd) 
  {
    theta <- par
    
    n <- length(Y)
    
    K <- FIGP.kernel(d, theta, nu, G, kernel=kernel)
    Ki <- solve(K+diag(nug,n))
    ldetK <- determinant(K+diag(nug,n), logarithm=TRUE)$modulus
    
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
  
    tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
    ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
    return(-ll)
  }
  
  tic <- proc.time()[3]
  
  if(kernel=="linear"){
    out <- optim(c(rep(theta.init, d)), nlsep,
                  method="L-BFGS-B", lower=theta.lower, upper=theta.upper, G=G, Y=y, d=d, rnd=rnd)
    theta <- out$par
  }else{
    out <- optimize(nlsep, lower=theta.lower, upper=theta.upper, G=G, Y=y, d=d, rnd=rnd)
    theta <- out$minimum
  }
  toc <- proc.time()[3]
  
  K <- FIGP.kernel(d, theta, nu, G, kernel=kernel)
  Ki <- solve(K+diag(nug,n))
  one.vec <- matrix(1,ncol=1,nrow=n)
  mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
  
  return(list(theta = theta, nu=nu, Ki=Ki, d=d, kernel=kernel, ElapsedTime=toc-tic, 
              nug = nug, G = G, y = y, mu.hat = mu.hat))
}

pred.FIGP <- function(fit, gnew, sig2.fg=TRUE){
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  nug <- fit$nug
  G <- fit$G
  y <- fit$y
  d <- fit$d
  mu.hat <- fit$mu.hat
  kernel <- fit$kernel
  
  n <- length(y)
  n.new <- length(gnew)
  
  KX <- FIGP.kernel(d, theta, nu, G, gnew, kernel=kernel)
  KXX <- FIGP.kernel(d, theta, nu, gnew, kernel=kernel)
  
  mup2 <- drop(mu.hat + KX %*% Ki %*% (y - mu.hat))
  if(sig2.fg){
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / n)
    Sigmap2 <- pmax(0,diag(tau2hat*(KXX + diag(nug, n.new) - KX %*% Ki %*% t(KX))))
  }else{
    Sigmap2 <- NULL
  }
  
  return(list(mu=mup2, sig2=Sigmap2))
}
