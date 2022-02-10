# FIGP with nonlinear kernel
FIGP_nl <- function(G, d, y, nu, nug, rnd=1e5,
                     init=1, lower=1e-2, upper=100){
  
  n <- length(y)
  nlsep <- function(par, G, d, Y, rnd) 
  {
    theta <- par
    
    n <- length(Y)
    X <- sobol(rnd, d)
    
    A <- matrix(0,ncol=n,nrow=rnd)
    for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])

    R <- sqrt(distance(t(A)/theta)/rnd)
    K <- matern.kernel(R, nu=nu)
  
    Ki <- solve(K+diag(nug,n))
    ldetK <- determinant(K+diag(nug,n), logarithm=TRUE)$modulus
    
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
    ll <- - (n/2)*log(t(Y-mu.hat) %*% Ki %*% (Y-mu.hat)) - (1/2)*ldetK
    return(-ll)
  }
  
  tic <- proc.time()[3]
  
  outg <- optim(init, nlsep,
                method="L-BFGS-B", lower=lower, upper=upper, G=G, Y=y, d=d, rnd=rnd)
  toc <- proc.time()[3]
  
  theta <- outg$par
  
  X <- matrix(runif(rnd*d), ncol=d)
  A <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
  
  R <- sqrt(distance(t(A)/theta)/rnd)
  K <- matern.kernel(R, nu=nu)
  
  Ki <- solve(K+diag(nug,n))
  one.vec <- matrix(1,ncol=1,nrow=n)
  mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
  
  return(list(theta = theta, nu=nu, Ki=Ki, A=A,
              nug = nug, G = G, y = y, X = X, rnd = rnd, mu.hat=mu.hat))
}

pred.FIGP_nl <- function(fit, gnew, sig2.fg=TRUE){
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  nug <- fit$nug
  G <- fit$G
  X <- fit$X
  y <- fit$y
  A <- fit$A
  rnd <- fit$rnd
  mu.hat <- fit$mu.hat
  
  n <- length(y)
  n.new <- length(gnew)
  

  # A <- matrix(0,ncol=n,nrow=rnd)
  # for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
  a <- matrix(0,ncol=n.new,nrow=rnd)
  for(i in 1:n.new)  a[,i] <- apply(X, 1, gnew[[i]])
  
  RX <- sqrt(distance(t(a)/theta,t(A)/theta)/rnd)
  KX <- matern.kernel(RX, nu=nu)
  
  RXX <- sqrt(distance(t(a)/theta)/rnd)
  KXX <- matern.kernel(RXX, nu=nu)
  
  mup2 <- drop(mu.hat + KX %*% Ki %*% (y - mu.hat))
  if(sig2.fg){
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / n)
    Sigmap2 <- pmax(0,diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))
  }else{
    Sigmap2 <- NULL
  }

  
  return(list(mu=mup2, sig2=Sigmap2))
}
