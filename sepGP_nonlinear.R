sepGP_nl <- function(G, d, y, nu, nug, rnd=1e5){
  
  n <- length(y)
  nlsep <- function(par, G, d, Y, rnd) 
  {
    theta <- par
    
    n <- length(Y)
    X <- sobol(rnd, d)
    
    A <- matrix(0,ncol=n,nrow=rnd)
    for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])

    R <- sqrt(distance(t(A)/theta))
    K <- matern.kernel(R, nu=nu)
  
    Ki <- solve(K+diag(nug,n))
    ldetK <- determinant(K+diag(nug,n), logarithm=TRUE)$modulus
    
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
    ll <- - (n/2)*log(t(Y-mu.hat) %*% Ki %*% (Y-mu.hat)) - (1/2)*ldetK
    return(-ll)
  }
  
  tic <- proc.time()[3]
  
  outg <- optim(1e2, nlsep,
                method="L-BFGS-B", lower=1, upper=1e5, G=G, Y=y, d=d, rnd=rnd)
  toc <- proc.time()[3]
  
  theta <- outg$par
  
  X <- matrix(runif(rnd*d), ncol=d)
  A <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
  
  R <- sqrt(distance(t(A)/theta))
  K <- matern.kernel(R, nu=nu)
  
  Ki <- solve(K+diag(nug,n))
  one.vec <- matrix(1,ncol=1,nrow=n)
  mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
  
  return(list(theta = theta, nu=nu, Ki=Ki,
              nug = nug, G = G, y = y, X = X, rnd = rnd, mu.hat=mu.hat))
}

pred.sepGP_nl <- function(fit, gnew){
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  nug <- fit$nug
  G <- fit$G
  X <- fit$X
  y <- fit$y
  rnd <- fit$rnd
  mu.hat <- fit$mu.hat
  
  n <- length(y)
  n.new <- length(gnew)
  

  A <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
  a <- matrix(0,ncol=n.new,nrow=rnd)
  for(i in 1:n.new)  a[,i] <- apply(X, 1, gnew[[i]])
  
  RX <- sqrt(distance(t(a)/theta,t(A)/theta))
  KX <- matern.kernel(RX, nu=nu)
  
  RXX <- sqrt(distance(t(a)/theta))
  KXX <- matern.kernel(RXX, nu=nu)
  
  mup2 <- drop(mu.hat + KX %*% Ki %*% (y - mu.hat))
  tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / n)
  Sigmap2 <- pmax(0,diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))
  
  return(list(mu=mup2, sig2=Sigmap2))
}
