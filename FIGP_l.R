# FIGP with linear kernel
FIGP_l <- function(G, d, y, nu, nug, rnd=1e5, 
                  init=0.01, lower=1e-6, upper=0.1){
                  # init=0.1, lower=1e-6, upper=5){
  
  n <- length(y)
  nlsep <- function(par, G, d, Y, rnd) 
  {
    theta <- par
    
    n <- length(Y)
    X <- sobol(rnd, d)
    R <- sqrt(distance(t(t(X)/theta)))
    Phi <- matern.kernel(R, nu=nu)
    
    A <- matrix(0,ncol=n,nrow=rnd)
    for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
    K <- t(A) %*% Phi %*% A

    Ki <- solve(K+diag(nug,n))
    ldetK <- determinant(K+diag(nug,n), logarithm=TRUE)$modulus
    
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
    ll <- - (n/2)*log(t(Y-mu.hat) %*% Ki %*% (Y-mu.hat)) - (1/2)*ldetK
    return(-ll)
  }
  
  tic <- proc.time()[3]
  
  outg <- optim(c(rep(init, d)), nlsep,
                method="L-BFGS-B", lower=lower, upper=upper, G=G, Y=y, d=d, rnd=rnd)
  toc <- proc.time()[3]
  
  theta <- outg$par
  
  X <- sobol(rnd, d)
  R <- sqrt(distance(t(t(X)/theta)))
  Phi <- matern.kernel(R, nu=nu)
  
  A <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
  K <- t(A) %*% Phi %*% A
  Ki <- solve(K+diag(nug,n))
  one.vec <- matrix(1,ncol=1,nrow=n)
  mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
  
  return(list(theta = theta, nu=nu, Ki=Ki, Phi=Phi, A=A,
              nug = nug, G = G, y = y, X = X, rnd = rnd, mu.hat=mu.hat))
}

pred.FIGP_l <- function(fit, gnew, sig2.fg=TRUE){
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  nug <- fit$nug
  G <- fit$G
  X <- fit$X
  y <- fit$y
  Phi <- fit$Phi
  A <- fit$A
  rnd <- fit$rnd
  mu.hat <- fit$mu.hat
  
  n <- length(y)
  n.new <- length(gnew)
  
  # R <- sqrt(distance(t(t(X)/theta)))
  # Phi <- matern.kernel(R, nu=nu)
  # A <- matrix(0,ncol=n,nrow=rnd)
  # for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
  a <- matrix(0,ncol=n.new,nrow=rnd)
  for(i in 1:n.new)  a[,i] <- apply(X, 1, gnew[[i]])
  KXX <- t(a) %*% Phi %*% a
  KX <- t(a) %*% Phi %*% A
  
  mup2 <- drop(mu.hat + KX %*% Ki %*% (y - mu.hat))
  if(sig2.fg){
    tau2hat <- drop(t(y - mu.hat) %*% Ki %*% (y - mu.hat) / n)
    Sigmap2 <- pmax(0,diag(tau2hat*(KXX + diag(nug, n.new) - KX %*% Ki %*% t(KX))))
  }else{
    Sigmap2 <- NULL
  }
  
  return(list(mu=mup2, sig2=Sigmap2))
}
