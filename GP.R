# convential GP
sepGP <- function(X, y, nu, g=eps, 
                  init=rep(1, ncol(X)),
                  lower=0.001, upper=100){
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)
  
  X <- scale(X)
  X.center <- attr(X,"scaled:center")
  X.scale <- attr(X,"scaled:scale")
  
  n <- length(y)
  nlsep <- function(par, X, Y) 
  {
    theta <- par
    R <- sqrt(distance(t(t(X)/theta)))
    K <- matern.kernel(R, nu=nu)
    Ki <- solve(K+diag(g,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
    return(-ll)
  }
  
  outg <- optim(init, nlsep, 
                method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y)
  
  theta <- outg$par
  
  R <- sqrt(distance(t(t(X)/theta)))
  K <- matern.kernel(R, nu=nu)
  Ki <- solve(K+diag(g,n))
  one.vec <- matrix(1,ncol=1,nrow=n)
  mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
  
  return(list(theta = theta, nu=nu, Ki=Ki, mu.hat=mu.hat,
              g = g, X = X, y = y, X.center=X.center, X.scale=X.scale))
}

pred.sepGP <- function(fit, xnew){
  
  xnew <- as.matrix(xnew)
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  g <- fit$g
  X <- fit$X
  y <- fit$y
  X.center <- fit$X.center
  X.scale <- fit$X.scale
  mu.hat <- fit$mu.hat
  
  xnew <- t((t(xnew)-X.center)/X.scale)
  
  tau2hat <- drop(t(y) %*% Ki %*% y / nrow(X))
  
  RXX <- sqrt(distance(t(t(xnew)/theta)))
  RX <- sqrt(distance(t(t(xnew)/theta), t(t(X)/theta)))
  KXX <- matern.kernel(RXX, nu=nu)
  KX <- matern.kernel(RX, nu=nu)
  
  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- pmax(0,diag(tau2hat*(KXX + diag(g,nrow(xnew)) - KX %*% Ki %*% t(KX))))
  
  return(list(mu=mup2, sig2=Sigmap2))
}
