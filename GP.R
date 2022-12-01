# convential GP
sepGP <- function(X, y, nu, nug=eps, 
                  init=rep(1, ncol(X)),
                  lower=0.001, upper=10, scale.fg=TRUE, iso.fg=FALSE){
  if(is.null(dim(X))) X <- matrix(X, ncol = 1)

  if(scale.fg){
    X <- scale(X, center = TRUE, scale = TRUE)
    X.center <- attr(X,"scaled:center")
    X.scale <- attr(X,"scaled:scale")
  }else{
    X.center <- rep(0,ncol(X))
    X.scale <- rep(1,ncol(X))
  }

  if(iso.fg) init <- init[1]
  
  n <- length(y)
  nlsep <- function(par, X, Y) 
  {
    theta <- par
    R <- sqrt(distance(t(t(X)/theta)))
    K <- matern.kernel(R, nu=nu)
    Ki <- solve(K+diag(nug,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))
    
    tau2hat <- drop(t(Y - mu.hat) %*% Ki %*% (Y - mu.hat) / n)
    ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
    return(-ll)
  }
  
  outg <- optim(init, nlsep, 
                method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y)
  
  theta <- outg$par
  
  R <- sqrt(distance(t(t(X)/theta)))
  K <- matern.kernel(R, nu=nu)
  Ki <- solve(K+diag(nug,n))
  one.vec <- matrix(1,ncol=1,nrow=n)
  mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
  
  return(list(theta = theta, nu=nu, Ki=Ki, mu.hat=mu.hat,
              nug = nug, X = X, y = y, X.center=X.center, X.scale=X.scale))
}

pred.sepGP <- function(fit, xnew){
  
  xnew <- as.matrix(xnew)
  
  Ki <- fit$Ki
  theta <- fit$theta
  nu <- fit$nu
  nug <- fit$nug
  X <- fit$X
  y <- fit$y
  X.center <- fit$X.center
  X.scale <- fit$X.scale
  mu.hat <- fit$mu.hat
  
  xnew <- t((t(xnew)-X.center)/X.scale)
  
  tau2hat <- drop(t(y-mu.hat) %*% Ki %*% (y-mu.hat) / nrow(X))
  
  RXX <- sqrt(distance(t(t(xnew)/theta)))
  RX <- sqrt(distance(t(t(xnew)/theta), t(t(X)/theta)))
  KXX <- matern.kernel(RXX, nu=nu)
  KX <- matern.kernel(RX, nu=nu)
  
  mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
  Sigmap2 <- pmax(0,diag(tau2hat*(KXX + diag(nug,nrow(xnew)) - KX %*% Ki %*% t(KX))))
  
  return(list(mu=mup2, sig2=Sigmap2))
}
