KL.expan <- function(d, G, fraction=0.99, rnd=1e3){
  
  X <- sobol(rnd, d)
  n <- length(G)
  Y <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n) Y[,i] <- apply(X,1,G[[i]])
  Y.center <- apply(Y, 1, mean)
  Y <- Y - Y.center
  
  R <- matrix(0,n,n)
  for(i in 1:n) {
    for(j in i:n){
      R[i,j] <- R[j,i] <- sum(Y[,i] * Y[,j])
    }
  }
  eig.out <- eigen(R)

  varphi <- Y %*% eig.out$vectors
  varphi <- t(t(varphi)/apply(varphi, 2, FUN = function(x) sqrt(sum(x^2))))   # normalization
  betai <- t(Y) %*% varphi
  
  # fraction of variance explained
  FoVE <- rep(0,n)
  for(m in 1:n) FoVE[m] <- sum((varphi[,1:m]%*%t(betai[,1:m]))^2)/sum(Y^2)
  M <- which(FoVE > fraction)[1]
  
  return(list(basis=varphi[,1:M], B=betai[,1:M], rnd=rnd,
              rndX=X, Y.center=Y.center, M=M, FoVE=FoVE))
}

KL.Bnew <- function(KL.ls, gnew){
  n <- length(gnew)
  Y <- matrix(0,ncol=n,nrow=KL.ls$rnd)
  for(i in 1:n) Y[,i] <- apply(KL.ls$rndX,1,gnew[[i]])
  Y <- Y - KL.ls$Y.center
  
  return(t(Y) %*% KL.ls$basis)
}
