KL.expan <- function(d, G, fraction=0.99, rnd=1e3){
  
  X <- sobol(rnd, d)
  n <- length(G)
  Y <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n) Y[,i] <- apply(X,1,G[[i]])
  Y.center <- apply(Y, 1, mean)
  Y <- Y - Y.center
  # R <- matrix(0,rnd,rnd)
  # for(i in 1:rnd) {
  #   for(j in i:rnd){
  #     R[i,j] <- R[j,i] <- sum(Y[i,] * Y[j,])
  #   }
  # }
  # eig.out <- eigen(R)
  # 
  # varphi <- eig.out$vectors[,1:M]
  # betai <- t(Y) %*% varphi
  # print(betai)
  
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

# KL.expan <- function(M=10, d, theta, nu=2.5, rnd=1e3){
#  
#   X <- sobol(rnd, d)
#   R <- sqrt(distance(t(t(X)/theta)))
#   Phi <- matern.kernel(R, nu=nu) 
#   
#   eig.out <- eigen(Phi)
#   
#   return(list(basis = eig.out$vector[,1:M] * matrix(sqrt(eig.out$value[1:M]),nrow=rnd,ncol=M,byrow=TRUE),
#               rndX = X))
# }
# 
# optim.theta <- function(G, M, d, nu=2.5, rnd=1e2){
#   
#   theta.fun <- function(theta){
#     KL.ls <- KL.expan(M, d, theta, nu, rnd=rnd)
#     recon.error <- 0
#     for(i in 1:length(G)){
#       lm.fit <- lm(apply(KL.ls$rndX,1,G[[i]]) ~ KL.ls$basis - 1)
#       lm.fit <- summary(lm.fit)
#       recon.error <- recon.error + mean(lm.fit$residuals^2)
#     }
#     return(recon.error)
#   }
#   opt.out <- optim(c(1,1), theta.fun, method="L-BFGS-B", lower=1e-6, upper=100)
# 
#   return(opt.out$par)
# }
# 
# extraEigV <- function(KL.ls, G, M){
#   eigval.mx <- matrix(0,nrow=length(G),ncol=M)
#   for(i in 1:length(G)) eigval.mx[i,] <- coef(lm(apply(KL.ls$rndX,1,G[[i]]) ~ KL.ls$basis - 1))
#   for(i in 1:length(G)) {
#     lm.fit <- lm(apply(KL.ls$rndX,1,G[[i]]) ~ KL.ls$basis - 1)
#     print(mean(lm.fit$residuals^2))
#   }
#   return(eigval.mx)
# }