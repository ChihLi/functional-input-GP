library(randtoolbox)
library(R.matlab)
library(cubature)
library(plgp)
source("FIGP_l.R")              # FIGP with linear kernel
source("FIGP_nl.R")             # FIGP with non-linear kernel
source("matern.kernel.R")       # matern kernel computation
source("loocv.R")               # LOOCV for FIGP
source("KL.expan.R")            # KL expansion for comparison
source("GP.R")                  # conventional GP

##### Section 4 #####
cat("Section 4...\n")
set.seed(1) #set a random seed for reproducing
eps <- sqrt(.Machine$double.eps) #small nugget for numeric stability

kernel.linear <- function(nu, theta, rnd=1000){
  x <- seq(0,2*pi,length.out = rnd)
  R <- sqrt(distance(x*theta))
  Phi <- matern.kernel(R, nu=nu)
  a <- seq(0,1,0.01)
  n <- length(a)
  A <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n)  A[,i] <- sin(a[i]*x)
  K <- t(A) %*% Phi %*% A
  return(K)
}
kernel.nonlinear <- function(nu, theta, rnd=1000){
  x <- seq(0,2*pi,length.out = rnd)
  a <- seq(0,1,0.01)
  n <- length(a)
  A <- matrix(0,ncol=n,nrow=rnd)
  for(i in 1:n)  A[,i] <- sin(a[i]*x)
  R <- sqrt(distance(t(A)*theta)/rnd)
  
  K <- matern.kernel(R, nu=nu)
  return(K)
}

theta <- 1
s2 <- 1
nu <- c(0.5,3,10)
K1 <- kernel.linear(nu=nu[1], theta=theta)
K2 <- kernel.linear(nu=nu[2], theta=theta) 
K3 <- kernel.linear(nu=nu[3], theta=theta) 

cat("   reproducing Figure 2...\n")
pdf("sample_path_linear.pdf", width = 7, height = 7)
par(mfrow=c(3,3), mar = c(4, 4, 2, 1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(nu==1/2))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(nu==3))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, xlab=expression(alpha), 
        ylab="y", main=expression(nu==10))

nu <- 2.5
theta <- c(0.01,1,100)
s2 <- 1
K1 <- kernel.linear(nu=nu, theta=theta[1])
K2 <- kernel.linear(nu=nu, theta=theta[2]) 
K3 <- kernel.linear(nu=nu, theta=theta[3])
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(theta==0.01))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(theta==1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, xlab=expression(alpha), 
        ylab="y", main=expression(theta==100))

nu <- 2.5
theta <- 1
s2 <- c(0.1,1,100)
K1 <- kernel.linear(nu=nu, theta=theta)
K2 <- kernel.linear(nu=nu, theta=theta) 
K3 <- kernel.linear(nu=nu, theta=theta) 

matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[1]*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==0.1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[2]*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[3]*K3)), type="l", col=3, lty=3, xlab=expression(alpha), 
        ylab="y", main=expression(sigma^2==100))
dev.off()

cat("   reproducing Figure 3...\n")
pdf("sample_path_nonlinear.pdf", width = 7, height = 7)
gamma <- 1
s2 <- 1
nu <- c(0.5,2,10)
K1 <- kernel.nonlinear(nu=nu[1], theta=gamma)
K2 <- kernel.nonlinear(nu=nu[2], theta=gamma) 
K3 <- kernel.nonlinear(nu=nu[3], theta=gamma) 

par(mfrow=c(3,3), mar = c(4, 4, 2, 1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(nu==1/2))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(nu==2))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, xlab=expression(alpha), 
        ylab="y", main=expression(nu==10))

nu <- 2.5
gamma <- c(0.1,1,10)
s2 <- 1
K1 <- kernel.nonlinear(nu=nu, theta=gamma[1])
K2 <- kernel.nonlinear(nu=nu, theta=gamma[2]) 
K3 <- kernel.nonlinear(nu=nu, theta=gamma[3])
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(gamma==0.1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(gamma==1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, xlab=expression(alpha), 
        ylab="y", main=expression(gamma==10))

nu <- 2.5
gamma <- 1
s2 <- c(0.1,1,100)
K1 <- kernel.nonlinear(nu=nu, theta=gamma)
K2 <- kernel.nonlinear(nu=nu, theta=gamma) 
K3 <- kernel.nonlinear(nu=nu, theta=gamma) 

matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[1]*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==0.1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[2]*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[3]*K3)), type="l", col=3, lty=3, xlab=expression(alpha), 
        ylab="y", main=expression(sigma^2==100))
dev.off()

# training functional inputs (G)
G <- list(function(x) x[1]+x[2],
          function(x) x[1]^2,
          function(x) x[2]^2,
          function(x) 1+x[1],
          function(x) 1+x[2],
          function(x) 1+x[1]*x[2],
          function(x) sin(x[1]),
          function(x) cos(x[1]+x[2]))
n <- length(G)
# y1: integrate g function from 0 to 1
y1 <- rep(0, n) 
for(i in 1:n) y1[i] <- hcubature(G[[i]], lower=c(0, 0),upper=c(1,1))$integral

# y2: integrate g^3 function from 0 to 1
G.cubic <- list(function(x) (x[1]+x[2])^3,
                 function(x) (x[1]^2)^3,
                 function(x) (x[2]^2)^3,
                 function(x) (1+x[1])^3,
                 function(x) (1+x[2])^3,
                 function(x) (1+x[1]*x[2])^3,
                 function(x) (sin(x[1]))^3,
                 function(x) (cos(x[1]+x[2]))^3)
y2 <- rep(0, n) 
for(i in 1:n) y2[i] <- hcubature(G.cubic[[i]], lower=c(0, 0),upper=c(1,1))$integral

# y3: integrate sin(g^2) function from 0 to 1
G.sin <- list(function(x) sin((x[1]+x[2])^2),
              function(x) sin((x[1]^2)^2),
              function(x) sin((x[2]^2)^2),
              function(x) sin((1+x[1])^2),
              function(x) sin((1+x[2])^2),
              function(x) sin((1+x[1]*x[2])^2),
              function(x) sin((sin(x[1]))^2),
              function(x) sin((cos(x[1]+x[2]))^2))
y3 <- rep(0, n) 
for(i in 1:n) y3[i] <- hcubature(G.sin[[i]], lower=c(0, 0),upper=c(1,1))$integral

cat("   reproducing Table 1...\n")
Y <- cbind(y1,y2,y3)
print(Y)

cat("   reproducing Table 2...\n")
# model fitting
loocv.l <- loocv.nl <- rep(0,3)
gp.fit <- gpnl.fit <- vector("list", 3)
set.seed(1)
for(i in 1:3){
  gp.fit[[i]] <- FIGP_l(G, d=2, Y[,i], nu=2.5, nug=eps, rnd=1e3)
  loocv.l[i] <- loocv(gp.fit[[i]])
  
  gpnl.fit[[i]] <- FIGP_nl(G, d=2, Y[,i], nu=2.5, nug=eps, rnd=1e3)
  loocv.nl[i] <- loocv(gpnl.fit[[i]])
}

# for comparison: basis expansion approach
set.seed(1)
KL.out <- KL.expan(d=2, G, fraction=0.99, rnd=1e3)
B <- KL.out$B
  
KL.fit <- vector("list", KL.out$M)
for(i in 1:KL.out$M) KL.fit[[i]] <- sepGP(B, Y[,i], nu=2.5, g=eps)
  
# run predictions
set.seed(1)
n.test <- 100

alpha1 <- runif(n.test,0,1)
alpha2 <- runif(n.test,0,1)
beta1 <- runif(n.test,0,1)
kappa1 <- runif(n.test,0,1)

mse.linear <- mse.nonlinear <- mse.kl <- 
  cvr.linear <- cvr.nonlinear <- cvr.kl <- 
  score.linear <- score.nonlinear <- score.kl <- rep(0,3)

# scoring rule function
score <- function(x, mu, sig2){
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -(x-mu)^2/sig2-log(sig2)
}

for(i in 1:3){
  mse.linear.i <- mse.nonlinear.i <- mse.kl.i <- 
    cvr.linear.i <- cvr.nonlinear.i <- cvr.kl.i <- 
    score.linear.i <- score.nonlinear.i <- score.kl.i <- rep(0, n.test)
  for(ii in 1:n.test){
    gnew <- list(function(x) sin(alpha1[ii]*x[1]+alpha2[ii]*x[2]),
                 function(x) beta1[ii]+x[1]^2+x[2]^3,
                 function(x) exp(-kappa1[ii]*x[1]*x[2]))    
    if(i==1){
      g.int <- gnew
    }else if(i==2){
      g.int <- list(function(x) (sin(alpha1[ii]*x[1]+alpha2[ii]*x[2]))^3,
                    function(x) (beta1[ii]+x[1]^2+x[2]^3)^3,
                    function(x) (exp(-kappa1[ii]*x[1]*x[2]))^3)
    }else if(i==3){
      g.int <- list(function(x) sin((sin(alpha1[ii]*x[1]+alpha2[ii]*x[2]))^2),
                    function(x) sin((beta1[ii]+x[1]^2+x[2]^3)^2),
                    function(x) sin((exp(-kappa1[ii]*x[1]*x[2]))^2))
    }
    
    n.new <- length(gnew)
    y.true <- rep(0,n.new)
    for(iii in 1:n.new) y.true[iii] <- hcubature(g.int[[iii]], lower=c(0, 0),upper=c(1,1))$integral
    
    ynew <- pred.FIGP_l(gp.fit[[i]], gnew)
    mse.linear.i[ii] <- mean((y.true - ynew$mu)^2)
    lb <- ynew$mu - qnorm(0.975)*sqrt(ynew$sig2)
    ub <- ynew$mu + qnorm(0.975)*sqrt(ynew$sig2)
    cvr.linear.i[ii] <- mean(y.true > lb & y.true < ub)
    score.linear.i[ii] <- mean(score(y.true, ynew$mu, ynew$sig2))
    
    ynew <- pred.FIGP_nl(gpnl.fit[[i]], gnew)
    mse.nonlinear.i[ii] <- mean((y.true - ynew$mu)^2)
    lb <- ynew$mu - qnorm(0.975)*sqrt(ynew$sig2)
    ub <- ynew$mu + qnorm(0.975)*sqrt(ynew$sig2)
    cvr.nonlinear.i[ii] <- mean(y.true > lb & y.true < ub)
    score.nonlinear.i[ii] <- mean(score(y.true, ynew$mu, ynew$sig2))
    
    B.new <- KL.Bnew(KL.out, gnew)
    ynew <- pred.sepGP(KL.fit[[i]], B.new)
    mse.kl.i[ii] <- mean((y.true - ynew$mu)^2)
    lb <- ynew$mu - qnorm(0.975)*sqrt(ynew$sig2)
    ub <- ynew$mu + qnorm(0.975)*sqrt(ynew$sig2)
    cvr.kl.i[ii] <- mean(y.true > lb & y.true < ub)
    score.kl.i[ii] <- mean(score(y.true, ynew$mu, ynew$sig2))
  }
  mse.linear[i] <- mean(mse.linear.i)*100
  mse.nonlinear[i] <- mean(mse.nonlinear.i)*100
  mse.kl[i] <- mean(mse.kl.i)*100
  cvr.linear[i] <- mean(cvr.linear.i)*100
  cvr.nonlinear[i] <- mean(cvr.nonlinear.i)*100
  cvr.kl[i] <- mean(cvr.kl.i)*100
  score.linear[i] <- mean(score.linear.i)
  score.nonlinear[i] <- mean(score.nonlinear.i)
  score.kl[i] <- mean(score.kl.i)
}

out <- rbind(format(loocv.l,digits=4),
             format(loocv.nl,digits=4),
             format(mse.linear,digits=4),
             format(mse.nonlinear,digits=4))
rownames(out) <- c("linear LOOCV", "nonlinear LOOCV", "linear MSE", "nonlinear MSE")
colnames(out) <- c("y1", "y2", "y3")
print(out)

select.idx <- apply(rbind(loocv.l, loocv.nl), 2, which.min)
select.mse <- diag(rbind(mse.linear, mse.nonlinear)[select.idx,])
select.cvr <- diag(rbind(cvr.linear, cvr.nonlinear)[select.idx,])
select.score <- diag(rbind(score.linear, score.nonlinear)[select.idx,])

cat("   reproducing Table 3...\n")
out <- rbind(format(select.mse,digits=4),
             format(mse.kl,digits=4),
             format(select.cvr,digits=4),
             format(cvr.kl,digits=4),
             format(select.score,digits=4),
             format(score.kl,digits=4))
rownames(out) <- c("FIGP MSE", "Basis MSE", "FIGP coverage", "Basis coverage", "FIGP score", "Basis score")
colnames(out) <- c("y1", "y2", "y3")
print(out)

##### Section 5 #####
cat("Section 5...\n")
func.title <- c("g(x1,x2)=1+x1","g(x1,x2)=1-x1","g(x1,x2)=1+x1x2","g(x1,x2)=1-x1x2",
                "g(x1,x2)=1+x2","g(x1,x2)=1-x2","g(x1,x2)=1+x1^2","g(x1,x2)=1-x1^2",
                "g(x1,x2)=1+x2^2","g(x1,x2)=1-x2^2")

cat("   reproducing Figure 4...\n")
output.mx <- matrix(0,nrow=10,ncol=32*32)
pdf("realcase_real.pdf", width = 12, height = 5)
par(mfrow=c(2,5))
par(mar = c(1, 1, 2, 1))
for(i in 1:10){
  g.out <- readMat(paste0("DATA/q_func",i,".mat"))$Ffem
  image(Re(g.out), zlim=c(0.064,0.103),yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main=func.title[i])
  contour(Re(g.out), add = TRUE, nlevels = 5)
  output.mx[i,] <- c(Re(g.out))
}
dev.off()

pca.out <- prcomp(output.mx, scale = TRUE, center = TRUE)
n.comp <- which(summary(pca.out)$importance[3,] > 0.999)[1]
print(n.comp)

cat("   reproducing Figure 5...\n")
pdf("realcase_pc.pdf", width = 7, height = 2.5)
par(mfrow=c(1,3))
par(mar = c(1, 1, 2, 1))
for(i in 1:n.comp){
  eigen.vec <- matrix(c(pca.out$rotation[,i]), 32, 32)
  image(eigen.vec,yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main=paste("PC",i))
  contour(eigen.vec, add = TRUE, nlevels = 5)
}
dev.off()

# training functional inputs (G)
G <- list(function(x) 1+x[1],
          function(x) 1-x[1],
          function(x) 1+x[1]*x[2],
          function(x) 1-x[1]*x[2],
          function(x) 1+x[2],
          function(x) 1-x[2],
          function(x) 1+x[1]^2,
          function(x) 1-x[1]^2,
          function(x) 1+x[2]^2,
          function(x) 1-x[2]^2)
n <- length(G)

set.seed(1)
gp.fit <- gpnl.fit <- vector("list",n.comp)
for(i in 1:n.comp){
  y <- pca.out$x[,i]
  # fit FIGP with a linear kernel  
  gp.fit[[i]] <- FIGP_l(G, d=2, y, nu=2.5, nug=eps, rnd=1e3)
  # fit FIGP with a nonlinear kernel    
  gpnl.fit[[i]] <- FIGP_nl(G, d=2, y, nu=2.5, nug=eps, rnd=1e3)
}

loocv.recon <- cbind(loocv.pred(gp.fit[[1]]), loocv.pred(gp.fit[[2]]), loocv.pred(gp.fit[[3]])) %*% 
  t(pca.out$rotation[,1:n.comp])
loocv.recon <- scale(loocv.recon, center = FALSE, scale = 1/pca.out$scale)
loocv.recon <- scale(loocv.recon, scale = FALSE, center = -1 * pca.out$center)
loocv.linear <- mean((loocv.recon - output.mx)^2)

loocv.nl.recon <- cbind(loocv.pred(gpnl.fit[[1]]), loocv.pred(gpnl.fit[[2]]), loocv.pred(gpnl.fit[[3]])) %*% 
  t(pca.out$rotation[,1:n.comp])
loocv.nl.recon <- scale(loocv.nl.recon, center = FALSE, scale = 1/pca.out$scale)
loocv.nl.recon <- scale(loocv.nl.recon, scale = FALSE, center = -1 * pca.out$center)
loocv.nonlinear <- mean((loocv.nl.recon - output.mx)^2)

cat("   computing LOOCVs...\n")
out <- c(loocv.linear, loocv.nonlinear)
names(out) <- c("linear", "nonlinear")
print(out)

# test functional inputs (gnew)
gnew <- list(function(x) 1-sin(x[2]))
n.new <- length(gnew)

# make predictions using a linear kernel
ynew <-  matrix(0,ncol=n.comp,nrow=n.new)
for(i in 1:n.comp){
  y <- pca.out$x[,i]
  # fit FIGP with a linear kernel  
  gp.fit[[i]] <- FIGP_l(G, d=2, y, nu=2.5, nug=eps, rnd=1e3)
  pred.out <- pred.FIGP_l(gp.fit[[i]], gnew)
  ynew[,i] <- pred.out$mu
}

# reconstruct the image
pred.recon <- ynew %*% t(pca.out$rotation[,1:n.comp])
pred.recon <- scale(pred.recon, center = FALSE, scale = 1/pca.out$scale)
pred.recon <- scale(pred.recon, scale = FALSE, center = -1 * pca.out$center)

# basis function method for comparison
KL.out <- KL.expan(d=2, G, fraction=0.99, rnd=1e3)
B <- KL.out$B
B.new <- KL.Bnew(KL.out, gnew)

ynew <-  matrix(0,ncol=n.comp,nrow=n.new)
KL.fit <- vector("list", 3)
for(i in 1:n.comp){
  KL.fit[[i]] <- sepGP(B, pca.out$x[,i], nu=2.5, g=eps)
  pred.out <- pred.sepGP(KL.fit[[i]], B.new)
  ynew[,i] <- drop(pred.out$mu)
}

# reconstruct the image
pred.KL.recon <- ynew %*% t(pca.out$rotation[,1:n.comp])
pred.KL.recon <- scale(pred.KL.recon, center = FALSE, scale = 1/pca.out$scale)
pred.KL.recon <- scale(pred.KL.recon, scale = FALSE, center = -1 * pca.out$center)


# true data on the test data
gnew.true <- matrix(0, ncol=n.new, nrow=32*32)
gnew.dat <- readMat(paste0("DATA/q_sine.mat"))$Ffem
gnew.true[,1] <- c(Re(gnew.dat))


# plot the result
cat("   reproducing Figure 6...\n")
pdf("realcase_prediction.pdf", width = 6, height = 2.3)

par(mfrow=c(1,3))
par(mar = c(1, 1, 2, 1))

mse <- mse.kl <- rep(0, n.new)
for(i in 1:n.new){
  image(matrix(gnew.true[,i],32,32), zlim=c(0.064,0.11),yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main=ifelse(i==1, "g(x1,x2)=1-sin(x2)", "g(x1,x2)=1"))
  contour(matrix(gnew.true[,i],32,32), add = TRUE, nlevels = 5)

  image(matrix(pred.recon[i,], 32, 32), zlim=c(0.064,0.11),yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main="FIGP prediction")
  contour(matrix(pred.recon[i,], 32, 32), add = TRUE, nlevels = 5)
  mse[i] <- mean((gnew.true[,i]-pred.recon[i,])^2)*100

  image(matrix(pred.KL.recon[i,], 32, 32), zlim=c(0.064,0.11),yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main="Basis prediction")
  contour(matrix(pred.KL.recon[i,], 32, 32), add = TRUE, nlevels = 5)
  mse.kl[i] <- mean((gnew.true[,i]-pred.KL.recon[i,])^2)*100
}
dev.off()
cat("   computing MSE...\n")

out <- cbind(mse, mse.kl)
colnames(out) <- c("FIGP", "Basis")
print(out)
