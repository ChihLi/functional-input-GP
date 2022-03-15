Functional-Input Gaussian Processes with Applications to Inverse
Scattering Problems (Reproducibility)
================
Chih-Li Sung
March 15, 2022

This instruction aims to reproduce the results in the paper
“*Functional-Input Gaussian Processes with Applications to Inverse
Scattering Problems*” by Sung et al. ([https://arxiv.org/abs/2201.01682](https://arxiv.org/abs/2201.01682)).  Hereafter, functional-Input
Gaussian Process is abbreviated by *FIGP*.

The following results are reproduced in this file

-   The sample path plots in Section 4.1 (Figures 2 and 3)
-   The prediction results in Section 4.2 (Tables 1, 2, and 3)
-   The plots and prediction results in Section 5 (Figures 4, 5, and 6)

##### Step 0.1: load functions and packages

``` r
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
```

##### Step 0.2: setting

``` r
set.seed(1) #set a random seed for reproducing
eps <- sqrt(.Machine$double.eps) #small nugget for numeric stability
```

## Reproducing Section 4.1: Sample Path

Set up the kernel functions introduced in Section 3. `kernel.linear` is
the linear kernel in Section 3.1, while `kernel.nonlinear` is the
non-linear kernel in Section 3.2.

``` r
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
```

##### Reproducing Figure 2

Consider a linear kernel with various choices of parameter settings,
including `nu`, `theta`, `s2`.

-   First row: Set `theta=1` and `s2=1` and set different values for
    `nu`, which are 0.5, 3, and 10.
-   Second row: Set `nu=2.5` and `s2=1` and set different values for
    `theta`, which are 0.01, 1, and 100.
-   Third row: Set `nu=2.5` and `theta=1` and set different values for
    `s2`, which are 0.01, 1, and 100.

``` r
theta <- 1
s2 <- 1
nu <- c(0.5,3,10)
K1 <- kernel.linear(nu=nu[1], theta=theta)
K2 <- kernel.linear(nu=nu[2], theta=theta) 
K3 <- kernel.linear(nu=nu[3], theta=theta) 

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
```

<img src="README_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

##### Reproducing Figure 3

Consider a non-linear kernel with various choices of parameter settings,
including `nu`, `gamma`, `s2`.

-   First row: Set `gamma=1` and `s2=1` and set different values for
    `nu`, which are 0.5, 2, and 10.
-   Second row: Set `nu=2.5` and `s2=1` and set different values for
    `gamma`, which are 0.1, 1, and 10.
-   Third row: Set `nu=2.5` and `gamma=1` and set different values for
    `s2`, which are 0.1, 1, and 100.

``` r
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
```

<img src="README_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

## Reproducing Section 4.2: Prediction Performance

Three different test functions are considered:

-   *f*<sub>1</sub>(*g*) = ∫∫*g*
-   *f*<sub>2</sub>(*g*) = ∫∫*g*<sup>3</sup>
-   *f*<sub>3</sub>(*g*) = ∫∫sin (*g*<sup>2</sup>)

Eight training functional inputs are

-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = *x*<sub>1</sub> + *x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = *x*<sub>1</sub><sup>2</sup>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = *x*<sub>2</sub><sup>2</sup>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>1</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>1</sub>*x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = sin (*x*<sub>1</sub>)
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = cos (*x*<sub>1</sub>+*x*<sub>2</sub>)

The domain space of *x* is \[0,1\]<sup>2</sup>.

Test functional inputs are

-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = sin (*α*<sub>1</sub>*x*<sub>1</sub>+*α*<sub>2</sub>*x*<sub>2</sub>)
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = *β* + *x*<sub>1</sub><sup>2</sup> + *x*<sub>2</sub><sup>3</sup>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = exp (−*κ**x*<sub>1</sub>*x*<sub>2</sub>)

with random *α*<sub>1</sub>, *α*<sub>2</sub>, *β* and *κ* from \[0,1\].

``` r
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
```

##### Reproducing Table 1

``` r
Y <- cbind(y1,y2,y3)
knitr::kable(round(t(Y),2))
```

|     |      |      |      |      |      |      |      |      |
|:----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| y1  | 1.00 | 0.33 | 0.33 | 1.50 | 1.50 | 1.25 | 0.46 | 0.50 |
| y2  | 1.50 | 0.14 | 0.14 | 3.75 | 3.75 | 2.15 | 0.18 | 0.26 |
| y3  | 0.62 | 0.19 | 0.19 | 0.49 | 0.49 | 0.84 | 0.26 | 0.33 |

Now we are ready to fit a FIGP model. In each for loop, we fit a FIGP
for each of `y1`, `y2` and `y3`. In each for loop, we also compute LOOCV
errors by `loocv` function.

``` r
loocv.l <- loocv.nl <- rep(0,3)
gp.fit <- gpnl.fit <- vector("list", 3)
set.seed(1)
for(i in 1:3){
  # fit FIGP with a linear kernel
  gp.fit[[i]] <- FIGP_l(G, d=2, Y[,i], nu=2.5, nug=eps, rnd=1e3)
  loocv.l[i] <- loocv(gp.fit[[i]])
  
  # fit FIGP with a nonlinear kernel
  gpnl.fit[[i]] <- FIGP_nl(G, d=2, Y[,i], nu=2.5, nug=eps, rnd=1e3)
  loocv.nl[i] <- loocv(gpnl.fit[[i]])
}
```

As a comparison, we consider a basis expansion approach using KL
expansion.

``` r
# for comparison: basis expansion approach
# KL expansion that explains 99% of the variance
set.seed(1)
KL.out <- KL.expan(d=2, G, fraction=0.99, rnd=1e3)
B <- KL.out$B
  
KL.fit <- vector("list", KL.out$M)
# fit a conventional GP on the scores
for(i in 1:KL.out$M) KL.fit[[i]] <- sepGP(B, Y[,i], nu=2.5, g=eps)
```

Let’s make predictions on the test functional inputs. We test `n.test`
times.

``` r
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
```

##### Reproducing Table 2

``` r
out <- rbind(format(loocv.l,digits=4),
             format(loocv.nl,digits=4),
             format(mse.linear,digits=4),
             format(mse.nonlinear,digits=4))
rownames(out) <- c("linear LOOCV", "nonlinear LOOCV", "linear MSE", "nonlinear MSE")
colnames(out) <- c("y1", "y2", "y3")
knitr::kable(out)
```

|                 | y1        | y2        | y3        |
|:----------------|:----------|:----------|:----------|
| linear LOOCV    | 7.327e-07 | 1.833e+00 | 4.584e-01 |
| nonlinear LOOCV | 3.764e-06 | 2.319e-01 | 1.656e-02 |
| linear MSE      | 3.226e-07 | 1.098e+02 | 1.401e+01 |
| nonlinear MSE   | 6.057e-06 | 1.206e+00 | 1.775e+00 |

##### Reproducing Table 3

``` r
select.idx <- apply(rbind(loocv.l, loocv.nl), 2, which.min)
select.mse <- diag(rbind(mse.linear, mse.nonlinear)[select.idx,])
select.cvr <- diag(rbind(cvr.linear, cvr.nonlinear)[select.idx,])
select.score <- diag(rbind(score.linear, score.nonlinear)[select.idx,])

out <- rbind(format(select.mse,digits=4),
             format(mse.kl,digits=4),
             format(select.cvr,digits=4),
             format(cvr.kl,digits=4),
             format(select.score,digits=4),
             format(score.kl,digits=4))
rownames(out) <- c("FIGP MSE", "Basis MSE", "FIGP coverage", "Basis coverage", "FIGP score", "Basis score")
colnames(out) <- c("y1", "y2", "y3")
knitr::kable(out)
```

|                | y1        | y2        | y3        |
|:---------------|:----------|:----------|:----------|
| FIGP MSE       | 3.226e-07 | 1.206e+00 | 1.775e+00 |
| Basis MSE      | 0.001739  | 8.870046  | 2.356393  |
| FIGP coverage  | 95.33     | 100.00    | 100.00    |
| Basis coverage | 75.33     | 79.00     | 49.67     |
| FIGP score     | 13.060    | 2.553     | 3.417     |
| Basis score    | 4.587     | -1.991    | -12.208   |

## Reproducing Section 5: Inverse Scattering Problems

Now we move to a real problem: inverse scattering problem. First, since
the data were generated through Matlab, we use the function `readMat` in
the package `R.matlab` to read the data. There were ten training data
points, where the functional inputs are

-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>1</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − *x*<sub>1</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>1</sub>*x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − *x*<sub>1</sub>*x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − *x*<sub>2</sub>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>1</sub><sup>2</sup>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − *x*<sub>1</sub><sup>2</sup>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + *x*<sub>2</sub><sup>2</sup>
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − *x*<sub>2</sub><sup>2</sup>

##### Reproducing Figure 4

The outputs are displayed as follows, which reproduces Figure 4.

``` r
func.title <- c("g(x1,x2)=1+x1","g(x1,x2)=1-x1","g(x1,x2)=1+x1x2","g(x1,x2)=1-x1x2",
                "g(x1,x2)=1+x2","g(x1,x2)=1-x2","g(x1,x2)=1+x1^2","g(x1,x2)=1-x1^2",
                "g(x1,x2)=1+x2^2","g(x1,x2)=1-x2^2")

output.mx <- matrix(0,nrow=10,ncol=32*32)
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
```

<img src="README_files/figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

We perform PCA (principal component analysis) for dimension reduction,
which shows that only three components can explain more than 99.9%
variation of the data.

``` r
pca.out <- prcomp(output.mx, scale = TRUE, center = TRUE)
n.comp <- which(summary(pca.out)$importance[3,] > 0.999)[1]
print(n.comp)
```

    ## PC3 
    ##   3

##### Reproducing Figure 5

Plot the three principal components, which reproduces Figure 5.

``` r
par(mfrow=c(1,3))
par(mar = c(1, 1, 2, 1))
for(i in 1:n.comp){
  eigen.vec <- matrix(c(pca.out$rotation[,i]), 32, 32)
  image(eigen.vec,yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main=paste("PC",i))
  contour(eigen.vec, add = TRUE, nlevels = 5)
}
```

<img src="README_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

Now we are ready to fit the FIGP model on those PC scores. Similarly, we
fit the FIGP with a linear kernel and a nonlinear kernel.

``` r
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
```

Perform a LOOCV to see which kernel is a better choice.

``` r
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

out <- c(loocv.linear, loocv.nonlinear)
names(out) <- c("linear", "nonlinear")
print(out)
```

    ##       linear    nonlinear 
    ## 9.617627e-08 3.115015e-05

We see linear kernel leads to a smaller LOOCV, which indicates that it’s
a better choice.

##### Reproducing Figure 6

Thus, we perform the predictions on a test input using the FIGP model
with the linear kernel, which is

-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − sin (*x*<sub>2</sub>)

``` r
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
```

<img src="README_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

The MSEs (%) for the two test data are given below.

``` r
out <- cbind(mse, mse.kl)
colnames(out) <- c("FIGP", "Basis")
print(out)
```

    ##              FIGP        Basis
    ## [1,] 0.0001309401 0.0006983163
