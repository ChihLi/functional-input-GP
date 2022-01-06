This instruction aims to reproduce the results in the paper
“*Functional-Input Gaussian Processes with Applications to Inverse
Scattering Problems*” by Sung et al. Hereafter, functional-Input
Gaussian Process is abbreviated by *FIGP*.

The following results are reproduced in this file

-   The sample path plots in Section 4.1 (Figures 2 and 3)
-   The prediction results in Section 4.2 (Tables 1 and 2)
-   The plots and prediction results in Section 5 (Figures 4, 5, and 6)

##### Step 0.1: load functions and packages

``` r
library(randtoolbox)
library(cubature)
library(plgp)
library(R.matlab)
source("sepGP.R")               # FIGP with linear kernel
source("sepGP_nonlinear.R")     # FIGP with non-linear kernel
source("matern.kernel.R")       # matern kernel computation
source("loocv.R")               # LOOCV for FIGP
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
  R <- sqrt(distance(t(A)*theta))
  
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
K2 <- kernel.linear(nu=nu[1], theta=theta) 
K3 <- kernel.linear(nu=nu[3], theta=theta) 

par(mfrow=c(3,3), mar = c(4, 4, 2, 1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(nu==1/2))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(nu==3))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, 
        xlab=expression(alpha), ylab="y", main=expression(nu==10))

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
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, 
        xlab=expression(alpha), ylab="y", main=expression(theta==100))

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
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[3]*K3)), type="l", col=3, lty=3, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==100))
```

<img src="README_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

##### Reproducing Figure 3

Consider a non-linear kernel with various choices of parameter settings,
including `nu`, `gamma`, `s2`.

-   First row: Set `gamma=0.01` and `s2=1` and set different values for
    `nu`, which are 0.5, 2, and 10.
-   Second row: Set `nu=2.5` and `s2=1` and set different values for
    `gamma`, which are 0.001, 0.01, and 0.1.
-   Third row: Set `nu=2.5` and `gamma=0.01` and set different values
    for `s2`, which are 0.01, 1, and 100.

``` r
gamma <- 0.01
s2 <- 1
nu <- c(0.5,2,10)
K1 <- kernel.nonlinear(nu=nu[1], theta=gamma)
K2 <- kernel.nonlinear(nu=nu[1], theta=gamma) 
K3 <- kernel.nonlinear(nu=nu[3], theta=gamma) 

par(mfrow=c(3,3), mar = c(4, 4, 2, 1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(nu==1/2))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(nu==2))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, 
        xlab=expression(alpha), ylab="y", main=expression(nu==10))

nu <- 2.5
gamma <- c(0.001,0.01,0.1)
s2 <- 1
K1 <- kernel.nonlinear(nu=nu, theta=gamma[1])
K2 <- kernel.nonlinear(nu=nu, theta=gamma[2]) 
K3 <- kernel.nonlinear(nu=nu, theta=gamma[3])
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(gamma==0.001))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(gamma==0.01))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2*K3)), type="l", col=3, lty=3, 
        xlab=expression(alpha), ylab="y", main=expression(gamma==0.1))

nu <- 2.5
gamma <- 0.01
s2 <- c(0.1,1,100)
K1 <- kernel.nonlinear(nu=nu, theta=gamma)
K2 <- kernel.nonlinear(nu=nu, theta=gamma) 
K3 <- kernel.nonlinear(nu=nu, theta=gamma) 

matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[1]*K1)), type="l", col=1, lty=1, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==0.1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[2]*K2)), type="l", col=2, lty=2, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==1))
matplot(seq(0,1,0.01), t(rmvnorm(8,sigma=s2[3]*K3)), type="l", col=3, lty=3, 
        xlab=expression(alpha), ylab="y", main=expression(sigma^2==100))
```

<img src="README_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

## Reproducing Section 4.2: Prediction Performance

Three different test functions are considered:

-   *f*<sub>1</sub>(*g*) = ∫∫*g*
-   *f*<sub>2</sub>(*g*) = ∫∫*g*<sup>2</sup>
-   *f*<sub>3</sub>(*g*) = ∫∫sin (*g*)

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

-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 + sin (*α*<sub>1</sub>*x*<sub>1</sub>+*α*<sub>2</sub>*x*<sub>2</sub>)
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

# y2: integrate g^2 function from 0 to 1
G.square <- list(function(x) (x[1]+x[2])^2,
           function(x) (x[1]^2)^2,
           function(x) (x[2]^2)^2,
           function(x) (1+x[1])^2,
           function(x) (1+x[2])^2,
           function(x) (1+x[1]*x[2])^2,
           function(x) (sin(x[1]))^2,
           function(x) (cos(x[1]+x[2]))^2)
y2 <- rep(0, n) 
for(i in 1:n) y2[i] <- hcubature(G.square[[i]], lower=c(0, 0),upper=c(1,1))$integral

# y3: integrate sin(g) function from 0 to 1
G.sin <- list(function(x) sin((x[1]+x[2])),
              function(x) sin((x[1]^2)),
              function(x) sin((x[2]^2)),
              function(x) sin((1+x[1])),
              function(x) sin((1+x[2])),
              function(x) sin((1+x[1]*x[2])),
              function(x) sin((sin(x[1]))),
              function(x) sin((cos(x[1]+x[2]))))
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
| y2  | 1.17 | 0.20 | 0.20 | 2.33 | 2.33 | 1.61 | 0.27 | 0.35 |
| y3  | 0.77 | 0.31 | 0.31 | 0.96 | 0.96 | 0.93 | 0.43 | 0.45 |

Now we are ready to fit a FIGP model. In each for loop, we fit a FIGP
for each of `y1`, `y2` and `y3`. In each for loop, we also compute LOOCV
errors by `loocv` function.

``` r
loocv.l <- loocv.nl <- rep(0,3)
gp.fit <- gpnl.fit <- vector("list", 3)
set.seed(1)
for(i in 1:3){
  # fit FIGP with a linear kernel
  gp.fit[[i]] <- sepGP(G, d=2, Y[,i], nu=2.5, nug=eps, rnd=1e3)
  loocv.l[i] <- loocv(gp.fit[[i]])
  
  # fit FIGP with a nonlinear kernel
  gpnl.fit[[i]] <- sepGP_nl(G, d=2, Y[,i], nu=2.5, nug=eps, rnd=1e3)
  loocv.nl[i] <- loocv(gpnl.fit[[i]])
}
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

mape.linear <- mape.nonlinear <- rep(0,3)
for(i in 1:3){
  mape.linear.i <- mape.nonlinear.i <- rep(0, n.test)
  for(ii in 1:n.test){
    gnew <- list(function(x) 1+sin(alpha1[ii]*x[1]+alpha2[ii]*x[2]),
                 function(x) beta1[ii]+x[1]^2+x[2]^3,
                 function(x) exp(-kappa1[ii]*x[1]*x[2]))    
    if(i==1){
      g.int <- gnew
    }else if(i==2){
      g.int <- list(function(x) (1+sin(alpha1[ii]*x[1]+alpha2[ii]*x[2]))^2,
                    function(x) (beta1[ii]+x[1]^2+x[2]^3)^2,
                    function(x) (exp(-kappa1[ii]*x[1]*x[2]))^2)
    }else if(i==3){
      g.int <- list(function(x) sin(1+sin(alpha1[ii]*x[1]+alpha2[ii]*x[2])),
                    function(x) sin(beta1[ii]+x[1]^2+x[2]^3),
                    function(x) sin(exp(-kappa1[ii]*x[1]*x[2])))
    }
    
    n.new <- length(gnew)
    y.true <- rep(0,n.new)
    for(iii in 1:n.new) y.true[iii] <- hcubature(g.int[[iii]], lower=c(0, 0),upper=c(1,1))$integral
    
    ynew <- pred.sepGP(gp.fit[[i]], gnew)
    mape.linear.i[ii] <- mean(abs((y.true - ynew$mu)/y.true))
    
    ynew <- pred.sepGP_nl(gpnl.fit[[i]], gnew)
    mape.nonlinear.i[ii] <- mean(abs((y.true - ynew$mu)/y.true))
  }
  mape.linear[i] <- mean(mape.linear.i)*100
  mape.nonlinear[i] <- mean(mape.nonlinear.i)*100
}
```

##### Reproducing Table 2

``` r
out <- rbind(format(loocv.l,digits=4),
             format(loocv.nl,digits=4),
             format(mape.linear,digits=4),
             format(mape.nonlinear,digits=4))
rownames(out) <- c("linear LOOCV", "nonlinear LOOCV", "linear MAPE", "nonlinear MAPE")
colnames(out) <- c("y1", "y2", "y3")
knitr::kable(out)
```

|                 | y1        | y2        | y3        |
|:----------------|:----------|:----------|:----------|
| linear LOOCV    | 1.486e-09 | 1.751e-01 | 3.084e-02 |
| nonlinear LOOCV | 3.764e-06 | 2.095e-02 | 1.390e-03 |
| linear MAPE     | 6.744e-04 | 3.688e+01 | 1.630e+01 |
| nonlinear MAPE  | 0.01599   | 5.69984   | 3.20862   |

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

<img src="README_files/figure-markdown_github/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

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

<img src="README_files/figure-markdown_github/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

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
  gp.fit[[i]] <- sepGP(G, d=2, y, nu=2.5, nug=eps, rnd=1e3)
  # fit FIGP with a nonlinear kernel    
  gpnl.fit[[i]] <- sepGP_nl(G, d=2, y, nu=2.5, nug=eps, rnd=1e3)
}
```

Perform a LOOCV to see which kernel is a better choice.

``` r
loocv.recon <- cbind(loocv.pred(gp.fit[[1]]), 
                     loocv.pred(gp.fit[[2]]), 
                     loocv.pred(gp.fit[[3]])) %*% t(pca.out$rotation[,1:n.comp])
loocv.recon <- scale(loocv.recon, center = FALSE, scale = 1/pca.out$scale)
loocv.recon <- scale(loocv.recon, scale = FALSE, center = -1 * pca.out$center)
loocv.linear <- mean((loocv.recon - output.mx)^2)

loocv.nl.recon <- cbind(loocv.pred(gpnl.fit[[1]]), 
                        loocv.pred(gpnl.fit[[2]]), 
                        loocv.pred(gpnl.fit[[3]])) %*% t(pca.out$rotation[,1:n.comp])
loocv.nl.recon <- scale(loocv.nl.recon, center = FALSE, scale = 1/pca.out$scale)
loocv.nl.recon <- scale(loocv.nl.recon, scale = FALSE, center = -1 * pca.out$center)
loocv.nonlinear <- mean((loocv.nl.recon - output.mx)^2)

out <- c(loocv.linear, loocv.nonlinear)
names(out) <- c("linear", "nonlinear")
print(out)
```

    ##       linear    nonlinear 
    ## 9.617627e-08 9.898094e-08

We see linear kernel leads to a smaller LOOCV, which indicates that it’s
a better choice.

##### Reproducing Figure 6

Thus, we perform the predictions on the test inputs using the FIGP model
with a linear kernel, which are

-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1 − sin (*x*<sub>2</sub>)
-   *g*(*x*<sub>1</sub>,*x*<sub>2</sub>) = 1

``` r
# test functional inputs (gnew)
gnew <- list(function(x) 1-sin(x[2]),
             function(x) 1)
n.new <- length(gnew)

# make predictions using a linear kernel
ynew <-  matrix(0,ncol=n.comp,nrow=n.new)
for(i in 1:n.comp){
  y <- pca.out$x[,i]
  # fit FIGP with a linear kernel  
  gp.fit[[i]] <- sepGP(G, d=2, y, nu=2.5, nug=eps, rnd=1e3)
  pred.out <- pred.sepGP(gp.fit[[i]], gnew)
  ynew[,i] <- pred.out$mu
}

# reconstruct the image
pred.recon <- ynew %*% t(pca.out$rotation[,1:n.comp])
pred.recon <- scale(pred.recon, center = FALSE, scale = 1/pca.out$scale)
pred.recon <- scale(pred.recon, scale = FALSE, center = -1 * pca.out$center)

# true data on the test data
gnew.true <- matrix(0, ncol=n.new, nrow=32*32)
gnew.dat <- readMat(paste0("DATA/q_sine.mat"))$Ffem
gnew.true[,1] <- c(Re(gnew.dat))
gnew.dat <- readMat(paste0("DATA/q_one.mat"))$Ffem
gnew.true[,2] <- c(Re(gnew.dat))

# plot the result
par(mfrow=c(1,4))
par(mar = c(1, 1, 2, 1))

mape <- rep(0, n.new)
for(i in 1:n.new){
  image(matrix(gnew.true[,i],32,32), zlim=c(0.064,0.11),yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main=ifelse(i==1, "g(x1,x2)=1-sin(x2)", "g(x1,x2)=1"))
  contour(matrix(gnew.true[,i],32,32), add = TRUE, nlevels = 5)
  
  image(matrix(pred.recon[i,], 32, 32), zlim=c(0.064,0.11),yaxt="n",xaxt="n",
        col=heat.colors(12, rev = FALSE),
        main="prediction")
  contour(matrix(pred.recon[i,], 32, 32), add = TRUE, nlevels = 5)
  mape[i] <- mean(abs((gnew.true[,i]-pred.recon[i,])/gnew.true[,i]))*100
}
```

<img src="README_files/figure-markdown_github/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

The MAPEs (%) for the two test data are given below.

``` r
print(mape)
```

    ## [1] 1.10113370 0.03438794
