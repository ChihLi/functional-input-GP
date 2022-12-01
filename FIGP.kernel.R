FIGP.kernel <- function(d, theta, nu, G, g=NULL, kernel, ...){
  if(kernel=="linear"){
    # hcubature could take a long time to converge (for some reasons)
    return(FIGP.kernel.linear(d, theta, nu, G, g, MC.fg=TRUE, rnd=5000, ...)) 
  }else{
    # hcubature is stable for nonlinear cases
    return(FIGP.kernel.nonlinear(d, theta, nu, G, g, MC.fg=TRUE, rnd=5000, ...))
  }
}

FIGP.kernel.linear <- function(d, theta, nu, G, g=NULL, MC.fg=FALSE, rnd=1e3, ...){
  
  n <- length(G)
  if(!is.null(g)) n.new <- length(g)
  
  if(is.null(g)){
    if(!MC.fg){ # numerical integration by the R function hcubature
      K <- matrix(0, n, n)
      for(i in 1:n){
        for(j in 1:i){
          int.fun <- function(x) G[[i]](x[1:d])*G[[j]](x[(d+1):(2*d)])*matern.kernel(sqrt(sum(((x[1:d]-x[(d+1):(2*d)])/theta)^2)), nu=nu)
          K[i,j] <- K[j,i] <- 
            hcubature(int.fun, rep(0,2*d),rep(1,2*d), ...)$integral
        }
      }
    }else{ # numerical integration by Monte-Carlo approximation
      X <- randtoolbox::sobol(rnd, d)
      R <- sqrt(distance(t(t(X)/theta)))
      Phi <- matern.kernel(R, nu=nu)
      A <- matrix(0,ncol=n,nrow=rnd)
      for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
      K <- (t(A) %*% Phi %*% A) / rnd # should be /rnd^2 but the values become too small, but it doesn't hurt without it because of scale parameter
      K <- (K+t(K))/2
    }
  }else{
    if(!MC.fg){ # numerical integration by the function integrate
      K <- matrix(0, nrow=n.new, ncol=n)
      for(i in 1:n.new){
        for(j in 1:n){
          int.fun <- function(x) g[[i]](x[1:d])*G[[j]](x[(d+1):(2*d)])*matern.kernel(sqrt(sum(((x[1:d]-x[(d+1):(2*d)])/theta)^2)), nu=nu)
          K[i,j] <- hcubature(int.fun, rep(0,2*d),rep(1,2*d), ...)$integral
        }
      }
    }else{ # numerical integration by Monte-Carlo approximation
      X <- randtoolbox::sobol(rnd, d)
      R <- sqrt(distance(t(t(X)/theta)))
      Phi <- matern.kernel(R, nu=nu)
      A <- matrix(0,ncol=n,nrow=rnd)
      a <- matrix(0,ncol=n.new,nrow=rnd)
      for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
      for(i in 1:n.new)  a[,i] <- apply(X, 1, g[[i]])
      K <- (t(a) %*% Phi %*% A) / rnd # should be /rnd^2 but the values become too small, but it doesn't hurt without it because of scale parameter
    }
  }
  
  return(K)
}



FIGP.kernel.nonlinear <- function(d, theta, nu, G, g=NULL, MC.fg=FALSE, rnd=1e3, ...){
  
  n <- length(G)
  if(!is.null(g)) n.new <- length(g)
  
  if(is.null(g)){
    if(!MC.fg){ # numerical integration by the R function hcubature
      K <- matrix(0, n, n)
      for(i in 1:n){
        for(j in 1:i){
          diffsq <- function(x) (G[[i]](x)-G[[j]](x))^2
          l2.sq <- hcubature(diffsq, lower=rep(0,d),upper=rep(1,d), ...)$integral
          K[i,j] <- K[j,i] <- matern.kernel(sqrt(l2.sq)/theta, nu=nu)
        }
      }
    }else{ # numerical integration by Monte-Carlo approximation
      X <- randtoolbox::sobol(rnd, d)
      A <- matrix(0, ncol=n, nrow=rnd)
      for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
      
      R <- sqrt(distance(t(A))/rnd)
      K <- matern.kernel(R/theta, nu=nu)
    }
  }else{
    if(!MC.fg){ # numerical integration by the function hcubature
      K <- matrix(0, nrow=n.new, ncol=n)
      for(i in 1:n.new){
        for(j in 1:n){
          diffsq <- function(x) (g[[i]](x)-G[[j]](x))^2
          l2.sq <- hcubature(diffsq, lower=rep(0,d),upper=rep(1,d), ...)$integral
          K[i,j] <- matern.kernel(sqrt(l2.sq)/theta, nu=nu)
        }
      }
    }else{ # numerical integration by Monte-Carlo approximation
      X <- randtoolbox::sobol(rnd, d)
      a <- matrix(0,ncol=n.new,nrow=rnd)
      for(i in 1:n.new)  a[,i] <- apply(X, 1, g[[i]])
      A <- matrix(0, ncol=n, nrow=rnd)
      for(i in 1:n)  A[,i] <- apply(X, 1, G[[i]])
      
      R <- sqrt(distance(t(a),t(A))/rnd)
      K <- matern.kernel(R/theta, nu=nu)
    }
  }
  
  return(K)
}


# library(lhs)
# function.dist.nl <- function(d, G, g=NULL, norm=c("L2", "uniform"), init.samples=15){
#   
#   n <- length(G)
#   if(!is.null(g)) n.new <- length(g)
#   
#   if(norm=="L2"){
# 
#   }else if(norm=="uniform"){
#     init.mx <- maximinLHS(init.samples, d)
#     
#     if(is.null(g)){
#       R <- matrix(0, n, n)
#       for(i in 1:n){
#         for(j in 1:i){
#           negdiff <- function(x) -abs(G[[i]](x)-G[[j]](x))
#           R[i,j] <- R[j,i] <- max(apply(init.mx, 1, function(x){
#             -optim(x,negdiff,lower=0,upper=1,method="L-BFGS-B")$value
#           }))
#         }
#       }
#     }else{
#       R <- matrix(0, nrow=n.new, ncol=n)
#       for(i in 1:n.new){
#         for(j in 1:n){
#           negdiff <- function(x) -abs(g[[i]](x)-G[[j]](x))
#           R[i,j] <- max(apply(init.mx, 1, function(x){
#             -optim(x,negdiff,lower=0,upper=1,method="L-BFGS-B")$value
#           }))
#         }
#       }
#     }
#   }
#   
#   return(R)
# }