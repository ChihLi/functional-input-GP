loocv <- function(fit){
  
  Ki <- fit$Ki
  theta <- fit$theta
  y <- fit$y
  mu.hat <- fit$mu.hat
  
  return(mean((diag(1/diag(Ki)) %*% Ki %*% (y-mu.hat))^2))
}

loocv.pred <- function(fit){
  
  Ki <- fit$Ki
  theta <- fit$theta
  y <- fit$y
  mu.hat <- fit$mu.hat
  
  out <- rep(0, length(y))
  for(i in 1:length(y)) out[i] <- mu.hat -  Ki[i,-i] %*% (y-mu.hat)[-i] / diag(Ki)[i]
  return(out)
}