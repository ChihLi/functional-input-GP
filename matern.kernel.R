matern.kernel <- function(r, nu) 
{
  if(nu==1.5){
    out <- (1+r*sqrt(3)) * exp(-r*sqrt(3))
  }else if(nu==2.5){
    out <- (1+r*sqrt(5)+5*r^2/3) * exp(-r*sqrt(5))
  }else{
    rat <- r*sqrt(2*nu)
    out <- (2^(1 - nu))/gamma(nu) * rat^nu * besselK(rat, nu)
    out[is.nan(out)] <- 1
  }
  return(out)
}
