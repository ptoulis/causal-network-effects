##   Responses.  
##  Creates the linear response model
## Yi = ZNi * γ + Zi * τ + ε  
linear.y = function(direct.tau=10, 
                    indirect.gamma=2.5, 
                    error.sigma=0.25) {
  fn = function(x,v) {
    zn = get.ZNi(x,v)
    if(sum(is.na(zn))>0) {
      stop("Zn should not have NA elements. Some neighbors are not assigned.")
    }
    direct.effect = get.Zi(x, v) * direct.tau
    indirect.effect = sum(zn * indirect.gamma)
    error = rnorm(1, sd=error.sigma)
    
    direct.effect + indirect.effect + error
  }
}