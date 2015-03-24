###  Interference

fi  = function(k, n) {
  error = (k>n);
  if(error)
    stop(sprintf("Error occured k=%d n=%d", k, n))
  if(k==0) return(1.0)
  if(k==1 && n==2) return(0.5)
  if(k==n) return(1/n)
  
  pK = k/n
  pN = (1-k/n)
  tmp = rbinom(1, size=n-k, p=pN)
  X = rbinom(1, size=k, p=pK)
  Y = tmp+X
  if(is.na(X) || is.na(Y) ) {
    print(sprintf("pK = %.1f pN=%.1f k=%d  n=%d  X=%d,   Y=%d",
                  pK, pN, k,n, X,Y))
    stop()
  }
  
  a = (n-k-1) / n;
  return( a* fi(k-X, n-Y) + 1/n)
}