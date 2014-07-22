##  Goal is to estimate ANY quadratic form.
rm(list=ls())
N = 10
Nt = N/2
Nc = N-Nt
true.var = 2.38^2
trials = 10


##  Draw a graph-like constraint between variables.
G = matrix(rbinom(N^2, size=1, prob=0.2), nrow=N)
G  = ceiling((G+t(G))/2)
diag(G)=1
get.neighbors= function(i)  which(G[i,]==1); 
Nei = sapply(1:N, function(i) length(get.neighbors(i)))

##  draw outcomes
print("1.  Sampling outcomes.")
yt = rnorm(N, mean=15, sd=sqrt(true.var))
yc = rnorm(N, mean=5, sd=1)
print("2.  Sampling a quadratic form.")
V = matrix(runif(N^2, min=0, max=1), nrow=N)


# Draw the assignment. Returns W, Q, vectors of units in treatment or control, or NA otherwise.
## Be careful of the network connections. 
draw.W.Q = function() {
  all.pairs = 1:nrow(G)
  candidate.set = 1:nrow(G)
  sample.treat.set = c()
  sample.control.set = c()
  while(length(candidate.set)>0) {
    s1 = ifelse(length(candidate.set)>1,
                sample(candidate.set, size=1, replace=F),
                candidate.set[1])
    if(runif(1)<0.5) {
      sample.treat.set = c(sample.treat.set, s1)
    } else {
      sample.control.set = c(sample.control.set, s1)
    }
    neighbors = get.neighbors(s1)
    #print(sprintf("Knocking out neighbors %d for unit %d", neighbors, s1))
    
    
    ## knock out the neighbors
    candidate.set = setdiff(candidate.set, 
                            intersect(candidate.set, c(s1, neighbors) ))
  }
  x = list(w= as.numeric(all.pairs %in% sample.treat.set),
              q = as.numeric(all.pairs %in% sample.control.set),
              na= as.numeric(all.pairs %in% setdiff(all.pairs, c(sample.treat.set, sample.control.set))))
  if(sum(x$q * x$w)>0)
    stop("SEVERE error. Should not have W=1 and Q=1 at the same time.")
  return(x)
}

##  The projection matrix
U = function(n) rep(1,n) %*% t(rep(1, n))
proj = function() (diag(N)-(1/N) * U(N))

## Compute the joint inclusion probabilities
print("4.  Computing joint inclusion probabilities.")
trials = 10000
W = matrix(NA, nrow=trials, ncol=N)
Q = matrix(NA, nrow=trials, ncol=N)
pb = txtProgressBar(style=3,)
for(i in 1:trials) {
  a = draw.W.Q()
  W[i,] = a$w
  Q[i,] = a$q
  setTxtProgressBar(pb, value=i/trials)
}

Pi = colMeans(W)
Qi = colMeans(Q)
Tij = (t(W) %*% Q)/trials
## to avoid diving by zero
#  the actual value does not matter since Vij=0 when Pij=0
Tij2 = Tij
Tij2[Tij2==0] = 10^5  

bound.f = function(x) x^2/2
bound.g = function(x) 0

############     ESTIMATION PART   #######################

###   Find estimators for the quadratic form  x' V y
##  Draw the assignment
quad.est= function(x, V, y) {
  ##  Compute the sample joint probabilities.
  ## split V into obs   and mis
  V.obs = V * (1-G)
  V.mis = V * G  ## this part will be missing.
  V.mis.plus = V.mis * (V.mis>0)
  V.mis.minus = V.mis * (V.mis<0)
  ## Goal:  Estimate  yt' V yc     quadratic form
  ## In general, V is positive definite, but it is general for now.
  #print("3. Computing estimands.")
  v.obs = t(x) %*% V.obs %*% y 
  v.mis.plus = t(x) %*% V.mis.plus %*% y
  v.mis.minus = t(x) %*% V.mis.minus %*% y
  #print(sprintf("True estimand %.3f  obs=%.3f  mis=%.3f", est, est1, est.mis))
  ## Draw the assignment
  a = draw.W.Q()
  W = diag(a$w)
  Q = diag(a$q)
  

  ##   B matrix
  ## B =  (Bij)   such that Bij = Vij / Pij   or 0 if Pij=0
  B = V.obs / Tij2
  
  v.obs.hat = t(W %*% x) %*% B %*% (Q %*% y)
  
  r1.f = a$w * bound.f(x) / Pi
  r2.f = a$q * bound.f(y) / Qi

  r1.g = a$w * bound.g(x) / Pi
  r2.g = a$q * bound.g(y) / Qi
  
  v.mis.plus.hat= t(r1.f) %*% colSums(V.mis.plus) +  t(r2.f) %*% rowSums(V.mis.plus)
  v.mis.minus.hat = t(r1.g) %*% colSums(V.mis.minus) +  t(r2.g) %*% rowSums(V.mis.minus)
  

  return(list(true.all=v.obs+v.mis.plus+v.mis.minus, 
              est.all = v.obs.hat+v.mis.plus.hat+v.mis.minus.hat, 
              v.obs=v.obs, v.obs.hat = v.obs.hat,
              v.mis.plus = v.mis.plus, v.mis.plus.hat = v.mis.plus.hat,
             v.mis.minus= v.mis.minus, v.mis.minus.hat = v.mis.minus.hat))
  
}

##  Causal estimand
causal.est = (1/N) * (t(yt) %*% Pi - t(yc) %*% Qi)

causal.estimator.ci = function() {
  a = draw.W.Q()
  causal.est.hat = (1/N) * (  t(yt) %*% a$w - t(yc) %*% a$q )

  var.true = (N^-2) * ( t(yt) %*% cov(W) %*% yt + t(yc) %*% cov(Q) %*% yc
          -2 * t(yt) %*% cov(W, Q) %*% yc )
  var.est = (N^-2) * ( quad.est(yt, cov(W), yt)$est.all + quad.est(yc, cov(Q), yc)$est.all
                       -2 * t(yt) %*% cov(W, Q) %*% yc )
  return( c(causal.est.hat-2* sqrt(var.est), causal.est.hat+2*sqrt(var.est)))
}


