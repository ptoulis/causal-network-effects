###  Synopsis.  Example use of code.
source("cii.R")
source("graphs.R")
source("estimands.R")
source("responses.R")
source("randomizations.R")
source("var-estim.R")
source("../../../code/R/libs/CausalInference.R")
## Load the saved object.
## This will load the    cii   object.
print(">   Loading cii object")
load(file="data/barabasi20.Rdata")
## This will load the   Y    potential outcomes.
print(">    Loading potential outcomes")
load(file="data/PotentialOutcomes.Rdata")
#  This will load Y.simple  = simple peer effects.
print(">    Loading simple potential outcomes")
load(file="data/PotentialOutcomes-simplify.Rdata")
print(">    Loading W and Q matrices")
load(file="data/WQ.Rdata")

plot.cii(cii)
#print(Y)

##  ALL potential outcomes for the units.
##  If we SIMPLIFY then all active treatments are the same.
generate.potential.outcomes = function(cii, mu.t=5, mu.c=0,
                                       sigma.t= 2.0,  sigma.c=0.3,
                                       simplify=F) {
  Vk = get.Vk(cii)  ##  has the ids of eligible nodes.
  k = cii$k
  Y = list()
  for(i in Vk) {
    ni = length(get.neighbors(cii, i))
    Nik = choose(ni, k) ##  the total no. of treatment arms
    Y[[i]] = rep(NA, Nik+1)
    if(simplify) {
      a =rnorm(1, mean=mu.t, sd=sigma.t)
      Y[[i]][1:Nik] = rep(a, Nik)  ## repeat the outcomes
    } else {
      for(j in 1:Nik) 
        Y[[i]][j] = rnorm(1, mean=mu.t, sd=sigma.t)
    }
    ##  Set the control response
    Y[[i]][Nik+1] = rnorm(1, mean=mu.c, sd=sigma.c)
  }
  ## Returns all potential outcomes as a big list.
  return(Y)
}
##  Given a vector  z= (0,0,0,1,1,0....)   return an index from 1 ... choose(ni,k)+1
## max is control.
index.ZNi = function(cii, i, z) {
  ni = length(get.neighbors(cii, i))
  nik = choose(ni, cii$k)
  k = cii$k
  if(sum(na.omit(z))==0)  {
   return(nik+1)
  } else{
    M = assign.complete.rand.matrix(ni,  k)
    for(j in 1:nrow(M)) 
      if(sum(M[j,] == z)==ni)
      {
       return(j)
      }
  }
  stop("Error. Could convert to Wi")
}




get.Y.outcome = function(i, treatment.index, simplify) {
  if(simplify) return(Y.simple[[i]][treatment.index])
  return(Y[[i]][treatment.index])
}
##  Given a cii do the following
##  1. Draw a sequential design sample
##  2. Compute the delta estimate
##  3.  Return   W = m x N = trials x units    where Wij = 1 if unit j is in trial i
##               Q = m x N = trials x units    where Qij = 1 iff unit j is controlled in trial i
many.estimates = function(cii, nestims, simplify) {
  
  pb = txtProgressBar(style=3)
  Vk = get.Vk(cii)
  W = matrix(NA, nrow=0, ncol=length(Vk))
  Q = matrix(NA, nrow=0, ncol=length(Vk))
  delta = c()
   
  for(n in 1:nestims) {
    x = design.simple.sequential(cii)
    ## Total number of treated/controls.
    w = sapply(Vk, function(v) as.numeric(v %in% x$ex$Exposed))  
    q = sapply(Vk, function(v) as.numeric(v %in% x$ex$Control))
    Nt= sum(w)
    Nc= sum(q)
    if(Nt * Nc ==0) next;
    yt = c()
    yc = c()
    # 1.  Get the treatment vector   Yt
    for(i in x$ex$Exposed)  {
      zni = get.ZNi(x, i)
      yti =  get.Y.outcome(i, index.ZNi(cii, i, z=zni), simplify)
      yt = c(yt, yti)
    }
    ## Good idea to delete iterators. Avoids bugs.
    rm(i)
    ## 2. Get the control vector   Yc
    for(j in x$ex$Control) {
      zni = get.ZNi(x, j)
      yci =  get.Y.outcome(j, index.ZNi(cii, j, z=zni), simplify)
      yc = c(yc, yci)
    }
    ## Remove this index.
    rm(j)
    
    
    val = mean(yt) - mean(yc)
    delta = c(delta, val)
    setTxtProgressBar(pb, value=n/nestims)
    
      
    W = rbind(W, w)
    Q = rbind(Q, q)
  }
  
  return(list(W=W, Q=Q,delta=delta))
  #return(mt-mc)  
}

###    Analysis of variance
##     We need to find out the sources of variance for our peer effects estimator
analysis.of.variance = function(cii, trials=100) {
  z = many.estimates(cii, nestims=trials, simplify=T)
  Nt = rowSums(z$W)
  Nc= rowSums(z$Q)

  vals = z$delta
  
 
  ##  Get Yt, Yc 
  Vk=  get.Vk(cii)
  yt = matrix(sapply(Vk, function(i) head(Y.simple[[i]], 1)), nrow=1)
  yc = matrix(sapply(Vk, function(i) tail(Y.simple[[i]],1)), nrow=1)
  
  
  Vw = yt %*% cov(z$W/Nt) %*% t(yt)
  Vq = yc %*% cov(z$Q/Nc) %*% t(yc)
  Vwq = 2 * yt %*% cov(z$W/Nt, z$Q/Nc) %*% t(yc)
  
  theor.var = Vw + Vq - Vwq
  
  print(sprintf("Variance %.3f", var(vals)))
  print(sprintf("Theoretical %.3f", theor.var))
  
}



neyman = function() {
  N = 100
  Nt = 40
  Nc = N-Nt
  Yt = rnorm(100, mean=4, sd=8)
  Yc = rnorm(100, mean=0, sd=1)
  U = rep(1, N) %*% t(rep(1,N))
  A = diag(N) - (1/N) * U
  
  theor.vars = function() {
    vt = (Nc/N) * var(Yt)/Nt
    vc = (Nt/N) * var(Yc)/Nc
    print("Theoretical variances")
    print(sprintf("VarT = %.3f  VarC=%.3f", vt, vc))
    cov.theor = -t(Yt) %*% A %*% Yc / (N*(N-1))
    print(sprintf("Covariance = %.3f", cov.theor))
    
  }
  sample.obs = function() {
    W = sample( c(rep(1, Nt), rep(0, Nc)))
    W0 = 1-W
    c( (1/Nt) * sum(W * Yt), (1/Nc) * sum(W0 * Yc) )
  }
  
  trials =10^5
  reps = replicate(trials, { sample.obs()})
  print("Empirical variances")
  print(apply(reps, 1, var))
  
  print("Empirical covariance")
  se = sd(replicate(1000, { i =sample(1:trials, replace=T); cov(reps[1,i], reps[2,i])}))
  cov.pt = cov(reps[1,], reps[2,])
  print(sprintf("Confidence is  %.5f, %.5f", cov.pt-2*se, cov.pt+2*se))
  
  theor.vars()
}


