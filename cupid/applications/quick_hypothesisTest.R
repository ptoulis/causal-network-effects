# Copyright (c) 2015
# Panos Toulis, ptoulis@fas.harvard.edu
# 
#
source("cupid_honduras.R")
source("wrapper_honduras.R")

#
# Code to perform hypothesis testing for causal effects 
# with partially-revealed interference.
#
# Example run.
#   s = sanitize.village.data(village.data(2))
#   with(s, test.H0(Y, E, Z, X, G, nreps=10))
#
#
plot.state <- function(saneData) {
  # s = sanitize.village.data(village.data(village.id = 2))
  # plot.villageData(s ,plot.coupon.only =T)
}

impute.E <- function(Yobs, Eobs, Z, X, G, verbose=F) {
  # Impute the missing edges.
  #
  # Fit a simple model, then impute. The model is like this 
  # p(E) = prod_ij Pr(e_ij=1)
  #  P(e_ij=1) = logistic(b0 + xi'B1 x_j + b2 * g_ij)
  N = length(Yobs)
   
  get.probs <- function(b) {
    b0 = head(b, 1)
    b1 = matrix(b[2:10], nrow=3)
    B1 = b1
    b2 = tail(b, 1)
    # X = 2xn, G = nxn
    
    W = b0 * matrix(1, nrow=N, ncol=N) + 
      X %*% B1 %*% t(X) + 
      b2 * G
    
    Z = matrix(rnorm(N**2, sd=1e-2), nrow=N)
    W[W < -5] <- -5 + Z[W < -5]
    W[W > 5] <- 5 + Z[W > 5]
    P = exp(W) / (1+exp(W))
    return(P)
  }

  log.lik <- function(b) {
    P = get.probs(b)
    # print(P)
    Q = 1-P
    x = sum(log(P[which(Eobs==1)])) + sum(log(Q[which(Eobs==0)]))
    if(x < 0 && is.infinite(x)) {
      return(-1e5) 
    } else if (x > 0 && is.infinite(x)) {
      return (1e5)
    }
    return(x)
  }
  beta = c(0, rep(0, 9), 0)
  
  # Fit a parametric model for the observed connections.
  fit = optim(par=beta, method = "L-BFGS", fn=
                function(b) {
                  ret = -log.lik(b)
                  if(verbose) print(ret)
                  return(ret)
                })
  
  P = get.probs(fit$par)
  # It holds: matrix(as.vector(X), nrow=nrow(X)) == X
  # Sample a new matrix of connectins.
  Eimp = rbinom(N**2, size=1, prob=as.vector(P))
  Eimp = matrix(Eimp, nrow=N)
  Eimp <- Eimp + t(Eobs) # we assume symmetry so Eobs^T is known.
  Eimp <- Eimp + t(Eimp) # force symmetry.
  Eimp[Eimp > 1] <- 1 
  diag(Eimp) <- 0  # no self-loops
  return(Eimp)
}

test.imputation <- function() {
  s = sanitize.village.data(village.data(2))
  Eimp = with(s, impute.E(Y, E, Z, X, G))
  par(mfrow=c(2, 2))
  # observed connections
  plot.villageData(s, T)
  
  s2 = s
  s2$E = Eimp - s$E - t(s$E)
  # imputed connections - obs  - T(obs)
  # those are the connections that were added from the model of impute.E
  plot.villageData(s2, T)
  
  # imputed  - observed  connections.
  s2$E = Eimp - s$E
  plot.villageData(s2, T)
  
  # imputed (complete connections)
  s2$E = Eimp
  plot.villageData(s2, T)
}

tau <- function(Yobs, Eobs, Z, X, G) {
  # Computes the value of the test statistic.
  #
  # Value: real number. A test statistic that is sensitive to violations 
  #   of the alternative.
  isolates = which(colSums(Eobs)==0)
  affected = which(colSums(Eobs) > 0)
  return(mean(Yobs[affected]) - mean(Yobs[isolates]))
}

rerandomize <- function(Z, Eobs, Eimp) {
  # Performs a rerandomization
  # TODO(ptoulis): Need to check that we are in the Z_0: set of 
  #   assignments defined by the null hypothesis.
  #
  # TODO(ptoulis): The following are important to clarify.
  # Rerandomization would be different under a different H0
  # For example if we were to test whether the outcomes 
  # of people in *exactly* two-step distance from seed 
  # differ from the outcomes of isolates, then we could 
  # not simply sample any 2 seeds. Instead we would need 
  # to restrict ourselves to randomizations where 
  # units change states. However, moving an isolate to being 
  # a seed, or a 1-step unit, should not happen, etc.
  new.seeds = sample(1:length(Z), size=sum(Z), replace=F)
  Znew = rep(0, length(Z))
  Znew[new.seeds] <- 1
  return(Znew)
}

new.outcomes <- function(Yobs, Zreal, Z.new) {
  # TODO(ptoulis): 
  #  Check whether the Z.new is ok with respect to Zreal
  # e.g., sum(Zreal) == sum(Z.new)
  return(Yobs)  # if the rerandomization is fine then the new outcomes are the same
                # because of the null hypothesis.
}

new.connections <- function(Znew, Eobs, Eimp) {
  # Compute new connections under the new treatment assignment.
  Ecom = Eimp
  if(any(Ecom > 1)) {
    stop("Missing E cannot overlap with observed (quick_hyp L:46). Was this intended?")
  }
  seeds = which(Znew==1)
  first.wave = seeds
  second.wave = c()
  third.wave = c()
  for(s in seeds) {
    # which units can be accessed in one-step  from <seed>?
    second.wave = union(second.wave, 
                        setdiff(which(Ecom[s, ]==1), first.wave))
  }
  for(s in second.wave) {
    third.wave = union(third.wave, 
                       setdiff(which(Ecom[s, ]==1), union(first.wave, second.wave)))
  }
  # Now include only those connections from seeds and immediates
  affected = union(first.wave, second.wave)
  others = setdiff(1:nrow(Ecom), affected)
  if(length(others) > 0) {
    Ecom[others, ] <- 0
  } else {
    warning("No other nodes, which is strange (quick_hypothesis:L57)")
  }
  Ecom[second.wave, first.wave] <- 0
  Ecom[third.wave, second.wave] <- 0
  return(Ecom)
}

test.newConnections <- function() {
  # TODO(ptoulis): Temp testing code. Should be moved elsewhere.
  s = sanitize.village.data(village.data(2))
  
  # 1 = original E
  par(mfrow=c(2, 2))
  par(mar=rep(0, 4))
  plot.villageData(s, T)
  
  # 2 = imputed.
  Eimp = with(s, impute.E(Y, E, Z, X, G))
  s2 = s
  s2$E = Eimp
  plot.villageData(s2, T)
  
  # 3 diff
  Ediff = Eimp - s$E - t(s$E)
  s2 = s
  s2$E = Ediff
  plot.villageData(s2, T)
  
  # 4 = new obs connections
  Z.new = rep(0, length(s$Y))
  Z.new[c(23, 28)] <- 1
  Enew = new.connections(Z.new, s$E, Eimp)

  s2$E = Enew
  s2$Z = Z.new
  plot.villageData(s2, T)
  
  # sanity check
  Z.new = rep(0, length(s$Y))
  Z.new[which(s$Z==1)] <- 1
  Enew = new.connections(Z.new, s$E, Eimp)
  print(all(Enew==s$E))
  i = which(Enew!=s$E)
  print(Enew[i])
  print(s$E[i])
}

test.H0 <- function(Yobs, Eobs, Zreal, X, G, nreps=100) {
  # Used to test one simple hypothesis: Whether the outcomes 
  #  of the units in the coupon tree are different from those 
  #   that are not. We call the former units "affected" and the latter 
  # "isolated".
  #
  # TODO(ptoulis): Do some sanity checks here.
  tobs <- NA
  tvals <- c()
  
  for(i in 1:nreps) {
    # 1. Impute the missing edges.
    Eimp = impute.E(Yobs, Eobs, Zreal, X, G)
    
    # 2. Calculate the test statistic
    Tobs = tau(Yobs, Eobs, Zreal, X, G)
    if(is.na(tobs)) {
      tobs <- Tobs
    } else {
      if(is.na(Tobs)) {
        warning("NA value in testing.")
      } else {
        tvals <- c(tvals, Tobs >= tobs)
        print(sprintf("Tobs = %.3f  initial=%.3f Current=%.2f%%", 
                      Tobs, tobs, 100 * mean(tvals)))
      }
    }
    
    # 3. Rerandomize
    Znew = rerandomize(Zreal, Eobs, Eimp)
    print("Rerandomizing...New seeds=")
    print(which(Znew==1))
    # 4. New outcomes
    Yobs.new = new.outcomes(Yobs, Zreal, Znew)
    
    # 5. New connections
    Eobs.new = new.connections(Znew, Eobs, Eimp)
    
    # 6. Repeat procedure
    Yobs <- Yobs.new
    Eobs <- Eobs.new
    Zreal <- Znew
  }
  
  return(tvals)
}