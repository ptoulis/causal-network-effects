# Copyright 2013 Panos Toulis, Donald B. Rubin
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Inference on treatment effects. Using hypothesis testing.
# par = parameters (df)
# mu.male, sigma.male, mu.female, sigma.female
#
# data = LIST(Yt, Yc, tc.true) -- same length.
#        tc.true = vector of true types.
source("population.R")
CHECK_par <- function(par) {
  CHECK_identical(names(par), c("mu.male", "sigma.male", 
                                "mu.female", "sigma.female"), msg="Valid par")
}

get.Y.types <- function(data, tc.imp) {
  # Will read the imputed vector and return 
  # a LIST(male, female) of the observed outcomes.
  #
  # f(Yt, Yc, C) = (male=..., female=...)
  males = which(tc.imp=="M")
  females = which(tc.imp=="F")
  CHECK_MEMBER(tc.imp, c("M", "F"))

  CHECK_TRUE(length(tc.imp)==length(data$Yc), msg="correct #imputations")
  
  return(list(male = c(data$Yt, data$Yc[males]),
              female = data$Yc[females]))
}

sample.par <- function(data, tc.imp, par.old) {
  # Will sample from  p(a | C, y)
  #
  CHECK_par(par.old)
  
  Y = get.Y.types(data, tc.imp)
  y.males = Y$male
  y.females = Y$female
  
  s2.males = var(y.males)
  s2.females = var(y.females)
  n.males = length(y.males)
  n.females = length(y.females)
  
  # Sample
  sigma2.male = 1 / rgamma(1, shape = n.males/2 + 1, rate = (n.males-1) * s2.males/2)
  mu.male = rnorm(1, mean=mean(y.males), sd=sqrt(sigma2.male / n.males))
  
  sigma2.female = 1 / rgamma(1, shape = n.females/2 + 1, rate = (n.females-1) * s2.females/2)
  mu.female = rnorm(1, mean=mean(y.females), sd=sqrt(sigma2.female/n.females))
  
  par.new = list(mu.male=mu.male, sigma.male=sqrt(sigma2.male),
                 mu.female=mu.female, sigma.female=sqrt(sigma2.female))
  CHECK_par(par.new)
  return(par.new)
}

complete.data.log.likelihood <- function(data, tc.imp, par) {
  # Given some model parameters "par" and data and an inputed type vector
  # compute the log-likelihood.
  #
  Y = get.Y.types(data, tc.imp)
  log.lik.male = sum(dnorm(Y$male, mean = par$mu.male, sd = par$sigma.male, log=T))
  log.lik.female = sum(dnorm(Y$female, mean=par$mu.female, sd = par$sigma.female, log=T))
  
  return(log.lik.male + log.lik.female)
}

complete.data.estimate <- function(data, tc.imp) {
  # Given an imputed type vector, compute the estimand.
  #
  Y = get.Y.types(data, tc.imp) # do this to check tc.imp
  males = which(tc.imp=="M")
  return(mean(data$Yt) - mean(data$Yc[males]))
}

init.par <- function(data) {
  # Initialize the model parameter.
  return(list(mu.male = mean(data$Yt), sigma.male=sd(data$Yt),
              mu.female = mean(data$Yc), sigma.female=sd(data$Yc)))
}

random.types <- function(Nc) {
  CHECK_TRUE(Nc %% 2==0, msg="even #units in the Control group")
  sample(c(rep("M", Nc/2), rep("F", Nc/2)))
}

perturb.types <- function(tc.imp, depth=2) {
  if(depth==0) return(tc.imp)
  m = sample(which(tc.imp=="M"), size=1)
  f = sample(which(tc.imp=="F"), size=1)
  # exchange between one male-female
  tc.imp[m] <- "F"
  tc.imp[f] <- "M"
  return(perturb.types(tc.imp, depth-1))
}

generate.data <- function(Nt, Nc, cupid.effect) {
  # Generates a data object
  #
  mu.male.c =  3.5  # male mean control
  mu.male.t = mu.male.c + cupid.effect  # male mean
 
  mu.f = 5.8  # female mean
  se.m = 2.0 ## keep it the same?
  se.f = 3.5
  
  par(mfrow=c(1, 1))
  # Sample male outcomes.
  Yt = rnorm(Nt, mean=mu.male.t, sd=se.m)
  Yc = rnorm(Nc, mean=mu.male.c, sd=se.m)
  # Sample female outcomes.
  tc.true = random.types(Nc)
  
  females = which(tc.true=="F")
  Yc[females] = rnorm(length(females), mean=mu.f, sd=se.f)

  return(list(Yt=Yt, Yc=Yc, tc.true=tc.true, true.CACE=mu.male.t-mu.male.c))
}

#### TESTING.
test.likelihood <- function() {
  # Shows the likelihood ratios around the true value.
  d = generate.data(300, 400)
  tc =  d$tc.true # true types.
  p = init.par(d)
  L0 = complete.data.log.likelihood(d, tc, p)
  print(sprintf("True-data likelihood = %.3f", L0))
  
  mh.ratio <- function(ll.values) {
    # need to show min(1, exp(ll-L0))  values
    sapply(ll.values, function(i) min(1, exp(i - L0)))
  }
  
  x = replicate(100, { complete.data.log.likelihood(d, perturb.types(tc, 1), p)})
  par(mfrow=c(2, 2))
  hist(mh.ratio(x), main="Histogram logL(perturb=1)", xlim=c(0, 1))
  x = replicate(100, { complete.data.log.likelihood(d, perturb.types(tc, 2), p) })
  hist(mh.ratio(x), main="Histogram logL(perturb=2)", xlim=c(0, 1))
  x = replicate(100, { complete.data.log.likelihood(d, perturb.types(tc, 5), p) })
  hist(mh.ratio(x), main="Histogram logL(perturb=5)", xlim=c(0, 1))
  x = replicate(100, { complete.data.log.likelihood(d, random.types(400), p) })
  hist(mh.ratio(x), main="Histogram logL(random)", xlim=c(0, 1))
}

## Inference methods start here.
#
# Simple MH
#
correct.types.ratio <- function(data, tc.imp) {
  # Computes the ratio of correct types in the imputed vector.
  #
  length(which(tc.imp==data$tc.true)) / length(data$Yc)
}

mcmc.mh <- function(data, niters=1000, verbose=T) {
  
  K = length(data$Yc)
  # Initialization
  types = random.types(K)  # holds the imputed vector of types.
  types.matrix = matrix(NA, nrow=niters, ncol=K)
  nacc = 0  # acceptance.
  par = init.par(data)
  est = c()
  plot.times = as.integer(seq(10, niters, 
                              length.out=max(1, as.integer(niters/1000))))
  
  true.est = complete.data.estimate(data, data$tc.true)
  print(sprintf("True complete-data estimate = %.3f", true.est))
  for(iter in 1:niters) {
    new.types <- perturb.types(types, depth=1)
    
    new.ll = complete.data.log.likelihood(d, new.types, par)
    old.ll = complete.data.log.likelihood(d, types, par)
    if(verbose) {
      print(sprintf("OLD: ll=%.2f CR=%.1f%% -- NEW: ll=%.2f CR=%.1f%% -- diff=%.2f ll.ratio=%.2f",
                    old.ll, 100 * correct.types.ratio(data, types),
                    new.ll, 100 * correct.types.ratio(data, new.types),
                    new.ll - old.ll,
                    exp(new.ll - old.ll)))
    }
    if(iter %in% plot.times) {
      plot(est, type="l", ylim=c(true.est-1, true.est+1))
      abline(h=true.est, col="red")
    }
    log.acc = min(0, new.ll - old.ll)
    if(runif(1) < exp(log.acc)) {
      nacc <- nacc + 1
      types <- new.types
    }
    # 1b. update the types matrix
    types.matrix[iter, ] <- types=="M"
    
    # 2. Sample parameter
    par = sample.par(data, tc.imp = types, par.old = par)
    
    # new estimate
    est = c(est, complete.data.estimate(data, types))
  }
  print("Chain information")
  avg.types = as.numeric(colMeans(types.matrix) > 0.5) + 1
  avg.types = c("F", "M")[avg.types]
  # print(table(avg.types))
  print(sprintf("Correct types in chain %.1f%%", 100 * correct.types.ratio(data, avg.types)))
  print(sprintf("MH acceptance ratio %.2f%%", 100 * nacc/ niters))
  return(tail(est, 0.2 * length(est)))
}

EM.data <- function(data, tol=1e-2, maxIters=1000, verbose=F) {
  ## Computes MLE estimate of CACE for a simple normal model.
  ## (assumes common variance)
  # 
  # data = (Yt, Yc) where Yt = outcomes of Treated males.
  #                       Yc = outcomes of males or females
  Yt = data$Yt
  Nt = length(Yt)
  Yc = data$Yc
  #        T:male   C:male-mean   C:male-sd   C:female     C:female-sd
  # par = ( mu1t,     mu1c,         ms1,        mu0, s0)
 
  get.Q <- function(par.old) {
    # Returns Q(theta, theta_t)
    mu.male.t = par.old[1]
    mu.male.c = par.old[2]
    sigma.male = par.old[3]
    mu.female.c = par.old[4]
    sigma.female.c = par.old[5]
    
    male.lik = dnorm(Yc, mean=mu.male.c, sd=sigma.male)
    female.lik = dnorm(Yc, mean=mu.female.c, sd=sigma.female.c)
    Ci.imp = male.lik / (male.lik + female.lik)  # impute expected.
    
    Q <- function(par.new) {
      mu.male.t = par.new[1]
      mu.male.c = par.new[2]
      sigma.male = par.new[3]
      mu.female.c = par.new[4]
      sigma.female.c = par.new[5]
      
      dnorm(mean(Yt), mean=mu.male.t, sd=sigma.male / sqrt(Nt), log=T) + 
        sum(Ci.imp * dnorm(Yc, mean=mu.male.c, sd=sigma.male, log=T)) + 
          sum((1-Ci.imp) * dnorm(Yc, mean=mu.female.c, sd=sigma.female.c, log=T))
    }
    return(Q)
  }
  
  ## Iterate.
  par = c(mean(Yt), mean(Yc), sd(Yt), mean(Yc), sd(Yc)) # current estimate
  
  for(i in 1:maxIters) {
    # E-step
    Q = get.Q(par)
    # M-step
    par.old = par
    par = optim(par = par.old, fn = Q, lower=c(-Inf, -Inf, 1e-5, -Inf, 1e-5),
                method="L-BFGS-B", control=list(fnscale=-1))$par
    if(verbose)
      print(sprintf("Iteration %d/%d : Params=%s", i, maxIters, paste(round(par, 3), collapse=", ")))
    
    if(sum(abs(par - par.old)) < tol)
      return(par)
  }
  return(par)
}

test.EM <- function() {
  data = generate.data(500, 900)
  par = EM.CACE(data, tol=1e-5, maxIters = 1000, verbose = T)
  print(sprintf("MLE of CACE=%.3f -- True = %.3f", 
                par[1]-par[2],
                data$true.CACE))
}

##  Do inference.
impute.types <- function(data) {

  par.mle = EM.data(data, verbose=F)
  
  mu.male.c = par.mle[2]
  sigma.male = par.mle[3]
  mu.female.c = par.mle[4]
  sigma.female.c = par.mle[5]
  
  Yc=  data$Yc
  male.lik = dnorm(Yc, mean=mu.male.c, sd=sigma.male)
  female.lik = dnorm(Yc, mean=mu.female.c, sd=sigma.female.c)
  Ci.prob = male.lik / (male.lik + female.lik)  # impute expected.
  Ci = as.numeric(runif(length(Yc)) < Ci.prob)
  tc.imp = c("F", "M")[Ci + 1]
  # print(correct.types.ratio(data, tc.imp))
  CHECK_TRUE(length(tc.imp)==length(Yc))
  return(tc.imp)
}

rerandomize.data <- function(data, tc.imp) {

  Yt.old = data$Yt
  Yc.old = data$Yc
  CHECK_MEMBER(tc.imp, c("F", "M"))
  males = which(tc.imp=="M")
  females = which(tc.imp=="F")
  Y.males = c(Yt.old, Yc.old[males])
  Y.females = Yc.old[females]
  
  # index = (1, 1, 1, ... 0, 0, 0...)  first are males in T and so on.
  index = c(rep(1, length(Yt.old)), rep(0, length(males)))
  # 1. Shuffle males.
  new.index = sample(index)
  treated.males = which(new.index==1)
  control.males = which(new.index==0)
  
  Yt.new <- c(Y.males[treated.males])
  Yc.new <- c(Y.males[control.males], Y.females)
  CHECK_TRUE(length(Yc.new)==length(Yc.old))
  CHECK_TRUE(length(Yt.new)==length(Yt.old))
  
  
  return(list(Yt=Yt.new, Yc=Yc.new))
}

run.most.powerful.test <- function(data, niters=100, verbose=F) {
  original.data = data
  get.effect <- function(par) par[1]-par[2]
  
  S.obs = get.effect(EM.data(original.data))
  if(verbose) {
    print(sprintf("Observed test statistic = %.3f", S.obs))
  }
  S <- c()
  info.times = seq(1, niters, length.out=10)
  for(i in 1:niters) {
    types.imputed = impute.types(data)
    data = rerandomize.data(data, types.imputed)
    S.new.obs = get.effect(EM.data(data))
    S <- c(S, S.new.obs >= S.obs)
    if(verbose)
      print(sprintf("i=%d/%d p-value=%.3f S.new=%.2f S.obs=%.2f", 
                    i, niters, mean(S), S.new.obs, S.obs))
  }
   return(mean(S))
}

plot.pvalues <- function(cupid.effect, num.pvalues=100) {
  pvals = c()
  for(j in 1:num.pvalues) {
    d = generate.data(Nt=200, Nc=350, cupid.effect = cupid.effect)
    pvals = c(pvals, run.most.powerful.test(d, niters = 100, verbose = F))
    hist(pvals)
  }
  hist(pvals)
}
