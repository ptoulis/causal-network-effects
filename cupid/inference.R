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

generate.data <- function(Nt, Nc) {
  # Generates a data object
  #
  mu.m = 3.5  # male mean
  mu.f = 5.5  # female mean
  se.m = 1.5
  se.f = 2.5
  
  par(mfrow=c(1, 1))
  # Sample male outcomes.
  Yt = rnorm(Nt, mean=mu.m, sd=se.m)
  Yc = rnorm(Nc, mean=mu.m, sd=se.m)
  # Sample female outcomes.
  tc.true = random.types(Nc)
  
  females = which(tc.true=="F")
  Yc[females] = rnorm(length(females), mean=mu.f, sd=se.f)

  return(list(Yt=Yt, Yc=Yc, tc.true=tc.true))
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


new.statistic.mcmc <- function(data, niters=1000) {
  Yt = data$Yt
  Yc = data$Yc
  K = length(Yc)
  
  types = random.types(length(Yc))
  par = init.par(data)
  est = c()
  types.matrix = matrix(NA, nrow=niters, ncol=K)
  
  for(iter in 1:niters) {
    #male probs = Prob(Data | Ci=m, C-i)
    Pd.male = dnorm(Yc, mean=par$mu.male, sd=par$sigma.male)
    Pd.female = dnorm(Yc, mean=par$mu.female, sd=par$sigma.female)
    male.probs = Pd.male / (Pd.male + Pd.female)
    male.ind = rbinom(K, size=1, prob=male.probs)
    # Relaxed Gibbs for types
    types = rep("F", K)
    types[which(male.ind==1)] <- "M"
    types.matrix[iter, ] <- types=="M"
    # Sample parameter
    par = sample.par(data, tc.imp = types, par.old = par)
    
    # new estimate
    est = c(est, complete.data.estimate(data, types))
  }
  
  avg.types = as.numeric(colMeans(types.matrix) > 0.5) + 1
  avg.types = c("F", "M")[avg.types]
  print(table(avg.types))
  print(sprintf("Correct %.1f%%", 100 * length(which(avg.types==data$tc.true)) / K))
  return(tail(est, 0.2 * length(est)))
}




