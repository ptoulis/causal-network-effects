# Copyright 2013 Panos Toulis, Donald B. Rubin
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Inference on treatment effects. Using hypothesis testing.
# par = parameters (df)
# mu.male, sigma.male, mu.female, sigma.female
#
source("population.R")
CHECK_par <- function(par) {
  CHECK_identical(names(par), c("mu.male", "sigma.male", 
                                "mu.female", "sigma.female"), msg="Valid par")
}

get.Y.types <- function(data, tc.imp) {
  # Will read the imputed vector and return 
  # a LIST(y.males, y.females) of the observed outcomes.
  males = which(tc.imp=="M")
  females = which(tc.imp=="F")
  CHECK_MEMBER(tc.imp, c("M", "F"))
  # CHECK_TRUE(length(males)==length(females), msg="#m= #f")
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

complete.data.estimate <- function(data, tc.imp) {
  Y = get.Y.types(data, tc.imp) # do this to check tc.imp
  males = which(tc.imp=="M")
  return(mean(data$Yt) - mean(data$Yc[males]))
}

init.par <- function(data) {
  return(list(mu.male = mean(data$Yt), sigma.male=sd(data$Yt),
              mu.female = mean(data$Yc), sigma.female=sd(data$Yc)))
}

random.types <- function(Nc) {
  sample(c(rep("M", Nc/2), rep("F", Nc/2)))
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


generate.dataset <- function(Nt, Nc) {
  mu.m = 3.5
  mu.f = 10.5
  se.m = 1.5
  se.f = 2.5
  
  par(mfrow=c(1, 1))
  
  Yt = rnorm(Nt, mean=mu.m, sd=se.m)
  Yc = rnorm(Nc, mean=mu.m, sd=se.m)
  
  tc.true = random.types(Nc)
  
  females = which(tc.true=="F")
  Yc[females] = rnorm(length(females), mean=mu.f, sd=se.f)
  return(list(Yt=Yt, Yc=Yc, tc.true=tc.true))
}

test.mcmc <- function(niters.mcmc=1000) {
  D = generate.dataset(Nt=250, Nc=500)
  p2 <- hist(D$Yc, freq=F, breaks=20, col=rgb(0, 0, 0.5, 0.1))      
  p1 <- hist(D$Yt, freq=F, breaks=20, col=rgb(1,0, 0, 0.2), add=T)   
  
  print(sprintf("True estimate %.3f", 
                complete.data.estimate(D, tc.imp = D$tc.true)))
  est = new.statistic.mcmc(D,  niters=niters.mcmc)
  print(summary(est))
  plot.new()
  frame()
  m = mcmc(est)
  plot(m)
  return(m)
}


test.statistic.nonparametric <- function(Yt, Yc) {
  y1 = mean(Yt)
  ordered.control.units = order(abs(Yc-y1))
  
}

cupid.effects <- function(pop) {
  
}




