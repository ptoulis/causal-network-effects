# Inference on treatment effects. Using hypothesis testing.

posterior.mean <- function(y1, y0) {
  # Given observed outcomes and male-females labels y1, C1
  # find the posterior mean of the set y0 where the labels are missing.
  CHECK_TRUE(length(y0) %% 2 ==0, msg="should be even because has only f-m")
  
  log.com.conditional <- function(t0.imp, par) {
    # Computes  p(types | y)
    CHECK_TRUE(length(t0.imp) == length(y0))
    males = which(t0.imp=="M")
    females = which(t0.imp=="F")
    if(length(males) != length(females)) {
      return(-Inf) # M-F females should have equal numbers.
    }
    
    y.males = c(y1, y0[males])
    y.females = y0[females]
    lik.male = dnorm(y.males, mean=par$mu.male, sd=par$sigma.male, log = T)
    lik.female = dnorm(y.females, mean=par$mu.female, sd=par$sigma.female, log = T)
    return(lik.male + like.female)
  }
  
  sample.par <- function(t0.imp) {
    # TODO
  }
  
  complete.data.estimate <- function(t0.imp) {
    males = which(t0.imp=="M")
    females = which(t0.imp=="F")
    if(length(males) != length(females)) {
      return(-Inf) # M-F females should have equal numbers.
    }
    mean(y1) - mean(y0[males])
  }
  
  sample.types <- function() {
    sample(c(rep("M", length(y0)/2), rep("F"), length(y0)/2))
  }
  
  t0.old = sample.types()
  nacc = 0
  estimates = c()
  for(i in 1:1000) {
    t0.new = sample.types()
    log.acc = min(0, log.com.conditional(t0.new) - log.com.conditional(t0.old))
    if(runif(1) < exp(log.acc)) {
      t0.old = t0.new
      nacc = nacc + 1
    }
    estimates = c(estimates, complete.data.estimate(t0.old))
  }
  return(tail(estimates, 0.9 * length(estimates)))
}


test.posterior.mean <- function() {
  mu = 0
}


cupid.effects <- function(pop) {
  
}