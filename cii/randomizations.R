## Randomizations
##  Defines different ways to sample a network in order to get 
##  estimates of peer effects.
##  Author: Panos Toulis (ptoulis@fas.harvard.edu)
sample.node <- function(nodes) {
  # Sample one node id
  #
  # Args: An ARRAY of ids
  # Returns: -1 if empty, or one node at random 
  if (length(nodes) == 0)
    return (-1)
  else if (length(nodes) == 1)
    return(nodes)  # beware of the silly R sample() bug.
  else sample(nodes, size=1)
}
get.balance <- function(cii.randomiz) {
  
}
SSR <- function(cii, verbose=F) {
  # Simple Sequential Randomization design.
  #
  # Args: A <cii> object
  # Returns: A modified <cii> object after the randomization
  still.more = function() (length(which.can.control(cii)) +
                           length(which.can.expose(cii))>0)
  while(still.more()) {
    cs = which.can.control(cii)
    ts = which.can.expose(cii)
    # -1 if no control, or a random one otherwise
    control.unit = sample.node(cs)
    expose.unit = sample.node(ts)
    if (verbose) {
      log.vector(cs, "Control units")
      loginfo(sprintf("Control.unit=%d", control.unit))
      log.vector(ts, "Treat units")
      loginfo(sprintf("expose.unit=%d", expose.unit))
    }
    if(control.unit < 0 && expose.unit > 0) {
      cii = cii.expose.unit(cii, expose.unit)
      if (verbose) loginfo(sprintf("** Exposing unit %d", expose.unit))
    } else if (expose.unit < 0 & control.unit > 0) {
      cii = cii.control.unit(cii, control.unit)
      if (verbose) loginfo(sprintf("** Controlling unit %d", control.unit))
    } else {
      decide.to.treat = runif(1) < 0.5
      if (decide.to.treat) {
        cii <- cii.expose.unit(cii, expose.unit)
        if (verbose) loginfo(sprintf("** Exposing unit %d", expose.unit))
      }
      else {
        cii <- cii.control.unit(cii, control.unit)
        if (verbose) loginfo(sprintf("** Controlling unit %d", control.unit))
      }
    }
    if (verbose) {
      readline()
      plot.cii(cii)
    }
  }
  return(cii)
}

INR <- function(cii, verbose=F) {
  # Implements Insulated Neighborhood Randomization
  shared.neighbors = get.cii.shared.neighbors(cii)
  control.shared = rep(0, length(shared.neighbors))
  cii <- set.units.Z(cii, shared.neighbors, control.shared)
  
  still.more = function() (length(which.can.control(cii)) +
                             length(which.can.expose(cii))>0)
  mode = 1
  while (still.more()) {
    ts = which.can.expose(cii)
    cs = which.can.control(cii)
    if (mode == 1 && length(ts) > 0) {
      cii <- cii.expose.unit(cii, sample.node(ts))
    } else if (mode == -1 && length(cs) > 0) {
      cii <- cii.control.unit(cii, sample.node(cs))
    }
    mode = mode * (-1)
  }
  #cii <- SSR(cii, verbose=verbose)
  return (cii)
}

RAND <- function(cii) {
  all.units = get.cii.units(cii)
  Z = sample(c(0,1), size=length(all.units), replace=T)
  cii <- set.units.Z(cii, units=all.units, Zs=Z)
}
# DEPRECATED
bootstrap = function(x) sd(replicate(1000, { y = sample(x, replace=T); mean(y)}))
prez.table <- function(n, power=0.8, k.array, trials) {
  results = list(ssr=list(k=c(), mean=c(), se=c(), fn=SSR),
                 inr=list(k=c(), mean=c(), se=c(), fn=INR),
                 rand=list(k=c(), mean=c(), se=c(), fn=RAND))
  # Create cii object
  G = barabasi.game(n=n, power=power,directed=F)
  print(summary(degree(G)))
  for (k in k.array) {
    args = list(G=G, Y=NULL, k=k)
    cii <- init.cii(args)
    save(cii, file="cii.Rdata")
    loginfo(sprintf("%d units in Vk, k=%d, trials=%d",
                    length(get.cii.Vk(cii)),
                    k, trials))
    for (f in names(results)) {
      loginfo(sprintf("Running %s", f))
      fn = results[[f]]$fn
      vals = replicate(trials, {get.experiment.balance(fn(cii))})
      results[[f]]$k = c(results[[f]]$k, k)
      results[[f]]$mean=c(results[[f]]$mean, mean(vals))
      results[[f]]$se =c(results[[f]]$se, bootstrap(vals))
    }
  }
  return(results)
}