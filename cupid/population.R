# Copyright 2013 Panos Toulis, Donald B. Rubin
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
rm(list=ls())
source("../../r-toolkit/checks.R")

#      
# Critical functions
# - new.population(N, singles) -- creates a new population
# - population.filter(pop) -- returns the units with the specified search criteria.
# - population.Yobs(pop, Y) -- returns what portion of the PO (Y) will be observed.
# - population.types.obs(pop) -- returns the vector of observed types

kTypes = c("F", "M", "S")
kTypesHaveMatches = c("F", "M")
other.type <- function(type) {
  if(type=="F") return("M")
  else if(type=="M") return("F")
  return("S")
}

new.population <- function(N, single.pct=0) {
  CHECK_TRUE(N%%2==0, msg="Better to use even number of units.")
  unit.ids = 1:N
  M = matrix(NA, nrow=N, ncol=6)
  colnames(M) <- c("id", "type", "z", "match.id", "match.type", "match.z")
  M = as.data.frame(M)
  M$id = unit.ids
  M$z = sample(c(rep(0, N/2), rep(1, N/2)))
  
  # 1. Update singles.
  nSingles = as.integer(single.pct * N)
  if(nSingles %% 2 == 1) nSingles = nSingles - 1  #make it even
  single.units = sample(unit.ids, size=nSingles)
  M$type[single.units] <- "S"
  M$match.id[single.units] <- single.units
  M$match.type[single.units] <- "S"
  M$match.z[single.units] <- M$z[single.units]
  
  # 2. Update the males/females
  unit.ids = setdiff(unit.ids, single.units)
  males = sample(unit.ids, size=length(unit.ids)/2)
  females = sample(setdiff(unit.ids, males)) # shuffle
  # 2a. Update males
  M$type[males] <- "M"
  M$match.id[males] <- females
  M$match.type[males] <- "F"
  M$match.z[males] <- M$z[females]
  
  # 2b. Update females.
  M$type[females] <- "F"
  M$match.id[females] <- males
  M$match.type[females] <- "M"
  M$match.z[females] <- M$z[males]
 return(M) 
}
population.no.singles <- function(pop) {
  nrow(subset(pop, type=="S"))==0
}

population.filter <- function(pop, 
                              is.type=c(kTypes), 
                              has.treatment=c(0, 1), 
                              match.treatment=c(0, 1)) {
  # Returns those unit ids that have the specified type,
  # and treatment and their match (spouse) has the specified treatment.
  #
  # WARNING: This is using information that might not be observed.
  # For example, the type or the match assignment are not always observed.
  if(identical(is.type, c("S")) & !identical(match.treatment, has.treatment)) {
    warning("Querying for single will sent match.z == z in the filter criteria")
    match.treatment = has.treatment
  }
  units = subset(pop, (type %in% is.type) & (z %in% has.treatment) & (match.z %in% match.treatment))$id 
  return(units)
}

population.obs <- function(pop) {
  # Returns the observed population.
  types.obs = population.types.obs(pop)
  pop.obs = pop
  pop.obs$type = types.obs
  na.types = which(is.na(types.obs))
  pop.obs[na.types, c("match.id", "match.type", "match.z")] <- NA
  # EXCEPTION:
  # If a male is treated and we know its type, but its female is not
  # we probably do not know its matched female
  treated.m.control.f = population.filter(pop, is.type="M", has.treatment = 1, match.treatment = 0)
  ided.treat.m.control.f = intersect(treated.m.control.f, which(!is.na(types.obs)))
  if(length(ided.treat.m.control.f) > 0) {
    pop.obs[ided.treat.m.control.f, c("match.id", "match.type", "match.z")] <- NA
  }
  return(pop.obs)
}

population.types.obs <- function(pop) {
  # Given the set of types and assignment, it defines what are the observed types.
  #  If there are singles, the only observed types are for those m/f pairs
  #  where the female is treated.
  #  If we assume no singles, we additionally observe males who are treated
  #   but their females are in control.
  #
  # TODO(ptoulis): Add checks?
  # TODO(ptoulis): This is incomplete. We can infer that some nodes are singles
  #   for example, if they are treated and there are no unidentified 
  #   nodes in the control.
  #
  # WARNING(ptoulis): Use this function with caution. It is not perfected yet.
  #   
  true.types = pop$type
  CHECK_notNA(true.types, msg="True types have no NAs")
  # vector of observed types.
  obs.t = rep(NA, length(true.types))
  # 1. Observe treated females.
  treated.f <- population.filter(pop, is.type = "F", has.treatment = 1)
  obs.t[treated.f] <- true.types[treated.f]
  # 2. Observe their matched males
  m.with.treated.f <- population.filter(pop, is.type = "M", match.treatment = 1)
  obs.t[m.with.treated.f] <- true.types[m.with.treated.f]
  
  if(population.no.singles(pop)) {
    treated.m.with.control.f = population.filter(pop, is.type = "M", 
                                                 has.treatment = 1,
                                                 match.treatment = 0)
    obs.t[treated.m.with.control.f] <- true.types[treated.m.with.control.f]
  }
  return(obs.t)
}

population.Yobs <- function(pop, Y) {
  # Explicitly computes the outcomes for the specified treatment vector.
  # Given the assignment vector (z) it will use Yall (potential outcomes)
  # to return the Nx1 vector of observed outcomes on the units.
  #
  # Args: pop = population, z=Nx1 treatment vector
  # Returns: Nx1 vector of Y outcomes.
  #
  CHECK_EQ(nrow(pop), nrow(Y)) 
  
  males = population.filter(pop, is.type = "M")
  matched.females = pop$match.id[males]
  other = population.filter(pop, is.type=c("F", "S"))
  N = nrow(pop)  # no.units
  # Will use linear indexing to access the elements of the matrix
  # of PO using (r, c) where r=row positions, c=column positions.
  male.c = 2 * pop$z[males] + pop$z[matched.females] + 1 # from 1, 4
  male.r = males
  Lindex.male = male.r + N * (male.c-1)
  # Other units.
  other.c = 3 * pop$z[other] + 1
  other.r = other
  Lindex.other = other.r + N * (other.c-1)
  
  y = rep(NA, N)
  y[males] <- Y[Lindex.male]
  y[other] <- Y[Lindex.other]

  CHECK_TRUE(sum(is.na(y))==0)
  return(y)
}

sample.Y <- function(pop) {
  # Sample potential outcomes
  # TODO(ptoulis): These need to be sampled according to a hypothesis.
  # Effects are listed as (baseline-primary)
  effects = list("F"=c(1, 2.75), "M"=c(0.5, 1.5), "S"=c(0.35, 0.1))
  effects.se = list("F"=c(1.0), "M"=c(0.7), "S"=c(0.25))
  kCupidEffect = 1.15
  
  Y = matrix(NA, nrow=nrow(pop), ncol=4)
  # columns correspond to (0, 0), (0, 1), (1, 0) and (1, 1) for males.
  #                        (0)       x      x         (1) for singles and females.
  
  for(type in kTypes) {
    units.type = population.filter(pop, is.type = type)
    nUnits = length(units.type)
    se = effects.se[[type]]
    baseline = effects[[type]][1]
    primary =  effects[[type]][2]
   # print(sprintf("Sampling Y: type=%s baseline=%.2f, primary=%.3f", type, baseline, primary))
    if(type=="M") {
      Y[units.type, 1] <- rnorm(nUnits, mean=baseline, sd=se)
      Y[units.type, 2] <- rnorm(nUnits, mean=baseline + kCupidEffect, sd=se)
      Y[units.type, 3] <- rnorm(nUnits, mean=baseline + primary, sd=se)
      Y[units.type, 4] <- rnorm(nUnits, mean=baseline + primary + kCupidEffect, sd=se)
    } else {
      Y[units.type, 1] <- rnorm(nUnits, mean=baseline, sd=se)
      Y[units.type, 4] <- rnorm(nUnits, mean=baseline + primary, sd=se)
    }
  }
  return(Y)
}

# Drawing stuff for visualization
population.colors <- function(pop, ref.types) {
  all.cols <- c("pink", "red", 
                rgb(0, 0.1, 0.7, 0.15), rgb(0, 0.1, 0.7, 0.65),
                rgb(0, 1, 0.1, 0.2), rgb(0, 1, 0.1, 0.8))
  cols <- rep(NA, nrow(pop))
  not.id.types = which(is.na(ref.types))
  cols[intersect(which(pop$z==0), not.id.types)] <- rgb(1, 1, 1, 1)
  cols[intersect(which(pop$z==1), not.id.types)] <- rgb(0.5, 0.5, 0.5, 0.8)
  for(z in c(0, 1)) {
    for(type in kTypes) {
      units = intersect(which(ref.types==type), population.filter(pop, has.treatment = z))
      col.index = (match(type, kTypes)-1) * 2 + z +1
      cols[units] <- all.cols[col.index]
    }
  }
  return(cols)
}
plot.types <- function(pop, plot.obs=F) {
  require(igraph)
  adjlist = list()
  units = pop$id
  
  # 1. Created the edge list.
  for(i in units) {
    adjlist[[i]] <- rep(1, 0)
    draw.edge = (pop$type[i] == "F") & (!plot.obs | (plot.obs & pop$z[i]==1))
    if(draw.edge) adjlist[[i]] <- pop$match.id[i]      
  }
  g = graph.adjlist(adjlist, mode = c("out"), duplicate = TRUE)
  V(g)$name = sapply(units, function(i) sprintf("%d(%d)", i, pop$z[i]))
  if(plot.obs) {
    V(g)$color =population.colors(pop, population.types.obs(pop))
    plot(g, vertex.size=300/length(units), main="Observed")  
  } else {
    V(g)$color = population.colors(pop, pop$type)
    plot(g, vertex.size=300/length(units), main="Complete")
  }
}
plot.population <- function(pop) {
  par(mfrow=c(2, 1))
  par(mar=rep(1, 4))
  plot.types(pop, plot.obs=F)
  plot.types(pop, plot.obs=T)
}
