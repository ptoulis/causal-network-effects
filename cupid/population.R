## population.R
rm(list=ls())
source("../../r-toolkit/checks.R")
source("types.R")
# 
#  Code for a "population" object. Terminology is as follows.
#    population = LIST(obs=LIST(types, Z, Y), com=LIST(types, Z, Yall)  where:
#
#     obs = the population as it is "observed"
#        types = Nx1 vector of types (see types.R) i.e., males, females, singles etc.
#                 Whatever is observed, so it contains NA values
#                 It holds obs.types = types.observed(com.types, Z)
#                i.e. the observed types are deterministically dependent 
#                on the assignment and the true types.
#         Z = Nx1 vector of treatment assignment (0,1) --always observed
#           So Zobs = Zcom
#        Y = Nx1 vector of observed outcomes on the units (fully observed)
#      
#      com = the population as it really "is"
#        types = Nx1 vector of true types.
#        Z = Nx1 vector of assignment (equal to Zobs)
#        Yall = LIST(Nx4, Nx2) of matrices containing the PO
#           The Nx4 matrix contains Yi(, ) the potential outcomes of males.
#           The column order is Y(0, 0), Y(0, 1), Y(1, 0), Y(1, 1)
#           The Nx2 matrix  contains Yi() the potential outcomes of singles or females.
#      
# Critical functions
# - population.treatment.outcomes(pop, z) : Given the specified assignment 
#     it will read the correct potential outcomes from the Yall matrix of POs.
# - population.rerandomize(pop) : It will draw a new Z for the same PO
#      and update the observed portion of the population.
# - population.filter(pop, ..) : Will return the units that match the criteria.

# Generic information for a population object
population.size <- function(pop) length(pop$obs$Z)
population.all.units <- function(pop) 1:population.size(pop)

population.check.units <- function(pop, units) {
  # Check whether these are valid units.
  CHECK_MEMBER(units, population.all.units(pop))
}

# Basic access functions (to reduce list-referencing notational clutter).
population.treatment <- function(pop) {
  # Returns the treatment vector(Z)
  pop$obs$Z
}

population.types <- function(pop, obs=T) {
  if(obs) return(pop$obs$types)
  return(pop$com$types)
}

population.Y <- function(pop, obs=T) {
  # Gets the (sub)vector of the observed outcomes Y for the units.
  #
  # If obs=T then returns the Nx1 vector of observed outcomes
  # If obs=F then returns the Yall object (LIST of Nx4, Nx2) matrices of PO
  if(obs) return(pop$obs$Y)
  return(pop$com$Yall)  # returns the potential outcomes matrix.
}


population.filter <- function(pop, is.type=NA, has.treatment=NA, match.treatment=NA,
                              obs=T) {
  # Returns those unit ids that have the specified type,
  # and treatment and their match (spouse) has the specified treatment.
  #
  # Args: pop = population, is.type=("s", "m", "f"), has.treatment=(0, 1), match=(0,1)
  # Returns:  vector of unit ids.
  #
  # Example:
  # pop = new.population(10)
  # population.filter(pop, is.type="m", has.treatment=1, match.treatment=0)
  #
  # = returns males that have z=1 and their spouses are not treated.
  if(!obs) {
    warning("Filtering population units based on non-observed types")
  }
  units = population.all.units(pop)
  z = population.treatment(pop)
  obs.types = population.types(pop, obs=obs)
  
  if(!is.na(is.type)) {
    CHECK_MEMBER(is.type, c("s", "m", "f"), msg="Valid type")
    if(is.type=="s") {
      units = intersect(units, types.singles(obs.types))
    } else if(is.type=="f") {
      units = intersect(units, types.females(obs.types))
    } else {
      units = intersect(units, types.males(obs.types))
    }
  }
  
  # 2. Keep those units under the treatment
  if(!is.na(has.treatment)) {
    units.z = population.all.units(pop)
    CHECK_MEMBER(has.treatment, c(0, 1), msg="Valid treatment")
    z.units = which(z==has.treatment)
    units = intersect(units, z.units)
  }
  
  # 3. Treatment of the match.
  if(!is.na(match.treatment)) {
    CHECK_MEMBER(match.treatment, c(0, 1), msg="Valid treatment")
    CHECK_MEMBER(is.type, c("m", "f"), msg="Only males/females have matches")
    
    other.units = which(z==match.treatment)
    if(is.type=="m") {
      other.units = intersect(other.units, types.females(types = obs.types))
      units = intersect(units, types.male.match(obs.types, other.units))
      # we need an exception when no singles.
      if(!is.na(has.treatment) & has.treatment==1 && 
          match.treatment==0 & !population.has.singles(pop)) {
        m11 = population.filter(pop, is.type="m", 
                                has.treatment=1, 
                                match.treatment = 1, 
                                obs=obs)
        m1 = population.filter(pop, is.type="m", has.treatment=1)
        return(union(units, setdiff(m1, m11)))
      }
    } else {
      other.units = intersect(other.units, types.males(types = obs.types))
      units = intersect(units, types.female.match(obs.types, other.units))
    }
  }
  
  return(units)
}

population.treatment.outcomes <- function(pop, z) {
  # Explicitly computes the outcomes for the specified treatment vector.
  # Given the assignment vector (z) it will use Yall (potential outcomes)
  # to return the Nx1 vector of observed outcomes on the units.
  #
  # Args: pop = population, z=Nx1 treatment vector
  # Returns: Nx1 vector of Y outcomes.
  #
  units = population.all.units(pop)
  N = length(units)
  CHECK_EQ(length(z), length(units)) # make sure z has the correct length
  CHECK_MEMBER(z, c(0, 1), msg="treatment should be binary") # check if treatment is binary
  
  all.types = population.types(pop, obs=F)
  males = types.males(all.types)
  match.females = types.female.match(all.types, males)
  other = setdiff(units, males)  # females and singles.
  
  # Will use linear indexing to access the elements of the matrix
  # of PO using (r, c) where r=row positions, c=column positions.
  male.c = 2 * z[males] + z[match.females] + 1 # from 1, 4
  male.r = males
  Lindex.male = male.r + N * (male.c-1)
  other.c = z[other] + 1
  other.r = other
  Lindex.other = other.r + N * (other.c-1)
  
  y = rep(NA, N)
  y[males] <- pop$com$Yall$males[Lindex.male]
  y[other] <- pop$com$Yall$other[Lindex.other]

  CHECK_TRUE(sum(is.na(y))==0)
  return(y)
}

# Sampling of the potential outcomes
# TODO(ptoulis): There are multiple models we can try.
sample.Yall <- function(types) {
  # Given the specified types and assignment, sample the potential outcomes.
  # TODO(ptoulis): It is important to define a NULL hypothesis to sample from.
  #  This is needed in order to check the validity of the methods.
  kBaselineEffect = 1.0
  kPrimaryEffect = 2.75

  kMalesPrimaryEffect = 1.5
  kMalesBaselineEffect = 0.5
  kCupidEffect = 4
  
  N = length(types) 
  y = list(males=matrix(NA, nrow=N, ncol=4), 
           other=matrix(NA, nrow=N, ncol=2))
  # For now we will assume no interference effect and *same* treatment effect 
  # for everyone.  TODO(ptoulis): Need to adapt this.
  m = types.males(types)
  f.and.s = union(types.females(types), types.singles(types))
  nFS = length(f.and.s)
  nM = length(m)
  # 1. Sample   Y_i(0)  for females and singles.
  #   E(Y(1)-Y(0)) = 2.75
  y$other[f.and.s, 1] <- rnorm(nFS, mean=kBaselineEffect, sd=1)
  # 2. Sample   Y_i(1)  for females and singles.
  y$other[f.and.s, 2] <- rnorm(nFS, mean=kBaselineEffect + kPrimaryEffect, sd=sqrt(2))
  
  mus = c(kMalesBaselineEffect, 
          kMalesBaselineEffect + kCupidEffect, 
          kMalesBaselineEffect + kMalesPrimaryEffect,
          kMalesBaselineEffect + kMalesPrimaryEffect + kCupidEffect)
  for(j in 1:4) {
    # 3. Sample Y_i(z1, z2)  for the males. Note that no interference effect
    #      means Y_i(0, 1) ~ Y_i(0, 0) and Y_i(1, 0) ~ Y_i(1, 1)
    #      so that, in the alt notation, Yi(z=1) ~ Yi(z=2) and  Yi(z=3) ~ Yi(z=4)
    # Assume treatment effect 1.9-0.6 = 1.3
    mu = mus[j]
    y$males[m, j] <- rnorm(nM, mean=mu, sd=0.6)
  }
  return(y)
}

sample.z <- function(nunits) {
  CHECK_TRUE(nunits%%2==0, msg="Better to have even units")
  return(sample(c(rep(0, nunits/2), rep(1, nunits/2))))
}

population.rerandomize <- function(pop) {
  # Samples a new treatment(z) for the population 
  # and the uses the potential outcomes matrix to update the 
  # observed values Yobs, and the true types to update the observed
  # unit types.
  #
  # Returns: A new population object.
  # Note: No sampling is performed for the new outcomes Y. These 
  # have been sampled fully when the population was created.
  #
  N = population.size(pop)
  # New randomization.
  z.new = sample.z(N)
  
  # 1. Update the Z vector
  pop$com$Z <- z.new
  pop$obs$Z <- z.new
  
  # 2. Update the observed types.
  true.types = population.types(pop, obs=F)
  pop$obs$types <- observed.types(true.types, z.new, 
                                  no.singles=pop$kNoSingles)
  
  # 3. Update the observed Y
  pop$obs$Y <- population.treatment.outcomes(pop, z.new)
  
  return(pop)
}
population.has.singles <- function(pop) {
  !pop$kNoSingles
}

new.population <- function(nUnits, singles.pct=0.0) {
  # Creates a new population object.
  #
  # 1. Sample the data
  types = sample.types(nUnits, singles.pct = singles.pct)
  Yall = sample.Yall(types=types)
  z = sample.z(nUnits)
  empty = rep(NA, nUnits)
  
  pop = list(com=list(types=types, Yall=Yall, Z=z),
             obs=list(types=empty, Z=empty, Y=empty))
  pop$kNoSingles = (singles.pct==0)
  pop <- population.rerandomize(pop)
  CHECK_population(pop)
  return(pop)
}

plot.population <- function(pop) {
  par(mfrow=c(2, 1))
  par(mar=rep(1, 4))
  plot.types(pop, plot.obs=F)
  plot.types(pop, plot.obs=T)
}
  

# This object represents a population of units.
CHECK_population <- function(pop) {
 
  # 1. Check the basics.
  CHECK_MEMBER(c("com", "obs", "kNoSingles"), names(pop), msg="Not the correct elements")
  CHECK_SETEQ(c("types", "Z", "Yall"), names(pop$com))
  CHECK_SETEQ(c("types", "Z", "Y"), names(pop$obs))
  # 2. Check Z (assignment) + lengths.
  N = length(pop$obs$types)
  CHECK_EQ(length(pop$com$Z), N)
  CHECK_EQ(length(pop$obs$Z), N)
  CHECK_EQ(length(pop$obs$Y), N)
  CHECK_EQ(nrow(pop$com$Yall$males), N)
  CHECK_EQ(nrow(pop$com$Yall$other), N)
  
  CHECK_SETEQ(pop$com$Z, pop$obs$Z, "All treatment is observed.")
  CHECK_MEMBER(pop$com$Z, c(0, 1), msg="Binary treatment")
  
  true.types = population.types(pop, obs=F)
  obs.types = population.types(pop, obs=T)
  Yall = population.Y(pop, obs=F)
  f = types.females(true.types)
  m = types.males(true.types)
  s = types.singles(true.types)
  
  # 3. Check types
  CHECK_types(true.types)
  # Check the consistency of obs vs complete types.
  z = population.treatment(pop)
  CHECK_consistent_types(types.com = true.types, types.obs = obs.types, 
                         z = z, no.singles = pop$kNoSingles)  
  
  # 4. Check Y (potential outcomes)
  CHECK_MEMBER(c("males", "other"), names(Yall))
  # only s+f are allowed NA values in the males potential outcomes matrix.
  CHECK_SETEQ(which(is.na(Yall$males[,1])), union(f, s), msg="No NA in complete data")
  # only m are allowed NA values in the (f+s) potential outcomes matrix.
  is.na.col = which(is.na(Yall$other[,1]))

  CHECK_MEMBER(m, is.na.col, msg="No NA in complete data")
  CHECK_TRUE(nrow(Yall$males)==N)
  CHECK_TRUE(ncol(Yall$males)==4)
  CHECK_TRUE(nrow(Yall$other)==N)
  CHECK_TRUE(ncol(Yall$other)==2)

  # 5. Check consistency of y.obs and y.com
  y.obs = population.Y(pop, obs=T)
  y.obs.valid = population.treatment.outcomes(pop, z=z)                                        
  CHECK_TRUE(all(y.obs==y.obs.valid))
}

test.population <- function() {
  # TODO(ptoulis): Add some tests here.
  pop = new.population(10)
  
}
