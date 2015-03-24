# Panos Toulis, ptoulis@fas.harvard.edu
rm(list=ls()) ####  remove workspace stuff
library(igraph)
library(logging)
basicConfig()
source("terminology.R")
source("randomizations.R")
init.cii  <- function(cii.args) {
  # Initializes the Causal Inference w/Interference object.
  # 
  # Args: A <cii.args> object. See "terminology"
  # Returns: A <cii> object. See terminology"
  cii <- empty.cii()
  # 1. Set the graph
  cii$G = cii.args$G 
  # 2. Response function Y(i, Zni)
  cii$Y = cii.args$Y 
  # 3.  Level of effects 
  cii$k = cii.args$k
  
  # 4. Update the "experiment" field (which are exposed, or can be etc)
  num.units = length(V(cii.args$G))
  cii$experiment <- empty.experiment(num.units)
  cii <- set.units.Z(cii, get.cii.units(cii), rep(NA, num.units))
  return(cii);
}

post.assignment.update <- function(cii) {
  # Updates can.expose/can.treat, exposed/control groups for the cii object.
  # Called after a new assignment is performed.
  #
  # Args: A <cii> object.
  # Returns: A modified <cii> object.
  eligible <- get.cii.Vk(cii)
  # Check all eligible units.
  for(unit in eligible) {
    ##  Check if already in  group or in treatment
    if(is.na(get.unit.Zi(cii, unit)) || is.unit.exposed(cii, unit) || is.unit.control(cii, unit)
         || is.unit.treated(cii, unit)) next;
    Zni <- get.unit.Znei(cii, unit)  # treatment of neighbors of i
    sum.z <- sum(Zni)  # sum of assignments
    sum.z.not.na <- sum(Zni[!is.na(Zni)])  # sum of the assigned
    # If exactly k neighbors treated, unit i is k-exposed
    if(sum.z.not.na == cii$k) cii$experiment$exposed = c(cii$experiment$exposed, unit)
    # 2. If < k in treatment and NAs exist, not sure
    if(is.na(sum.z)) next  ## there are neighbors in NA status.
    # 3. If all neighbors in control, unit is in k-control
    if(sum.z == 0) cii$experiment$control = c(cii$experiment$control, unit)
  }

  cii$experiment$can.expose <- which.can.expose(cii)
  cii$experiment$can.control <- which.can.control(cii)
  check.sets.disjoint(cii$experiment$exposed, cii$experiment$control,
                      str="Exposed and control groups should be disjoint")
  return(cii)
}

which.can.expose  <- function(cii) {
  # Calculates which units can be k-exposed in the <cii> object.
  #
  # Args: A <cii> object.
  # Returns: A vector of unit ids (in 1..N)
  CHECK_cii(cii)
  eligible <- get.cii.Vk(cii);
  ret.set <- rep(0,0) # to-be returned.
  # A node is eligible to be k-exposed iff
  # 1. It is eligible, unit in Vk
  # 2. Unit is not already assigned(control, k-exposure or treatment)
  # 3. #treated neighbors <= k AND #unassigned >= k - #treated
  for(unit in eligible) {
    zi <- get.unit.Zi(cii, unit)
    # Check if already exposed or treated
    # Condition 1.
    if(is.unit.exposed(cii, unit) || is.unit.treated(cii, unit) ||
       is.unit.control(cii, unit)) 
      next;
    # Get neighborhood assignment
    Zni <- get.unit.Znei(cii, unit)  # neighborhood assignment of unit
    num.unassigned <- length(which(is.na(Zni)));  # count units in point-control
    n1 <- length(which(Zni==1))  # which units are in treatment
    # Condition 2.
    if(n1 <= cii$k && num.unassigned >= (cii$k-n1))
      ret.set = c(ret.set, unit)
  }
  return(ret.set)
}
which.can.control <- function(cii) {
  # Calculates which units can be k-controlled in the <cii> object.
  #
  # Args: A <cii> object.
  # Returns: A vector of unit ids (in 1..N)
  eligible <- get.cii.Vk(cii);
  ret.set <- rep(0,0)
  # Conditions to be eligible for control
  # 1. Be in Vk
  # 2. Not in control already, or not assigned to Zi=1
  # 3. No neighbor has Zj=1
  # Condition 1. Take all eligible units.
  for(unit in eligible) {
    zi <- get.unit.Zi(cii, unit);
    ## Condition 2. Check if already assigned (see above)
    if(is.unit.treated(cii, unit) || is.unit.exposed(cii, unit) ||
       is.unit.control(cii, unit))
      next;
    Zni <- get.unit.Znei(cii, unit)
    # Condition 3. Not neighbor in Zj=1
    n1 <- length(which(Zni==1));  # count units in point-control
    if(n1 == 0) ret.set = c(ret.set, unit)
  }
  return(ret.set)
}

k.ones.of.n = function(k, n) {
  # Draw k "1"'s out of n total elements.
  if(k > n) 
    stop(sprintf("k=%d has to be <= n=%d", k,n))
  if(k == 0) return(rep(0, n))
  if(k == n) return(rep(1, n))
  sample(c(rep(1, k), rep(0, n-k)))  # sample a random permutation of k "1" and (n-k) "0"
}

cii.expose.unit <- function(cii, unit) {
  # Sets unit to be in k-exposure
  #
  # Args: <cii> object and the unit to be k-exposed
  # Returns: Modified <cii> object.
  CHECK_cii_units(cii, unit)
  ## Check to see if v is in treated.
  if(!(unit %in% cii$experiment$can.expose))
    stop(sprintf("Error: Node %d cannot be k-exposed..", unit))
  Ni <- get.unit.neighbors(cii, unit)
  Zni<- get.unit.Znei(cii, unit)  # get neighborhood assignment.
  Ni.nas <- Ni[which(is.na(Zni))]  # take the unassigned neighbors as a vector.
  n1 <- length(which(Zni==1))  # total #neighbors in Zj=1
  
  ##  Now we wish to treat  (k-n1) neighbors.
  size.nas = length(Ni.nas)
  to.treat = cii$k-n1
  # Throw errors if not enough unassigned neighbors
  # or if > k neighbors already assigned.
  # This is because we expect $experiment$can.expose to be valid.
  # Draw one random assignment (this throws the errors, if sizes not correct)
  assg <- k.ones.of.n(to.treat, size.nas)
  # Now change the Z assignment. 
  # Need to update can.expoase, can.control vector
  cii <- set.units.Z(cii, c(Ni.nas, unit), c(assg, 0))
  # Check whether all neighbors have been assigned.
  check.true((sum(is.na(get.unit.Znei(cii, unit)))==0), "All neighbors should be set")
  return(cii);
}

cii.control.unit <- function(cii, unit) {
  # Sets unit to be in k-control
  #
  # Args: <cii> object and the unit to be k-exposed
  # Returns: Modified <cii> object.
  ## Check to see if v is in can-control  
  if(!(unit %in% cii$experiment$can.control))     
    stop(sprintf("Error: Unit %d cannot be in control..", unit))
  
  Ni <- get.unit.neighbors(cii, unit)
  assignment <- rep(0, length(Ni))  # set all neighbors to control(Zj=0)
  cii <- set.units.Z(cii, c(Ni, unit), c(assignment, 0))
  # Check whether all neighbors are set to Zj=0
  check.lists.eq(get.unit.Znei(cii, unit), rep(0, length(Ni)), "All neighbors should be zero")
  return(cii);
}

## Response of a unit.
get.unit.response = function(cii, unit, Znei) {
  Ni.size <- length(get.unit.neighbors(cii, unit))
  check.true(Ni.size==length(Znei), "Treatment assignments in unit response should match.")
  cii$Y(unit, Znei)
}
# Return the units in the CII object
get.cii.units <- function(cii) as.numeric(V(cii$G))

get.unit.neighbors <- function(cii, unit) {
  # Returns neighbors of unit as a vector
  #
  # Args: <cii> object and unit id (integer)
  # Returns: Vector of unit ids j such that j->i is an edge
  nhood <- neighborhood(cii$G, 1, nodes=c(unit), mode="in")[[1]];
  return(nhood[-which(nhood==unit)]);
}

get.cii.Vk <- function(cii) {
  # Returns Vk = nodes that have >=k neighbors as a vector.
  # Also called the set of "eligibile" nodes.
  CHECK_cii(cii)
  all.units <- get.cii.units(cii);
  ind <- sapply(all.units, function(i) as.numeric(length(get.unit.neighbors(cii, i)) >= cii$k))
  return(all.units[ind==1]);
}

## Returns units assignment
get.unit.Zi <- function(cii, unit) {
  CHECK_cii_units(cii, unit)
  cii$experiment$Z[unit]
}

set.unit.Zi <- function(cii, unit, zi) {
  set.units.Z(cii, c(unit), c(zi))
}

is.unit.exposed <- function(cii, unit) (unit %in% cii$experiment$exposed)
is.unit.control <- function(cii, unit) (unit %in% cii$experiment$control)
is.unit.treated <- function(cii, unit) {
  zi <- get.unit.Zi(cii, unit)
  return(!is.na(zi) && zi==1)
}

get.cii.shared.neighbors <- function(cii) {
  all.units = get.cii.units(cii)
  which(sapply(all.units, function(i) is.shared.neighbor(cii, i)))
}
is.shared.neighbor <- function(cii, unit) {
  member = sapply(get.cii.Vk(cii),
                  function(i) unit %in% get.unit.neighbors(cii, i))
  return (sum(member) > 1)
}
can.expose.unit <- function(cii, unit) (unit %in% cii$experiment$can.expose)
can.control.unit <- function(cii, unit) (unit %in% cii$experiment$can.control)

get.units.Z <- function(cii, units) {
  CHECK_cii_units(cii, units)
  cii$experiment$Z[units]
}

set.units.Z <- function(cii, units, Zs) {
  CHECK_Z(Zs)
  current.Z <- get.units.Z(cii, units)
  check.lists.eq(length(current.Z), length(Zs), str="equal-length assignment")
  # which nodes are already assigned
  assigned.index = which(!is.na(current.Z))
  n = length(assigned.index)
  check.lists.eq(current.Z[assigned.index],
                 Zs[assigned.index],
                 str="Respecting prior assignments")
  cii$experiment$Z[units] <- Zs;
  cii <- post.assignment.update(cii)
  return(cii)
}

get.unit.Znei = function(cii, unit) {
  # Returns the assignment vector for the neighbors of unit i
  #
  # Args: <cii> object and the unit id.
  # Returns: An assignment vector Zni
  CHECK_cii_units(cii, unit)
  Ni <- get.unit.neighbors(cii, unit)
  Znei <- get.units.Z(cii, Ni)
  return(get.units.Z(cii, Ni))
}

plot.cii <- function(cii, layout=layout.circle) {
  # Plots the cii object. Different colors are used for the assignments.
  # blue=unassigned(NA), green=assigned, red=control
  par(mar=rep(0, 4))
  all.units <- get.cii.units(cii)
  toCodes <- sapply(get.units.Z(cii, all.units), function(i) ifelse(is.na(i), 3, i+1))
  num.units <- length(all.units)
  V(cii$G)$color <- c("red", "green", "cyan")[toCodes];
  plot.igraph(cii$G, layout=layout,
              vertex.size=floor(300/num.units))      ###  plot the graph
}

get.experiment.balance <- function(cii) {
  n = length(get.cii.Vk(cii))
  Nt = length(cii$experiment$exposed)
  Nc = length(cii$experiment$control)
  return (min(Nt, Nc) / (n/2))
}
