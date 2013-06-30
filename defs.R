# Panos Toulis, ptoulis@fas.harvard.edu
#
# Object definitions. 
# We assume we have 1...N experimental units (e.g. individuals)
# (0)  A "unit" is just the experimental unit (e.g. individual in a medical test)
#      We assume we have 1...N units and each one is represented by an integer in that range.
#      A "neighbor" is also a unit which is neighbor to the node (thus represented by id)
#      Unit j is a neighbor of i iff there is an incoming link j->i
#      "Vk"=set of units that have >=k neighbors.
#      An assignment "Zi" is in {0,1} and represents the treatment status of a unit i.
# (1) "experiment": list(Z=c(), exposed=c(), control=c(), can.expose=c(), can.control)
#     -------------------------------------------------------------------------------
#     Z=vector of assignments, Zi in {0,1}, len(Zi)=N (total #units) for every unit.
#     exposed/control=vector of unit ids which has have been exposed/controled.
#     exposed=all neighboring units of unit i have Zj=1
#     control=for neighboring units Zj=0 and also Zi=0
#     can.expose=vector of unit ids that can be exposed (e.g. neighboring nodes not assigned yet)
#     can.control=vector of unit ids that can be in control.
# (2) "cii": Causal inference with interference object. It is comprised of
#      {G=<igraph object>, Y=function(), k, experiment=<experiment>}
#      Y(i, Zni)=reponse of unit i when neighborhood has treatent Zni
#      k=level of effects we are considering as integer. When a unit i is exposed
#        to k-level effects, exactly k-neighbors are treated and the unit i is in control.
# (3) "cii.args": Arguments to construct a <cii> object (as a list). See above for the arguments.
source("tests.R")
empty.experiment <- function(num.units) {
  Zo <- rep(NA, num.units)
  obj <- list(Z=Zo, exposed=c(), control=c(), can.expose=c(), can.control=c())
  class(obj) <- "experiment"
  obj
}
empty.cii <- function() {
  experiment <- empty.experiment(0)
  graph <- graph.adjacency(matrix(0, nrow=0, ncol=0))
  obj <- list(G=graph, k=0, Y=function(i, Zni) {}, experiment=experiment)
  class(obj) <- "cii"
  obj
}
example.cii <- function() {
  args <- list()
  args$G <- watts.strogatz.game(dim=1, size=10, nei=2, p=0.3)
  args$Y <- function(i, Zi) { return(2 + i + mean(Zi, na.rm=T)) }
  args$k = 3
  cii <- init.cii(args)
  cii
}

# CHECK_* functions make sure these objects have the correct specs
CHECK_experiment <- function(experiment) {
  check.object(experiment, class.name="experiment",
               attrs=c("Z", "exposed", "control", "can.control", "can.expose"),
               str="Checking experiment object")
}
CHECK_cii <- function(cii) {
  # Check whether this is valid cii object
  check.object(cii, class.name="cii", attrs=c("G", "Y", "k", "experiment"))
  CHECK_experiment(cii$experiment)
}
CHECK_cii_units <- function(cii, units) {
  CHECK_cii(cii)
  if (length(units) ==0) 
    stop("Units list is empty. Possible error.")
  cii.units <- get.cii.units(cii)
  for (unit in units) {
    if (!(unit %in% cii.units))
      stop(sprintf("Unit %d does not belong in cii object (only %d units)", unit, length(cii.units)))
  }  
}
CHECK_Z <- function(z.list) {
  if (length(z.list) == 0)
    stop("Empty assignment vector")
  for (z in z.list) {
    if (!(z %in% c(0,1, NA))) {
      stop("Assignment z has to be in {0,1}")
    }
  }
}