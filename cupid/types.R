# Functions for types (we are using gender terminology)
#
# The current convention is as follows:
# types = numeric vector such that if:
#     t[i] = j  and i!=j, then i=female and j=male
#     t[i] = i then i = single
#     t[i] = NA, we "do not know the type"
#     t[i] = 0, this is a male
sample.types <- function(nunits, singles.pct=0) {
  # Samples a vector of types following the convention shown above.
  #
  # Args: nunits = #units(=N)
  #       singles.pct = proportion (0-1) of the singles in the population
  # Returns: Vector (Nx1) of types.
  all = 1:nunits
  CHECK_INTERVAL(singles.pct, min=0, max=1)
  nsingles = as.integer(nunits * singles.pct)
  is.even = function(x) x%%2==0
  nsingles = ifelse(is.even(nunits-nsingles), nsingles, nsingles+1)
  # 1. Sample single units
  singles = sample(all, size=nsingles, replace=F)
  
  left = setdiff(all, singles)
  CHECK_TRUE(is.even(length(left)), msg="m/f population should be even.")
  # 2. Sample females and males (half=females)
  females = sample(left, size=length(left)/2, replace=F)
  males = setdiff(left, females)
  
  # 3. Define the vector of types (no "NA")
  types = 1:nunits
  types[females] <- males
  types[males] <- 0
  return(types)
}

types.singles <- function(types) which(1:length(types)==types)
types.males <- function(types) which(types==0)
types.na <- function(types) which(is.na(types))
types.females <- function(types) {
  not.na = setdiff(1:length(types), types.na(types))
  # females = all_not_NA - males - singles
  possible.female = which(types!=1:length(types) & types!=0)
  intersect(not.na, possible.female)
}

types.male.match <- function(types, females) {
  # Returns the males that match the specified females
  #
  if(length(females)==0) return(c())
  all.f <- types.females(types)
  CHECK_MEMBER(females, all.f, msg="Should be subset of all females")
  types[females]
}

types.female.match <- function(types, males) {
  # Returns the males that match the specified females
  #
  if(length(males)==0) return(c())
  all.m <- types.males(types)
  CHECK_MEMBER(males, all.m, msg="Should be subset of all males")
  f = match(males, types)
  CHECK_TRUE(sum(is.na(f))==0, msg="No NAs")
  return(f)
}

observed.types <- function(true.types, z, no.singles=F) {
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
  CHECK_types(true.types)
  CHECK_notNA(true.types, msg="True types have no NAs")
  # vector of observed types.
  obs.t = rep(NA, length(true.types))
  # 1. Observe treated females.
  treated.f <- intersect(types.females(true.types), which(z==1))
  obs.t[treated.f] <- true.types[treated.f]
  # 2. Observe their matched males
  m.with.treated.f <- types.male.match(true.types, treated.f)
  obs.t[m.with.treated.f] <- 0
  
  if(no.singles) {
    ## more types are observed under the assumption of no singles.
    CHECK_TRUE(length(types.singles(true.types))==0, msg="Should have no singles.")
    treated = which(z==1)
    control.females = intersect(which(z==0), types.females(true.types))
    treated.m.with.control.f = intersect(which(z==1),
                                         types.male.match(true.types, control.females))
    obs.t[treated.m.with.control.f] <- 0
  }
  return(obs.t)
}

type.colors <- function(types, z=rep(NA, length(types))) {
  cols = list(f="red", m="cyan", s="gray")
  light.cols = list(f="pink", m="lightcyan", s="lightgray")
  f = types.females(types)
  m = types.males(types)
  s = types.singles(types)
  types.str = sapply(1:length(types), function(i) {
    if(i %in% m) return("m")
    else if(i %in% f) return("f")
    else return("s")
  })
  
  node.cols = sapply(types.str, function(s) cols[[s]])
  
  if(sum(!is.na(z)) > 0) {
    control.nodes = which(z==0)
    new.cols = sapply(types.str[control.nodes], function(s) light.cols[[s]])
    node.cols[control.nodes] = new.cols
  }
  return(node.cols)
}

plot.types <- function(pop, plot.obs=F) {
  require(igraph)
  adjlist = list()
  types = population.types.com(pop)
  units = population.all.units(pop)
  z = population.treatment(pop)
  
  if(plot.obs) types = population.types.obs(pop)
  
  f = types.females(types)
  m = types.males(types)
  for(i in 1:length(types)) {
    adjlist[[i]] <- rep(1, 0)
    if(i %in% f) adjlist[[i]] <- c(types.male.match(types, i))
  }
  g = graph.adjlist(adjlist, mode = c("out"), duplicate = TRUE)
  
  if(plot.obs) {
    cols = type.colors(types, z=z)
    V(g)$color = cols
    V(g)$name = sapply(units, function(i) sprintf("%d", i, z[i]))
    plot(g, vertex.size=250/length(types), main="Observed")  
  } else {
    V(g)$n
    cols = type.colors(types, z=NA)
    V(g)$color = cols
    plot(g, vertex.size=250/length(types), main="Complete")
  }
  
  
}

CHECK_types <- function(types) {
  # Checks whether the values in the vector are valid.
  CHECK_MEMBER(types[!is.na(types)], 0:length(types), msg="Correct values")
  
  # Checks whether this type object is self-consistent.
  m = types.males(types)
  f = types.females(types)
  s = types.singles(types)
  # f/m/s sets should be disjoint
  CHECK_DISJOINT(m, f)
  CHECK_DISJOINT(m, s)
  CHECK_DISJOINT(f, s)
  # same number of males and females.
  CHECK_EQ(length(m), length(f))
  
  # male indicator is 0
  CHECK_SETEQ(types[m], c(0))
  # check that matches(all f) = all m
  CHECK_SETEQ(types.male.match(types, f), m)
  # check match(match(f)) = f
  CHECK_identical(types.female.match(types,  types.male.match(types, f)), f)
}

CHECK_consistent_types <- function(types.com, types.obs, z, no.singles) {
  valid.types.obs = observed.types(types.com, z, no.singles)
  CHECK_TRUE(identical(types.obs, valid.types.obs))
}

