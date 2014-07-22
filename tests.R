# Panos Toulis, ptoulis@fas.harvard.edu
# Unit-testing routines.. 


throw.error <- function(str.many) {
  for (str in str.many) 
    print(sprintf("[TEST-FAIL] %s", str))
  stop("Quitting.")
}

are.equal.sets = function(x,y) {
  if(length(x) != length(y)) return(F)
  return(length(setdiff(x,y))==0)
}

is.subset = function(bigger, smaller) {
  xandy = intersect(bigger, smaller)
  return(are.equal.sets(xandy, smaller))
}

check.has.key <- function(x, key, str="n/a") {
  if(! key %in% names(x) )
    throw.error(c(sprintf("Object x does not have key %s", key), str))
  return(T)
}

check.sets.disjoint <- function(x, y, str="n/a") {
  if(length(intersect(x,y))>0)
    throw.error(c(sprintf("[TEST FAIL]...Sets x,y not disjoint. %s"), str))
}

## Testing set-equality
check.sets.eq <- function(x, y, str="n/a") {
  x = unique(x)
  y = unique(y)
  if(length(x) != length(y))
    throw.error(c(sprintf("Lists not equal length."), str))
  if(length(setdiff(x,y))>0)
    throw.error(c(sprintf("Lists of equal length but difference is not equal."), str))
}

check.lists.eq <- function(x,y, str="n/a", tol=0) {
  s = sum(which(abs(x-y)>tol))
  if(s > 0) 
    throw.error(c(sprintf("Lists x,y not equal at tol=%.3f", tol), str))
}

check.lists.geq <- function(x, y, str="n/a") {
  s = sum(which(x-y<0))
  if(s > 0)
    throw.error(c(sprintf("Does not hold that x >= y"), str))
}

check.strings.eq <- function(x, y, str="n/a") {
  if(x != y) 
    throw.error(c(sprintf("Strings %s != %s", x, y), str))
}

check.is.subset <- function(smaller, bigger, str="n/a)") {
  smaller = unique(smaller)
  bigger = unique(bigger)
  if(!is.subset(bigger=bigger, smaller=smaller))
    throw.error(c(sprintf("x is not subset of y"), str))
}

check.mean.eq <- function(x, mu0, str="n/a") {
  sd = bootstrap.mean(x)
  mu = mean(x)
  ci = c(mu - 2*sd, mu + 2*sd)
  
  if(mu0 < ci[1] || mu0 > ci[2])
    stop(sprintf("[TEST FAIL] Hypothesis test mu=mu0  mu0=%.3f  CI=[%.3f, %.3f]: %s", 
                 mu0, ci[1], ci[2], str))
  if(sd> mu/2)
    warning("SE is probably too high. Try increasing the sample size?")
}

check.true <- function(x, str="n/a") {
  if(!x) 
   throw.error(c("x != TRUE", str))
}

check.object <- function(x, class.name, attrs, str="n/a") {
  check.strings.eq(class(x), class.name, str=sprintf("%s/Class name check.", str))
  check.is.subset(names(x), attrs, str=sprintf("%s/Object attributes check", str))
}

check.exception <- function(fn, args, str="n/a") {
  error.occured <<- F
  tryCatch(do.call(fn, args=args), 
           error=function(error) { error.occured <<- T },
           finally={ if (error.occured) {
                        print(sprintf("[OK] Exception was produced. More:%s", str))
                    } else {
                       stop(sprintf("Exception was not produced. More: %s", str))
                     }
                   })
}

# Project-specific tests.
test.cii <- function() {
  # This constructs an ad-hoc cii object and checks all the functions.
  A <- matrix(0, nrow=10, ncol=10)
  A[c(1,3,5,7, 10), 4] <- 1
  A[c(2,3,6), 5] <- 1
  A[c(1,7,8), 6] <- 1
  A[c(1,3,4,6,8), 10] <- 1
  A[c(1,2,3,7), 8] <- 1
  g <- graph.adjacency(A, mode="directed")
  args = list(G=g, Y=function(i, Zni) {i + mean(Zni, na.rm=T)}, k=3)
  cii <- init.cii(args)
  plot.cii(cii)
  # neighbors of 4=(1,3,4,5,7,10)
  check.sets.eq(get.unit.neighbors(cii, 4), c(1,3,5,7,10))
  # neighbors of 8
  check.sets.eq(get.unit.neighbors(cii, 8), c(1,2,3,7))
  # Vk set = {4,5,6,8,10} (those with >=3 incoming links). Recall cii$k=3
  check.sets.eq(get.cii.Vk(cii), c(4,5,6,8,10))
  # No unit is exposed yet.
  check.true(!is.unit.exposed(cii, 2))
  # all initial Z=NA
  check.true(sum(is.na(cii$experiment$Z))==10, str="initial Z=all NAs")
  # In the beginning, can-expose=can-control
  check.sets.eq(cii$experiment$can.expose, cii$experiment$can.control)
  # Units in cii are represented by {1..10}
  check.sets.eq(get.cii.units(cii), 1:10, str="Units from 1..10")
  # Randomization test.
  check.lists.eq(k.ones.of.n(0, 10), rep(0, 10), "Randomiz test")
  check.exception(k.ones.of.n, args=list(k=5, n=2),
                  str="k <=n when sampling k out of n")
  
  # Start a randomization.
  cii <- cii.control.unit(cii, 6)
  plot.cii(cii)
  # Unit 6 in k-control -> units 1,7,8 have Zj=0
  check.sets.eq(get.units.Z(cii, c(1,7,8)), rep(0, 3), str="neighbors of 6 are in control")
  check.true(get.unit.Zi(cii, 6)==0, str="Z6=0 because in k-control")
  # Expose 4
  cii <- cii.expose.unit(cii, 4)
  plot.cii(cii)
  # Exactly k=3 neighbors are set Zj=1
  check.true(sum(get.unit.Znei(cii, 4))==3, str="Should treat exacty 3 neighbors")
  check.sets.eq(cii$experiment$exposed, c(4), str="Only 4 is exposed.")
  # Only one control so far
  check.true(length(cii$experiment$control)==1, str="One control")
  # check out the response definition above. Y=id+mean(neighborhood)
  check.true(get.unit.response(cii, 4, Znei=get.unit.Znei(cii, 4)) == (4+0.6),
             str="Response of unit4 should be 4.6")
  # cannot control unti 10 because neighborhood of 4 interferes
  # This function should throw an exception.
  check.exception(cii.control.unit, args=list(cii=cii, unit=10), str="Cannot control unit 10")
  check.true(length(cii$experiment$can.expose)==0, str="No more to expose")
  check.true(length(cii$experiment$can.control)==0, str="No more to control")
}





