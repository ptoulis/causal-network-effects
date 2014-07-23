## Testing for types.R moodule
rm(list=ls())
source("../types.R")
source("../../../r-toolkit/checks.R")

test.sample.types <- function(niters=1) {
  k = sample(5, size=1)
  N = 10^k
  prop = sample(10, size=1)
  x = sample.types(N, singles.pct = prop/10)
  CHECK_notNA(x)
  
  s = which(x==1:N)
  m = which(x==0)
  f = setdiff(1:N, union(m, s))
  CHECK_TRUE(length(s)==prop * 10^(k-1), msg = "Correct % singles.")
  CHECK_TRUE(length(m)==(N - prop * 10^(k-1))/2, msg = "Correct % singles.")
  CHECK_TRUE(length(m)==length(f), msg="same #of of females and males")
}

test.filter.types <- function() {
  # Test the function that filters the types to males/females etc.
  # Units from 1...N
  all = c(4, 5, 0, 0, 0, 3, 7, 8, NA, NA)
  
  N = length(all)
  # print(all)
  CHECK_SETEQ(c(1,2,6), types.females(all), msg="Correct females")
  CHECK_SETEQ(c(4,5,3), types.males(all))
  CHECK_SETEQ(c(4, 3), types.male.match(all, females = c(1, 6)))
  CHECK_SETEQ(c(7, 8), types.singles(all))
  print(sprintf("Checked population of %d Units: [OK]", N))
  
  # Check match functions.
  N  = 4 * rpois(1, lambda=8)
  all = sample.types(N, singles.pct = 0.3)
  f = types.females(all)
  m.match = types.male.match(all, f)
  # Checks whether   all males == match(all females)
  CHECK_SETEQ(m.match, types.males(all))
  f.sample = sample(f, size=length(f)/2, replace=F)
  # Checks the identity match(match(females)) = females
  CHECK_SETEQ(f.sample, types.female.match(all, types.male.match(all, f.sample)))
  print(sprintf("Checked matches consistency on  %d Units: [OK]", N))
}

test.observed.types <- function() {
  x = c(1, 2, 0, 0, 3, 0, 4, 6)
  # (1, 2) singles
  # 5 -> 3, 7->4, 8->6
  CHECK_types(x)
  # treat f(7), s(2), m(3)
  z = c(0, 1, 1, 0, 0, 0, 1, 0)
  CHECK_EXCEPTION(observed.types(x, z, no.singles = T), msg="there are singles")
  t.obs = observed.types(x, z, no.singles=F)
  t.obs.theoretical = c(NA, NA, NA, 0, NA, NA, 4, NA)
  CHECK_identical(t.obs, t.obs.theoretical)
}

test.checks <- function() {
  types = c(0, 0, 1, 3)
  CHECK_EXCEPTION(CHECK_types(types))
  types = c(0, 0, 1, NA, 2, 100)
  # error: 2 males but one female
  CHECK_EXCEPTION(CHECK_types(types))
}