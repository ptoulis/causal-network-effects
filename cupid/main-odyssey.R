source("inference.R")
# First argument determines how many iterations(imputations) to perform for ech p-value
args <- commandArgs(trailingOnly = TRUE)
true.effect = as.numeric(args[1])
niters = 10
if(length(args) == 2) {
 niters = as.numeric(args[2])
}

data = generate.data(Nt=200, Nc=300, cupid.effect=true.effect);
p = NA
if(length(args)==0) {
  p = run.most.powerful.test(data, niters=10);
} else {
  p = run.most.powerful.test(data, niters=niters)
}
print(p)
