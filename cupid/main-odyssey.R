source("inference.R")
# First argument determines how many iterations(imputations) to perform for ech p-value
args <- commandArgs(trailingOnly = TRUE)
data = generate.data(Nt=200, Nc=300, cupid.effect=0);
p = NA
if(length(args)==0) {
  p = run.most.powerful.test(data, niters=10);
} else {
  p = run.most.powerful.test(data, niters=as.numeric(args[1]))
}
print(p)
