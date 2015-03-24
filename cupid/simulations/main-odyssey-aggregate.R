
read.pvalues <- function() {
files = list.files(path="odyssey", full.names=T, pattern=".out$")
pvalues = c()
# print(sprintf("I have %d files to read", length(files)));

for(f in files) {
  x = readLines(f)[1];
  m = gregexpr(pattern="[\\d\\.]+$", perl=T, text=x)
  # print(sprintf("Input=%s", x))
  value = as.numeric(regmatches(x, m)[[1]]);
  ## print(sprintf("File = %s  Input = %s  value=%.3f", f, x, value))
  pvalues = c(pvalues, value)
}
  return(pvalues);
}

cleanup.pvalues <- function() {
  out.files = list.files(path="odyssey", full.names=T, pattern=".out$");
  for(f in out.files) {
  	unlink(f)
  }
  err.files = list.files(path="odyssey", full.names=T, pattern=".err$");
  for(f in err.files) {
    unlink(f)
 }
}

