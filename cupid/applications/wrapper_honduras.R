# Panos Toulis
# ptoulis@fas.harvard.edu
#
# Contains wrapper functions to access/slice the Honduras dataset.
#
#
rm(list=ls())
load("honduras.RData")

library(igraph)

kCouponType = "REAL_MVI_TICKET"
# kCouponType = "REAL_CHLORINE_TICKET"
kOutcome = c()
kOutcomeName = c()
if("REAL_MVI_TICKET" == kCouponType) {
  kOutcome = "QuizScore.M.PCA1.nonseed"
  kOutcomeName <- "score.mvi"
} else if("REAL_CHLORINE_TICKET" == kCouponType) {
  kOutcome = "QuizScore.C.PCA1.either"
  kOutcomeName <- "score.cl"
} else {
  stop("Not correct treatment.")
}
# What data to grab.
kCovariates = c("CODE", "Age", "Male", "Persons_in_house")
kNetworkFeatures = c("SENDER", "RECEIVER", "Nt.type")
kRelations = c("Spouse", "Sibling", "Friend")

village.data <- function(village.id) {
  # Returns the units, the social network, and the arrows.
  # 
  # Value: list(U, Y, G, A, seeds)  where
  #   U = units, Y=outcomes, G=network, A=coupon arrows.
  units = subset(surveys.ep2, Aldea.no==village.id, select=kCovariates)
  rownames(units) <- 1:nrow(units)
  Y = subset(surveys.ep2, Aldea.no==village.id, select=kOutcome)
 rownames(Y) <- 1:nrow(Y)
  G = subset(all.nets[[village.id]], Nt.type %in% kRelations, select=kNetworkFeatures)
  A = subset(all.nets[[village.id]], Nt.type == kCouponType & Wave <= 2, select=c(kNetworkFeatures))
  
  seeds = real.all.seeds.MVI.CODE[[village.id]]
  if(kCouponType == "REAL_CHLORINE_TICKET")
    seeds = real.all.seeds.Cloro.CODE[[village.id]]

  names(units) <- c("id", "age", "male", "household")
  names(G) <- c("sender", "receiver", "arrow.type")
  names(A) <- c("sender", "receiver", "arrow.type")
  names(Y) <- kOutcomeName
  return(list(U=units, Y=as.numeric(Y[,1]), G=G, A=A, seeds=seeds))  
}

plot.villageArrows <- function(vData) {
  # Assumes it is part of the village data.
  #
  # green=seed, color=score in the test.  gray=NA (missing value.)
  A = subset(vData$A, select=c("sender", "receiver"))
  u = vData$U$id # unit ids
  m = apply(A, 2, function(col) match(col, u))
  # m = (.., ..)   has the id's in the U vector of unit ids.
  Adj = matrix(0, nrow=length(u), ncol=length(u))
  Adj[m] <- 1
  g = graph.adjacency(Adj, mode = "directed")
  par(mar=rep(0, 4))
  igraph.options(plot.layout=layout.auto)
  scores = as.vector(unlist(vData$Y))
  sort.scores = sort(scores,  na.last = NA)
  n = length(sort.scores)
  cols = rev(heat.colors(n))
  # Given scores, choose the color.
  vertex.colors = sapply(scores, function(x) {
    if(is.na(x)) return("gray")
    return(cols[match(x, sort.scores)])
  })
  vertex.colors[match(vData$seeds, u)] <- "green"
  print("Seeds are")
  print(match(vData$seeds, u))
  plot.igraph(g, vertex.size=8, vertex.color=vertex.colors, edge.arrow.size=0.4)
}





