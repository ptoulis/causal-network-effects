# Copyright (c) 2015
# Panos Toulis, ptoulis@fas.harvard.edu
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
  # Value: list(X, Y, G, A, seeds)  where
  #   X = units (and covariates), Y=outcomes, G=network, A=coupon arrows.
  #
  # X = (ID, age, male, household)
  # Y  
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
  return(list(X=units, Y=as.numeric(Y[,1]), G=G, A=A, seeds=seeds))  
}

sanitize.village.data <- function(vData) {
  # Will get the village data and then make a transformation 
  # to handle it more easily.
  # e.g. E = network of connections, will become a binary sparse matrix.
  #
  out = list()
  id = vData$X$id
  out$X = as.matrix(vData$X[,2:4])
  E = matrix(0, nrow=length(id), ncol=length(id))
  matches_E = apply(vData$A, 2, function(col) match(col, id))
  for(i in 1:nrow(matches_E)) {
    E[matches_E[i,1], matches_E[i, 2]] <- 1
  }
  out$E <- E
  
  # Sanitize G
  G = matrix(0, nrow=length(id), ncol=length(id))
  matches_G = apply(vData$G, 2, function(col) match(col, id))
  for(i in 1:nrow(matches_G)) {
    G[matches_G[i,1], matches_G[i, 2]] <- 1
  }
  out$G <- G
  out$Y = vData$Y
  out$id = id
  out$Z <- rep(0, length(id))
  out$Z[match(vData$seeds, id)] <- 1
  # Deal with missing data.
  warning("Ad-hoc imputation (wrapper_honduras:L78). Please fix.")
  miss = which(is.na(out$Y))
  out$Y[miss] <- 
    sapply(miss, function(i) {
      if(sum(out$E[, i]) > 0) {
        sample(c(0.05, 0.1, 0.2), size=length(1), prob=c(5, 3, 2))
      } else {
        sample(c(.05, .1, .2), size=length(1), prob=c(5, 3, 2))
      }})
  print(sprintf("Missing imputed but pointeed to = %.2f not pointed to=%.3f", 
                mean(out$Y[intersect(miss, which(colSums(out$E) > 0))]),
                mean(out$Y[intersect(miss, which(colSums(out$E) == 0))])))
  out$Y[which(out$Z==1)] <- 2.5
  return(out)
}

plot.villageData <- function(saneData, plot.coupon.only=F, str="") {
  # saneData = sanitize.village.data(vData)
  # Assumes it is part of the village data.
  # Data needs to be sanitized (see sanitize.village.data())
  # green=seed, color=score in the test.  gray=NA (missing value.)
  #
  A = saneData$E
  G = saneData$G
  u = saneData$id
  
  Adj = matrix(0, nrow=length(u), ncol=length(u))
  g = graph.adjacency(Adj, mode = "directed")
  
  par(mar=rep(0, 4))
  # Scores -> defines vertex colors.
  scores = as.vector(unlist(saneData$Y))
  sort.scores = sort(scores,  na.last = NA)
  n = length(sort.scores)
  cols = rev(heat.colors(n))
  # Given scores, choose the color.
  vertex.colors = sapply(scores, function(x) {
    if(is.na(x)) return("gray")
    return(cols[match(x, sort.scores)])
  })
  seeds = which(saneData$Z==1)
  vertex.colors[seeds] <- "green"
  print("Seeds are")
  print(seeds)

  # Add Edges.
  for(i in 1:nrow(saneData$E)) {
    for(j in which(saneData$E[i, ]==1)) {
      g <- add.edges(g, c(i, j))
    }
  }
  for(i in 1:nrow(saneData$G)) {
    for(j in which(saneData$G[i, ]==1)) {
      g <- add.edges(g, c(i, j))
    }
  }
   
  V(g)$color = vertex.colors
  nE = sum(saneData$E)
  nG = sum(saneData$G)
  E(g)$arrow.size = c(rep(.5, nE), rep(0.5, nG))
  E(g)$width = c(rep(1, nE), rep(0.3, nG))

  E(g)$color = c(rep("orange", nE), rep("white", nG))
      
  par(mar=c(0, 0, 0, 0))
  if(plot.coupon.only) {
    # par(mfrow=c(1, 1))
    plot(g, vertex.size=12, layout=layout.auto)
    return()
  }
  par(mfrow=c(2, 2))
  plot(g, vertex.size=12, layout=layout.auto)
  
  
  E(g)$color = c(rep("orange", nE), rep("white", nG))
  plot(g, vertex.size=12, layout=layout.sphere)
  
  E(g)$color = c(rep("white", nE), rep("blue", nG))
  plot(g, vertex.size=12, layout=layout.sphere)
  
  E(g)$color = c(rep("orange", nE), rep("blue", nG))
  plot(g, vertex.size=12, layout=layout.sphere)
  
  
  print("Quick description")
  print("Color of node=Test score. red=good. green=seed. gray=missing., light=bad score")
  print("Blue edge=social connection.  orange edge=coupon passing")
}





