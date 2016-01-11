# Copyright (c) 2015
# Panos Toulis, ptoulis@fas.harvard.edu
# 
#
# Split-plot design with interference within plots but not between plots.
#
# Notation:
# plot  = whole unit of the split-plot design
# unit = subunit in the design (has unique ID)
#
# X = (x1, x2, ...) = covariates
# G = (snd, rec, val) = social network (SNA edge list)
# Y = observed data (numerical col. vector)
# Z = (Z_plot, Z_unit) = assignment vector (factors usually)
# Mz = (snd, rec, val) = observed interference (SNA edge list)
# S = unit x unit binary matrix of possible I/F  (not obs.)
#
# mechanism Q: Mz = Q(S, Z) -- all vars within-plot (wp)
#
# Abbreviation:
#  wp = within-plot
#  if or i/f = interference
#  pid = plot id 
#  
rm(list=ls())
library(network)
library(sna)
library(mvtnorm)
source("../../../r-toolkit/checks.R")

## Manipulations of the dataset
get_plot_X <- function(pid, dataset) {
  return(as.matrix(dataset[["X"]][[pid]]))
}
get_plot_G <- function(pid, dataset) {
  return(as.matrix(dataset[["G"]][[pid]]))
}
get_plot_Mz <- function(pid, dataset) {
  return(as.matrix(dataset[["Mz"]][[pid]]))
}
get_plot_Z <- function(pid, dataset) {
  return(dataset$Z[[pid]][,2])
}
get_plot_Y <- function(pid, dataset) {
  return(dataset$Y[[pid]])
}


plot_S_wp <- function(pid, dataset) {
  par(mar=c(1, 1, 1, 1))
  gplot(dataset$S[[pid]], displayisolates = T, displaylabels = T, mode="circle")
}
plot_G_wp <- function(pid, dataset) {
  par(mar=c(1, 1, 1, 1))
  gplot(dataset$G[[pid]], displayisolates = T, displaylabels = T, mode="circle")
}
plot_Mz_wp <- function(pid, dataset) {
  par(mar=c(1, 1, 1, 1))
  Z = get_plot_Z(pid, dataset)
  col_v = sapply(Z, function(i) ifelse(i==1, "green", "red"))
  gplot(dataset$Mz[[pid]], displayisolates = T, displaylabels = T, vertex.col = col_v, mode="circle")
}


example_if_mechanism <- function(S, Z) {
  # Example interference mechanism.
  #
  # Given source (S) and assignment (Z) compute the realized i/f.
  # In this example, when treated one unit makes all possible i/fs.
  Mz = S
  treated_units = which(Z==1)
  Mz[-c(treated_units), ] <- 0
  for(u in treated_units) {
    out_edges = Mz[u, ]
    if(sum(out_edges) > 0) {
      active = which(out_edges==1)
      choose_units = NA
      if(length(active) < 3) { 
        choose_units = active
      } else {
        choose_units = sample(active, size=3, replace=F)
      }
      Mz[u, -c(choose_units)] <- 0
    }
  }
  return(Mz)
}

synthetic_dataset <- function() {
  df = list()
  P = 2 # 2 plots
  N = 30 # 10 units per plot
  all_units = expand.grid(unit=seq(1, N), plot=seq(1, P))
  d0 = data.frame(plot=all_units$plot, unit=all_units$unit)
  ## 1. Sample covariates
  nx = 3 # no. of covariates
  df = list(X=list(), G=list(), S=list(), Z=list(), Mz=list(), Y=list())
  
  for(plot_id in 1:P) {
    # 1. Sample subunit covariates for plot_id
    X = matrix(sample(c(-1, 0, 1), size=N * nx, replace=T, prob=c(2, 5, 1)), ncol=nx)
    colnames(X) = c("x1", "x2", "x3")
    df[["X"]][[plot_id]] = X
    
    ## 2. Sample network for plot_id
    G = rgnm(1, N, as.integer(0.5 * N), mode="graph")
    df[["G"]][[plot_id]] = as.edgelist.sna(G, )
    
    ## 3. Sample assignment Z
    # two factors (A, B); A=plot treatment, B=unit treatment
    # TODO(ptoulis): Implement multilevel factors.
    Z_plot = ifelse(plot_id%%2==0, rep(0, N), rep(1, N))
    Z_unit = sample(c(rep(0, 0.8*N), rep(1, 0.2 * N)))
    M = matrix(0, nrow=N, ncol=2)
    colnames(M) <- c("Z_plot", "Z_unit")
    M[,1] <- Z_plot
    M[, 2] <- Z_unit
    df[["Z"]][[plot_id]] <- M
    #
    # 4. Sample S - source of interference (possible IF)
    #
    S <- matrix(0, nrow=0, ncol=4)
    x1 = X[, 1] # only the first column matters.
    dist = x1 %*% t(x1)
    diag(dist) <- 0
    dist[upper.tri(dist)] <- 0
    dist[dist < 0] <- 0
    S <- dist
    A = as.edgelist.sna(S)
    df[["S"]][[plot_id]] <- A
    # 4b. Generate Mz according to mechanism
    #
    df[["Mz"]][[plot_id]] <- as.edgelist.sna(example_if_mechanism(S, Z_unit))
    
    # 5. Sample Y
    Y = matrix(0, ncol=1, nrow=N)
    for(u in 1:N) {
      s = df[["Mz"]][[plot_id]][, 2]
      unit_if = sum(s==u)
      Y[u] = rnorm(1, mean=5 * unit_if, sd=0.1)
    }

    df[["Y"]][[plot_id]] <- Y
  }
  return(df)
}

CHECK_dataset <- function(dataset) {
  ## Makes a quick check whether the data are valid.
  CHECK_MEMBER(names(dataset), c("X", "G", "Y", "Z", "Mz"))
}

test_effect_wp <- function(pid,dataset) {
  # Null hypothesis of no effect within a plot (wp)
  # Ho: Y_i(Z) = Y_i(Z') if levels are a or b in Z, Z'
  #
  # TODO(ptoulis): When we implement multi-level factors we can revise this one 
  # to accept factor levels as argument; the levels should have the same response under Ho.
  # 
  data_0 = dataset # save the current dataset
  T_values = c() ## test statistic values
  Tobs = test_statistic_wp(pid, dataset)
  for(i in 1:100) {
    S_imp = impute_sourceIF_wp(pid, dataset)
    Z_new = resample_Z_Ho(S_imp, dataset)
    #
    dataset = update_dataset_wp(pid, level_a, level_b, Z_new, Smis, dataset)
    #
    T_values = c(test_statistic_wp(pid, dataset))
  }
}

