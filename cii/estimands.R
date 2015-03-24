# getSharedNeighbors(x,v)            -- #shared neighbors
# getNodeResponse(x,v)               -- Returns the Yi= response of node i
# getDeltaEstimand(x)                -- Returns δk   The main estimand we define
# getDeltaEstimandNoRandomization(x) -- Returns δk   The main estimand we define
# getDeltaEstimandMonteCarlo(x)      -- Returns a Monte-Carlo approximation of δk
###  RET set of neighbors which are shared with other nodes in Vk


getSharedNeighbors <- function(x, v) {
    shared.set <- c();
    Nv <- get.neighbors(x,v)
    if(! (v%in% x$Vk) )
        stop("This function is supposed to be called for members of Vk")
    
    
    for(u in Nv) {
        if(u==v) next;
        Nu <- get.neighbors(x, u);
        m <- match(Nu, x$Vk); ##     how many of the neighbors of u are in Vk ?
        ###    v is already in that set
        if(length(which(!is.na(m)))>=2) {
            # print(m)
            shared.set[length(shared.set)+1] <- u;
        }
    }
    return(shared.set);
}

##  Returns Yi = response of node v
##   Measuring is possible if the node is either in T or in C
getNodeResponse  <- function(x, v) {
    return(x$response.fn(x,v))
}

###  Take S= sum of treatment of neighbors
##  Return NA if  S==NA or S*(S-k) !=0 or Zi==NA or Zi==1 or i not in Vk (already covered)
canMeasureResponse <- function(x,v) {
    
    Nv <- get.neighbors(x,v);          ##   Neighbors of v
    subZ <- getNodesTreatment(x,Nv);  ##   the assignment of sub-vector
    zi <- getNodeTreatment(x,v);
    sum_z <- sum(subZ);
    
    if(is.na(sum_z) || sum_z * (sum_z-x$k)!=0  || is.na(zi) || zi==1 ||
           !  (v %in% x$Vk) ) 
        return(FALSE)
    
    return(TRUE)
}

###
###    Use Monte-Carlo   to get an estimate
###
getDeltaEstimandMonteCarlo <- function(y, trials=500, verbose=T) {
    
    if(sum(is.na(y$Z))!= length(y$Z)) 
        stop("Nodes have been assigned. Cannot compute delta estimand!");
    if(length(y$Vk)==0) {
        warning("No nodes in Vk. Quitting")
        return(NA)
    }
    
    Y1 <- c();  ## samples in treatment
    Y0 <- c();  ## samples in control
    
    # flip some coins
    coins <- rbinom(trials, size=1, prob=0.5)
    
    # samples nodes from Vk with replacement
    if(length(y$Vk)==1) 
        sampled.nodes <-  rep(y$Vk[1], trials)
    else 
        sampled.nodes <-  sample(y$Vk, trials, replace=T)
    
    pb <- NA
    if(verbose) {
        print("Monte Carlo delta estimate. Please wait..")
        flush.console()
        pb <- txtProgressBar(title="Monte Carlo delta estimate",char='*',   style=3)
    }
    
    for(i in 1:trials) {
        x <- y
        v <- sampled.nodes[i]; ### sample a node
        
        if(verbose) setTxtProgressBar(pb=pb, value=i/trials)
        
        if(coins[i]==1) {
            x <- treat.node(x, v);
            yi <- getNodeResponse(x,v)$val;
            Y1[length(Y1)+1] <- yi
            #print(sprintf("treat=%.3f", yi))
        } else {
            x <- control.node(x,v);
            yi <- getNodeResponse(x,v)$val;
            Y0[length(Y0)+1] <- yi
            # print(sprintf("control=%.3f", yi))
        }
        
    }
    if(verbose) close(pb)
    if(length(Y0) * length(Y1)==0) 
        stop("Monte Carlo failed.")
    #hist(Y1-Y0)
    #summary(Y1-Y0)
    #readline("done")
    return(mean(Y1)-mean(Y0));
}

# Computes the δκ  by brute force.
getDeltaEstimand <-function(x) {
    if(sum(is.na(x$Z))!= length(x$Z)) 
        stop("Nodes have been assigned. Cannot compute delta estimand!");
    
    delta <- 0;
    vk.size <- length(x$Vk);
    
    for(v in x$Vk) {
        Nv <- get.neighbors(x,v);
        Mv <- assign.complete.rand.matrix(length(Nv), x$k);
        rho.i <- 1/nrow(Mv);
        for(i in 1:nrow(Mv)) {
            y <- x
            y <- setNodesPointTreatment(y,Vs=Nv,Zs=Mv[i,]);
            y <- setNodePointTreatment(y,v,0);
            y$TGroup[1] <- v;
            val <- getNodeResponse(y,v);
            
            if(val$code!=1)
                stop("This should be a node in active treatment!")
            
            delta <- delta + (1/vk.size) * rho.i * val$val;
        }
        
        y2 <- x
        y2 <- setNodesPointTreatment(y2,Vs=Nv,Zs=rep(0, length(Nv)));
        y2 <- setNodePointTreatment(y2,v,0);
        y2$CGroup[1] <- v;
        val <- getNodeResponse(y2,v);
        if(val$code!=0)
            stop("This should be a node in  control!")
        
        delta <- delta - (1/vk.size) * val$val
        
        #  print(delta)
        #  flush.console()
    }
    
    return(delta)
}

### This calculates the delta estimator, given that 
##   we have finished assigning treatment.
delta.estimator <- function(x) {
  no.causation = F
  if(length(x$ex$Exposed)+length(x$ex$Control)==0)
    stop("You haven't assign assigned anyone papi! Try again.")
  
  if( length(x$ex$Can.Expose) +length(x$ex$Can.Control) >0 )
    stop("Further randomization is possible.")
  
  if(length(x$ex$Exposed) * length(x$ex$Control)==0)
  {
    warning("One group of T,C is empty. No causation without manipulation!!")
    no.causation=T
  }
  ##  Calls response function.
  delta_1 <- sapply(x$ex$Exposed, function(v) unit.response(x, v));
  delta_0 <-  sapply(x$ex$Control, function(v) unit.response(x, v) );
  
  ##  Returns NA if no causation
  return( mean(delta_1) - mean(delta_0) )
}
