
# load the posterior distributions of the model parameters
load("HondurasResults20150308_SocialNetOnly_PropModel_logistic_noBackProp/CausalRun_MVITestScore_propagationModel_fwdProp_logistic_indVary_smGammas_varyBase_var1_pass3.RData")

# allocate memory
expectedAvgPeerEffect1 = array(0,numVillage)
expectedAvgPeerEffect2 = array(0,numVillage)
expectedAvgPeerEffect3 = array(0,numVillage)
expectedAvgPeerEffectTotal = array(0,numVillage)
expBaselineEffect_Draws.byVillage = list()
expPrimaryEffect_Draws.byVillage = list()
expPeerEffect1_Draws.byVillage = list()
expPeerEffect2_Draws.byVillage = list()
expPeerEffect3_Draws.byVillage = list()
expAvgPeerEffect1_Draws.byVillage = list()
expAvgPeerEffect2_Draws.byVillage = list()
expAvgPeerEffect3_Draws.byVillage = list()
expResponsePi_Draws.byVillage = list()
expResponse_Draws.byVillage = list()
numTreatment = array(0,numVillage)
population = array(0,numVillage)

#note: number of draws in the imputation
numDraws = 1000

# note: Change here to focus on one village.
for (vil in 1:numVillage)
{
  print(paste("Imputing potential outcome for village #", vil))
  
  # note: treatment = vector of assignment (can be changed).
  treatment = treatment.byVillage[[vil]] # or choose a treatment vector other than the one from the experiment, as long as the length of the vector is the population size
  
  G.edges.lambda = G.edges.lambda.byVillage[[vil]]
  coefficientDraws = sample(1:length(mu_Hist.byVillage[[vil]]),numDraws)
  population[vil] = length(treatment.byVillage[[vil]])
  numTreatment[vil] = sum(treatment)
  
  # allocate memory
  muDraws = mu_Hist.byVillage[[vil]][coefficientDraws]
  tauDraws = tau_Hist.byVillage[[vil]][coefficientDraws]
  S1CoefficientDraws = tau_Hist.byVillage[[vil]][coefficientDraws]*gamma1_Hist.byVillage[[vil]][coefficientDraws]
  S2CoefficientDraws = tau_Hist.byVillage[[vil]][coefficientDraws]*gamma1_Hist.byVillage[[vil]][coefficientDraws]*gamma2_Hist.byVillage[[vil]][coefficientDraws]
  S3CoefficientDraws = tau_Hist.byVillage[[vil]][coefficientDraws]*gamma1_Hist.byVillage[[vil]][coefficientDraws]*gamma2_Hist.byVillage[[vil]][coefficientDraws]*gamma3_Hist.byVillage[[vil]][coefficientDraws]
  expPrimaryEffectExponent = matrix(0,population[vil],numDraws)
  expPeerEffect1Exponent = matrix(0,population[vil],numDraws)
  expPeerEffect2Exponent = matrix(0,population[vil],numDraws)
  expPeerEffect3Exponent = matrix(0,population[vil],numDraws)
  expPrimaryEffect_Draws.byVillage[[vil]] = matrix(0,numDraws,population[vil])
  expPeerEffect1_Draws.byVillage[[vil]] = matrix(0,numDraws,population[vil])
  expPeerEffect2_Draws.byVillage[[vil]] = matrix(0,numDraws,population[vil])
  expPeerEffect3_Draws.byVillage[[vil]] = matrix(0,numDraws,population[vil])
  expBaselineEffect_Draws.byVillage[[vil]] = matrix(0,numDraws,population[vil])
  expAvgPeerEffect1_Draws.byVillage[[vil]] = array(0,numDraws)
  expAvgPeerEffect2_Draws.byVillage[[vil]] = array(0,numDraws)
  expAvgPeerEffect3_Draws.byVillage[[vil]] = array(0,numDraws)
  expResponsePi_Draws.byVillage[[vil]] = matrix(0,numDraws,population[vil])
  expResponse_Draws.byVillage[[vil]] = matrix(-1,numDraws,population[vil])
  
  if (forwardPropOnly) # remove edges for each propagation step so there is no back-feeding into the source nodes (i.e. A influences B in step i, B can not influence A back in step i+1)
  {
    G.edges.lambda.S1 = G.edges.lambda
    treatedNodes = which(treatment>0)
    G.edges.lambda.S2 = G.edges.lambda.S1
    S1Nodes = c() 
    for (i in 1:length(treatedNodes))
    {
      receivers = which(G.edges.lambda.S1[treatedNodes[i],]>0) 
      G.edges.lambda.S2[receivers, treatedNodes[i]] = 0 # zero out the back-propagating edges for step 2 propagation
      S1Nodes = union(S1Nodes, receivers) # step 1 receivers
    }
    G.edges.lambda.S3 = G.edges.lambda.S2
    S2Nodes = c()
    for (i in 1:length(S1Nodes))
    {
      receivers = which(G.edges.lambda.S2[S1Nodes[i],]>0) 
      G.edges.lambda.S3[receivers, S1Nodes[i]] = 0 # zero out the back-propagating edges for step 3 propagation
      S2Nodes = union(S2Nodes, receivers) # step 2 receivers
    }
  } else # allow backward propagation so all edges are on at each step
  {
    G.edges.lambda.S1 = G.edges.lambda
    G.edges.lambda.S2 = G.edges.lambda
    G.edges.lambda.S3 = G.edges.lambda
  }
  
  # compute expected peer influence over some draws of the network sufficient statistics
  for (draw in 1:numDraws) 
  {
    S1.lambda <- t(G.edges.lambda.S1)%*%treatment 
    S1 <- rpois(population[vil],S1.lambda)
    S2.lambda <- t(G.edges.lambda.S2)%*%S1 
    S2 <- rpois(population[vil],S2.lambda)
    S3.lambda <- t(G.edges.lambda.S3)%*%S2 
    S3 <- rpois(population[vil],S3.lambda)  
    expPrimaryEffectExponent[,draw] = treatment*tauDraws[draw]
    expPeerEffect1Exponent[,draw] = S1 * S1CoefficientDraws[draw]
    expPeerEffect2Exponent[,draw] = S2 * S2CoefficientDraws[draw]
    expPeerEffect3Exponent[,draw] = S3 * S3CoefficientDraws[draw]
  }
  
  # impute the potential outcomes
  for (i in 1:population[vil])
  {
    expBaselineExponent = muDraws
    if (exists("epsilon_Hist.byVillage"))
    {
      expBaselineExponent = expBaselineExponent + epsilon_Hist.byVillage[[vil]][coefficientDraws,i]
    }
    expBaselineEffect_Draws.byVillage[[vil]][,i] = 1 / (1+exp(-expBaselineExponent))
    expPrimaryEffect_Draws.byVillage[[vil]][,i] = (1 / (1+exp(-(expPrimaryEffectExponent[i,]+expBaselineExponent)))) - expBaselineEffect_Draws.byVillage[[vil]][,i] 
    expPeerEffect1_Draws.byVillage[[vil]][,i] = (1 / (1+exp(-(expPeerEffect1Exponent[i,]+expPrimaryEffectExponent[i,]+expBaselineExponent)))) - expPrimaryEffect_Draws.byVillage[[vil]][,i] - expBaselineEffect_Draws.byVillage[[vil]][,i] 
    expPeerEffect2_Draws.byVillage[[vil]][,i] = (1 / (1+exp(-(expPeerEffect2Exponent[i,]+expPeerEffect1Exponent[i,]+expPrimaryEffectExponent[i,]+expBaselineExponent)))) - expPeerEffect1_Draws.byVillage[[vil]][,i] - expPrimaryEffect_Draws.byVillage[[vil]][,i] - expBaselineEffect_Draws.byVillage[[vil]][,i]
    expPeerEffect3_Draws.byVillage[[vil]][,i] = (1 / (1+exp(-(expPeerEffect3Exponent[i,]+expPeerEffect2Exponent[i,]+expPeerEffect1Exponent[i,]+expPrimaryEffectExponent[i,]+expBaselineExponent)))) - expPeerEffect2_Draws.byVillage[[vil]][,i] - expPeerEffect1_Draws.byVillage[[vil]][,i] - expPrimaryEffect_Draws.byVillage[[vil]][,i] - expBaselineEffect_Draws.byVillage[[vil]][,i]
    # note: Here we impute the outcomes Y_i
    expResponsePi_Draws.byVillage[[vil]][,i] = (1 / (1+exp(-(expPeerEffect3Exponent[i,]+expPeerEffect2Exponent[i,]+expPeerEffect1Exponent[i,]+expPrimaryEffectExponent[i,]+expBaselineExponent))))
    expResponse_Draws.byVillage[[vil]][,i] = rbinom(numDraws, 1, expResponsePi_Draws.byVillage[[vil]][,i])
  }
  # note: waves.
  expAvgPeerEffect1_Draws.byVillage[[vil]] = rowMeans(expPeerEffect1_Draws.byVillage[[vil]])
  expAvgPeerEffect2_Draws.byVillage[[vil]] = rowMeans(expPeerEffect2_Draws.byVillage[[vil]])
  expAvgPeerEffect3_Draws.byVillage[[vil]] = rowMeans(expPeerEffect3_Draws.byVillage[[vil]])
  
  expectedAvgPeerEffect1[vil] = mean(expAvgPeerEffect1_Draws.byVillage[[vil]])
  expectedAvgPeerEffect2[vil] = mean(expAvgPeerEffect2_Draws.byVillage[[vil]])
  expectedAvgPeerEffect3[vil] = mean(expAvgPeerEffect3_Draws.byVillage[[vil]])    
  expectedAvgPeerEffectTotal[vil] = expectedAvgPeerEffect1[vil] + expectedAvgPeerEffect2[vil] + expectedAvgPeerEffect3[vil]
}