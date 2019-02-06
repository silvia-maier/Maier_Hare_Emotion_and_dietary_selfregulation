# Jags-MultivariateNormal.R 
# John Kruschke, November 2015 - June 2017.
# For further info, see:
# Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
# A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

bayesRankCor <-function(y) {
  
  # takes as input a two-column vector, y
  # returns a vector of MCMC samples for the correlation
  
  #--------------------------------------------------------------------------
  # begin 
  #--------------------------------------------------------------------------
  
  # transform y into ranks
  y = apply(y,2,rank) #rank each column to compute rank-based correlations
  
  # Load some functions used below:
  source("DBDA2E-utilities.R") # path must be correct.
  # Install the ellipse package if not already:
  want = c("ellipse")
  have = want %in% rownames(installed.packages())
  if ( any(!have) ) { install.packages( want[!have] ) }
  
  
  # Standardize the data:
  sdOrig = apply(y,2,sd)
  meanOrig = apply(y,2,mean)
  zy = apply(y,2,function(yVec){(yVec-mean(yVec))/sd(yVec)})
  # Assemble data for sending to JAGS:
  dataList = list(
    zy = zy ,
    Ntotal =  nrow(zy) ,
    Nvar = ncol(zy) ,
    # Include original data info for transforming to original scale:
    sdOrig = sdOrig ,
    meanOrig = meanOrig ,
    # For wishart (dwish) prior on inverse covariance matrix:
    zRscal = ncol(zy) ,  # for dwish prior
    zRmat = diag(x=1,nrow=ncol(zy))  # Rmat = diag(apply(y,2,var))
  )
  
  # Define the model:
  modelString = "
model {
for ( i in 1:Ntotal ) {
zy[i,1:Nvar] ~ dmnorm( zMu[1:Nvar] , zInvCovMat[1:Nvar,1:Nvar] ) 
}
for ( varIdx in 1:Nvar ) { zMu[varIdx] ~ dnorm( 0 , 1/2^2 ) }
zInvCovMat ~ dwish( zRmat[1:Nvar,1:Nvar] , zRscal )
# Convert invCovMat to sd and correlation:
zCovMat <- inverse( zInvCovMat )
for ( varIdx in 1:Nvar ) { zSigma[varIdx] <- sqrt(zCovMat[varIdx,varIdx]) }
for ( varIdx1 in 1:Nvar ) { for ( varIdx2 in 1:Nvar ) {
zRho[varIdx1,varIdx2] <- ( zCovMat[varIdx1,varIdx2] 
/ (zSigma[varIdx1]*zSigma[varIdx2]) )
} }
# Convert to original scale:
for ( varIdx in 1:Nvar ) { 
sigma[varIdx] <- zSigma[varIdx] * sdOrig[varIdx] 
mu[varIdx] <- zMu[varIdx] * sdOrig[varIdx] + meanOrig[varIdx]
}
for ( varIdx1 in 1:Nvar ) { for ( varIdx2 in 1:Nvar ) {
rho[varIdx1,varIdx2] <- zRho[varIdx1,varIdx2]
} }
}
" # close quote for modelString
writeLines( modelString , con="Jags-MultivariateNormal-model.txt" )

# Run the chains:
nChain = 3
nAdapt = 500
nBurnIn = 500
nThin = 10
nStepToSave = 20000
require(rjags)
jagsModel = jags.model( file="Jags-MultivariateNormal-model.txt" , 
                        data=dataList , n.chains=nChain , n.adapt=nAdapt )
update( jagsModel , n.iter=nBurnIn )
codaSamples = coda.samples( jagsModel , 
                            variable.names=c("mu","sigma","rho") ,
                            n.iter=nStepToSave/nChain*nThin , thin=nThin )

# Convergence diagnostics:
parameterNames = varnames(codaSamples) # get all parameter names
# for ( parName in parameterNames ) {
#   diagMCMC( codaObject=codaSamples , parName=parName ) 
# }

# Examine the posterior distribution:
mcmcMat = as.matrix(codaSamples)

#rm JAGS model file
file.remove("Jags-MultivariateNormal-model.txt")

return(mcmcMat[,4])
}
#---------------------------------------------------------------------------

# how to read out the output:
# mean_rho = mean(mcmcMat[,4])
# HDI = HDIofMCMC(mcmcMat[,4])
# post_prob = mean(mcmcMat[,4] > 0) #Test greater than 0 if we expect it to be positive


