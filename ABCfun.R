#ABC Machinery used in Estimating Heterogeneity in Mosquito-Bird Contact Using Approximate Bayesian Computation
#Author: Jordyn Walton


#===============================================================================
#Distance function between two Occupancy spectrums, for use with Occupancy spectrums formatted in matrix
#see hetmosqfuns.R
SpectrumDifference <-function(S_t = S_True,S_sample){
  #get rid of s_0
  s1<-S_t[,-1]
  s2<-S_sample[-1]
  
  #Distance metric is the euclidean distance
  sqrt(sum((s1-s2)^2))
  
}


#===============================================================================
#Weighted Distance function between two Occupancy spectrums, for use with Occupancy spectrums formatted in matrix
#not used in final draft
SpectrumWDifference <-function(S_t = S_True,S_sample){
  #get rid of s_0
  s1<-S_t[,-1]
  s2<-S_sample[-1]
  w <- 1.1^(1:length(s2))
  #Distance metric is a weighted euclidean distance
  sqrt(sum(w*(s1-s2)^2))
  
}

#===============================================================================
#Priors
#===============================================================================

#===============================================================================
# Generates Negative Binomial prior for B
GenPriorBnBinom <-function(NoOfValues,size=6,p=0.1){
  B<-rnbinom(NoOfValues,size=size,p=p)+1
}

#===============================================================================
#Generates matrix. Row 1 is negative Binomial prior for B, Row 2 is Gamma prior for theta
GenPriorBTheta<-function(NoOfValues,Bsize=20,Bp=0.4,TShape = 5,TScale = 10){
  B <- rnbinom(NoOfValues,size=Bsize,p=Bp)+Bsize
  mu <- 1/B
  theta <-rgamma(NoOfValues,shape=TShape,scale = TScale)
  All_Pos_Parms<-cbind(B,theta)
  
  All_Pos_Parms

}




#===============================================================================
#Computes Euclidean distance from simulated Occupancies (inputted in matrix) and the true

AllDistance <- function(S_t,AllOccupSpec){
  Distance<- vector(mode="numeric",length=0)
  NoOfValues <-nrow(AllOccupSpec)
  for(trial in 1:NoOfValues){
    Distance<- c(Distance,SpectrumDifference(S_t=S_t,S_sample = AllOccupSpec[trial,]))
  }
  
  
  Distance
}

#===============================================================================
#Computes weighted euclidean distance from from simulated Occupancies (inputted in matrix) and the true
#unused in final draft

AllWDistance <- function(S_t,AllOccupSpec){
  Distance<- vector(mode="numeric",length=0)
  NoOfValues <-nrow(AllOccupSpec)
  for(trial in 1:NoOfValues){
    Distance<- c(Distance,SpectrumWDifference(S_t=S_t,S_sample = AllOccupSpec[trial,]))
  }
  
  
  Distance
}


#===============================================================================
#Returns Posterior after ABC rejection procedure
GetPosterior <-function(prior,distances,cutoff){
  
  acceptance = sum(distances <= cutoff)
  
  #If all distances are above cutoff, Posterior is NAs
  if(acceptance == 0){
    if(is.vector(prior)) Posterior=c(B=NA)
    
    if(is.matrix(prior)) Posterior = cbind(B=c(NA,NA),theta=c(NA,NA))
  }else{ #If some distances are below cutoff
    
    #Which indices are less than cutoff
    accepted_Parms_ind <- which(distances<= cutoff)
    
    #If estimating one value (B)
    if(is.vector(prior)){
      Posterior =prior[accepted_Parms_ind]
    }
    if(is.matrix(prior)){ #If estimating both B and theta
      Posterior = prior[accepted_Parms_ind,]
    }
    
    #Output: Posterior (matrix or a vector)
    Posterior
  }
  
}
#===============================================================================
#Mirror function of GetPosterior, but returns all of the rejected values
#Unused in final draft
GetRejections <-function(prior,distances,cutoff){
  
  rejections = sum(distances > cutoff)
  
  #if nothing is rejected, return NAs
  if(rejections == 0){
    if(is.vector(prior)) Rejected=c(B=NA)
    
    if(is.matrix(prior)) Rejected = cbind(B=c(NA,NA),theta=c(NA,NA))
  }else{  #if some are rejected
    
    #Which indices are rejected
    rejected_Parms_ind <- which(distances > cutoff)
    
    #if estimating one value (B)
    if(is.vector(prior)){
      Rejected =prior[rejected_Parms_ind]
    } 
    #if estimating both B and theta
    if(is.matrix(prior)){
      Rejected = prior[rejected_Parms_ind,]
    }
    
    #Output: the matrix or vector of rejected
    Rejected
  }
  
}



#===============================================================================
# Simulate Occupancy fuction
# 2. Simulate data from samples from prior 
# Input:
#   - No of Mosquitoes
#   - Prior
#     - B is in column 1, theta in 2 
#   - Is ABC assuming homogeneous probabilities?
#     - if false, Prior must have entries for B and theta
#     - if true, Prior only needs entries for B
#===============================================================================

SimulateOccupancy <-function(M,Prior,HomogAssump = TRUE){
  
  #case: Heterogeneous assumption
  if (!HomogAssump){
    stopifnot(is.matrix(Prior))
    NoOfParmsEst = ncol(Prior)
    NoOfValues = nrow(Prior)
    B = Prior[,1]; theta = Prior[,2]
  }
  #case: Homogeneous
  if (is.vector(Prior)){
    NoOfValues = length(Prior)
    B = Prior
  }
  
  #Setup Matrix containing all occupancy spectrums
  AllOccupSpecMat<- matrix(data=0,NoOfValues,M+1,dimnames = list(Trials = 1:NoOfValues,Occupancy = 0:M))
  
  #Iterate through Prior values
  for (trial in (1:NoOfValues)){
    B_trial = B[trial]
    
    #Homogeneous assumption will generate occupancies with B
    if (HomogAssump){
      V_trial <- vfun(B_trial,M)
      SWith0 <- table(V_trial)
      AllOccupSpecMat[trial,names(SWith0)]<-SWith0
    }
    #Assuming Heterogeneous will generate occupancies with B and theta from the prior
    if (!HomogAssump){
      theta_trial = theta[trial]
      
      #Calculate Shapes of beta distribution
      Shape1_trial = theta_trial /B_trial
      Shape2_trial = theta_trial-Shape1_trial
      
      #Generate Probabilities
      prob_trial = rbeta(B_trial,Shape1_trial,Shape2_trial)
      prob_trial = prob_trial/sum(prob_trial)
      
      #Generate Occupancy spectrum
      V_trial <- vHetfun(M,prob_trial)
      SWith0 <- table(V_trial)
      AllOccupSpecMat[trial,names(SWith0)]<-SWith0
    }
  }
  #Output full matrix
  AllOccupSpecMat
}



#===============================================================================
#gets distance vector and returns the 10th percentile value by default,
#constant mean to use to generate the cutoff value
#can be run with the best guess values, aka the prior mean

#Default theta == Homogeneous assumption
genABCcutoff <- function(B,M,level = 0.1,theta_true = 0){
  
  #Case: Homogeneous assumption == theta_true = 0
  if (theta_true==0){
    #An example of a true occupancy
    s_true <- sHetfun2(M,rep(1,B))
    All <- SimulateOccupancy(M,rep(B,1000),HomogAssump = TRUE)
    #Get vector of a thousand Differences
    Dist <- AllDistance(s_true,All)
  }
  #Case: Heterogeneous assumption == theta_true >0
  if(theta_true > 0){
    #Generate a beta distributed probability vector
    mu = 1/B
    Shape1_True = theta_true/B
    Shape2_True = theta_true - Shape1_True
    prob_true <- rbeta(B,Shape1_True,Shape1_True)
    
    #True occupancy
    s_true <-sHetfun2(M,prob_true)
    #the "prior"
    BTheta <- cbind(rep(B,1000),rep(theta_true,1000))
    #Get vector of a thousand differences
    All<- SimulateOccupancy(M,BTheta,HomogAssump = FALSE)
    Dist <- AllDistance(s_true,All)
    
  }
  
  #Output the desired quantile of the distance vector
  quantile(Dist,probs = level)
}



#===============================================================================
#Similar function as the second but takes the Occupancy spectrum as an input, will
#allow for estimating B with an occupancy spectrum generated heterogeneously or
#spatially
EstHomogPosterior2 <- function(B_true,M_fix,S_true,BPrior){
  #true value
  
  #determine cutoff value from constant
  cutoff <- genABCcutoff(B_true,M_fix,level =0.90)
  
  #Generate prior
  
  All <- SimulateOccupancy(M_fix,BPrior,HomogAssump = TRUE)
  
  #get distances from all occupancy spectrums
  Dist <-AllDistance(S_true,All)
  
  #Rejects all values > epsilon
  posterior <- GetPosterior(BPrior,Dist,cutoff)

  #Outputs list containing all information  
  me <- list(
    B_true = B_true,
    S_true = S_true,
    cutoff =cutoff,
    BPrior = BPrior,
    Posterior = posterior,
    Distances =Dist
    
    
  )
  me
}




#===============================================================================
# Estimates posterior of B and theta when heterogeneity is assumed
# - inputs: B and theta (for observed occupancy), M (used to generate all occupancies (both observed and ABC))
#           the Prior samples and the cutoff values

thetaBEstFun <- function(B, theta, M, Prior,ABCcutoff){
  
  #Generate Beta probability vector
  Shape1 = theta / B
  Shape2 = theta - Shape1
  prob_true <- rbeta(B,Shape1,Shape2)
  
  #Generate the "observed" occupancy
  S_true <- sHetfun2(M,prob_true)
  
  #Simulate Occupancy given prior
  AllOccup <- SimulateOccupancy(M,Prior,HomogAssump = FALSE)
  #Calculate distances
  Dist <- AllDistance(S_true,AllOccup)
  #Accept prior values with distance < epsilon (cutoff)
  Posterior <- GetPosterior(Prior,Dist,ABCcutoff)
  
  #avoid dimension issue
  if (is.vector(Posterior)){
    BPosterior <- Posterior[1];ThetaPosterior <- Posterior[2]
  }
  if (is.matrix(Posterior)){
    BPosterior <- Posterior[,1];ThetaPosterior <- Posterior[,2]
  }
  
  ##Calculate Acceptance
  if(sum(is.na(BPosterior))==length(BPosterior)){
    Accepted = 0
  }else{
    Accepted = length(BPosterior)/nrow(Prior)*100
  }
  
  #Output Result: point estimate and CI of B and theta, Acceptance rate of ABC
  B_estimate <- mean(BPosterior)
  B_CI <- quantile(BPosterior,c(0.025,0.975),na.rm=TRUE)
  theta_estimate <- mean(ThetaPosterior)
  theta_CI <- quantile(ThetaPosterior,c(0.025,0.975),na.rm = TRUE)
  
  c(B_est = B_estimate,
    B_CI = B_CI,
    theta_est = theta_estimate,
    theta_CI = theta_CI,
    Acceptance = Accepted
  )
  
}
#===============================================================================
## Coefficient of Variation in a beta distribution

CoVariationBeta <- function(b,theta){
  ((b-1)/(theta + 1))^0.5
}
#==================================================================================================
#Functions for the Cov with Occupancy

#Function: CovEstFun
# - Gives a point estimate and a credible interval of the coefficient of variation
#   in mosquito feeding probabilities by the occupancy spectrum generated within the 
#   function
# Input: The true B and theta of the observed occupancy spectrum, Number of mosquitoes,
#       Values drawn from priors of B and theta (matrix) and a predetermined cutoff 
#       for ABC
# Output: Vector with Point estimate and endpoints of the CI

CovEstFun <-function(B,theta,M,Prior,ABCcutoff){
  stopifnot(is.matrix(Prior))
  NoOfTrials <- nrow(Prior)
  
  #Generate beta probabilties
  mu = 1/B
  Shape1 <-theta*mu
  Shape2 <- theta- Shape1
  prob_true <-rbeta(B,Shape1,Shape2)
  
  #Generate true Occupancy
  S_true <-sHetfun2(M,prob_true)
  
  #Simulate Occupancies from the prior B and thetas
  AllOccup<-SimulateOccupancy(M,Prior,HomogAssump = FALSE)
  
  #Obtain Euclidean distances between true Occupancy and simulated Occupancies
  Dist<-AllDistance(S_true,AllOccup)
  
  #Reject all distances greater than ABC cutoff, Gets prior of joint distribution
  # of B and Theta
  Posterior <-GetPosterior(Prior,Dist,ABCcutoff)
  
  #Accounting for dimensional issue
  if (is.vector(Posterior)){
    BPosterior <- Posterior[1]
    ThetaPosterior <- Posterior[2]
  }
  if (is.matrix(Posterior)){
    BPosterior <- Posterior[,1]
    ThetaPosterior <- Posterior[,2]
  }
  
  #Calculate Coefficient of Variation from B and theta joint posterior
  CovPosterior <- CoVariationBeta(BPosterior,ThetaPosterior)
  
  
  #Output Results: Point estimate, CI  
  Cov_estimate = mean(CovPosterior)
  Cov_CI = quantile(CovPosterior,c(0.025,0.975),na.rm = TRUE)
  
  c(CoVestimate = Cov_estimate
    ,Cov_CI =Cov_CI)
}

#==================================================================================================
#Functions for the Moments


#===============================================================================
##Function that gets the moments of the occupancy spectrum
#Assumes s_occupancy has format given with sHetfun2
getMoments<-function(s_occupancy){
  i = 1:(length(s_occupancy)-1)
  i2 = i^2
  i3 = i^3
  
  # eta_n = \sum (i^n s_i/0 / \sum s_i/0 )
  eta1 = sum(i*s_occupancy[-1]/sum(s_occupancy[-1]))
  eta2 = sum(i2*s_occupancy[-1]/sum(s_occupancy[-1]))
  eta3 = sum(i3*s_occupancy[-1]/sum(s_occupancy[-1]))
  
  #Output vector of Moments
  c(eta1,eta2,eta3)
}

#===============================================================================
#Simulate Occupancy Moments
# Essentially the same as SimulateOccupancy function except with moments

SimulateOccupancyMom <-function(M,Prior,HomogAssump = TRUE){
  
  #Case Heterogeneous assumption (Prior needs to be a matrix)
  if (!HomogAssump){
    stopifnot(is.matrix(Prior))
    NoOfParmsEst = ncol(Prior)
    NoOfValues = nrow(Prior)
    B = Prior[,1]; theta = Prior[,2]
  }
  #Case Homogeneous assumption
  if (is.vector(Prior)){
    NoOfValues = length(Prior)
    B = Prior
  }
  
  #Setup Matrix containing all occupancy spectrum moments
  AllOccupSpecMom<- matrix(data=0,NoOfValues,3,dimnames = list(Trials = 1:NoOfValues,Moments = 1:3))
  
  #Iterate through prior, Simulate occupancy moments
  for (trial in (1:NoOfValues)){
    B_trial = B[trial]
    
    #Homogeneous assumption
    if (HomogAssump){
      #Simulate occupancy
      SWith0<- sHetfun2(M,rep(1,B_trial))
      #Get moments of occupancy simulated
      SMom <-getMoments(SWith0)
      #Enter into matrix
      AllOccupSpecMom[trial,]<-SMom
    }
    #Heterogeneity assumption
    if (!HomogAssump){
      theta_trial = theta[trial]
      
      #Generate beta probability vector
      Shape1_trial = theta_trial /B_trial
      Shape2_trial = theta_trial-Shape1_trial
      prob_trial = rbeta(B_trial,Shape1_trial,Shape2_trial)
      prob_trial = prob_trial/sum(prob_trial)
      
      #Simulate occupancy and get moments
      SWith0<-sHetfun2(M,prob_trial)
      SMom<-getMoments(SWith0)
      
      #Enter into matrix row
      AllOccupSpecMom[trial,]<-SMom
    }
  }
  
  #Output Occupancy moments (in matrix containing 1000 simulations)
  AllOccupSpecMom
}

#===============================================================================
#Calculate Distances of Moment between observed and simulated
AllMomDistance<- function(true_moments,sample_moments){
  Distance <-vector(mode = "numeric",length=0)
  NoOfValues<-nrow(sample_moments)
  #Iterate through moments and collect euclidean distance
  for(trial in 1:NoOfValues){
    Distance<- c(Distance,sqrt(sum((true_moments-sample_moments[trial,])^2)))
  }
  #Output Distance vectors of all moments
  Distance
}

#===============================================================================
#Generate a suggested ABC cutoff for moments
genMomentCutoff <- function(B,M,level = 0.95,theta_true = 0){
  
  #Case: Homogeneous Assumption, theta = 0
  if (theta_true==0){
    #Simulate observed occupancy and get moments
    s_true <- sHetfun2(M,rep(1,B))
    Mom_true <-getMoments(s_true)
    
    #Simulate Occupancy moments with the same values of B and theta
    All <- SimulateOccupancyMom(M,rep(B,1000),HomogAssump = TRUE)
    #Compute distances between true and simulated moments
    Dist <- AllMomDistance(Mom_true,All)
  }
  #Case: Heterogeneity assumption, theta >0
  if(theta_true > 0){
    #Generate beta probability vector
    mu = 1/B
    Shape1_True = theta_true/B
    Shape2_True = theta_true - Shape1_True
    prob_true <- rbeta(B,Shape1_True,Shape1_True)
    
    #Simulate true Occupancy and get moments
    s_true <-sHetfun2(M,prob_true)
    Mom_true <-getMoments(s_true)
    
    #BTheta a "prior"
    BTheta <- cbind(rep(B,1000),rep(theta_true,1000))
    #Simulate occupancy moments with the same values
    All<- SimulateOccupancyMom(M,BTheta,HomogAssump = FALSE)
    #Compute distances between true and simulated
    Dist <- AllMomDistance(Mom_true,All)
    
  }
  #Output desired quantile of Distances
  quantile(Dist,probs = level)
}

#===============================================================================
#Returns The occupancy moments that were not rejected
#Input: The observed moments, the simulated and the ABC cutoff

AcceptedMomOccupancy <- function(True_Moments, AllSimOccupMoments,ABCcutoff){
  
  #Calculate distances
  Dist<-AllMomDistance(True_Moments,AllSimOccupMoments)
  #Collect indices have distances less than < = cutoff
  accepted_Parms_ind <- which(Dist<= ABCcutoff)
  
  #Return Simulated moments that are accepted
  AllSimOccupMoments[accepted_Parms_ind,]
}

#===============================================================================
#Returns The occupancy moments that were rejected
RejectedMomOccupancy <- function(True_Moments, AllSimOccupMoments,ABCcutoff){
  
  #Calculate distances
  Dist<-AllMomDistance(True_Moments,AllSimOccupMoments)
  #Collect indices that have distances greater than cutoff
  rejected_Parms_ind <- which(Dist> ABCcutoff)
  
  #Return simulated moments that are rejected
  AllSimOccupMoments[rejected_Parms_ind,]
}


#===============================================================================
#Function: CovMomEstFun
# - Gives a point estimate and a credible interval of the coefficient of variation
#   in mosquito feeding probabilities by the occupancy spectrum's moments generated within the 
#   function
# Input: The true B and theta of the observed occupancy spectrum, Number of mosquitoes,
#       Values drawn from priors of B and theta (matrix) and a predetermined cutoff 
#       for ABC
# Output: Vector with Point estimate and endpoints of the CI, acceptance rate

CovMomEstFun<- function(B,theta,M,Prior,ABCcutoff){
  NoOfTrials <- nrow(Prior)
  
  #Generates mosquito feeding probability vector
  mu = 1/B
  Shape1 <-theta*mu
  Shape2 <- theta- Shape1
  prob_true <-rbeta(B,Shape1,Shape2)
  
  #Generates Occupancy spectrum from probability vector and calculates its moments
  S_true <-sHetfun2(M,prob_true)
  true_Mom <-getMoments(S_true)
  
  #Generates Occupancy spectrum moments from prior
  AllOccup <- SimulateOccupancyMom(M,Prior,HomogAssump = FALSE)
  
  #Calculates distance from true moments to simulated moments
  Dist <- AllMomDistance(true_Mom,AllOccup)
  
  # Returns posterior after theta's and B's are rejected for distance > ABCcutoff
  Posterior <- GetPosterior(Prior,Dist,ABCcutoff)
  
  #Avoids dimension issue
  if (is.vector(Posterior)){
    BPosterior <- Posterior[1]
    ThetaPosterior <- Posterior[2]
  }
  if (is.matrix(Posterior)){
    BPosterior <- Posterior[,1]
    ThetaPosterior <- Posterior[,2]
  }
  
  #Calculate CoV from joint posterior of B and theta
  CovPosterior <- CoVariationBeta(BPosterior,ThetaPosterior)
  
  
  ##if all entries of the posterior are NA, then Acceptance is 0
  #else, it is the length of the Posterior vector
  if (sum(is.na(CovPosterior)) == length(CovPosterior)){
    Accepted = 0
  }else{
    Accepted = length(CovPosterior)/nrow(Prior)*100
  }
  
  
  #Outputs result: Point estimate and endpoints of interval, and acceptance rate
  Cov_estimate = mean(CovPosterior)
  Cov_CI = quantile(CovPosterior,c(0.025,0.975),na.rm = TRUE)
  
  c(CoVestimate = Cov_estimate
    ,Cov_CI =Cov_CI
    ,Acceptance = Accepted)
  
  
  
}


