#Posterior of Coefficient of Variation - Heterogeneity Definition v2
# large scale test - shows Bias, MSE,CI width, Coverage, Mean estimate from estimation with several combinations
#                    of B, theta and M

source("hetmosqfuns.R")
source("ABCfun.R")
library(emdbook)
library(coda)
library(ggplot2)

##Preliminary Functions

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


#PRIORS
BThetaPrior2 <- GenPriorBTheta(1000,Bsize = 20,Bp=0.95,TShape = 5,TScale = 2)
#hist(BThetaPrior2[,"theta"])
#hist(BThetaPrior2[,"B"])
BThetaPrior2.data = as.data.frame(BThetaPrior2)



#######Trial Keeping B constant, with a thin prior, aim of estimating Coefficient of Variation


#factorial experiment: Setting values theta and M to repeat ABC with
uniqueTheta =c(2,5,10,15,20)
uniqCombo =length(uniqueTheta)
theta_true = rep(uniqueTheta,6)
M_fix = c(rep(25,uniqCombo),rep(50,uniqCombo),rep(100,uniqCombo),rep(200,uniqCombo),rep(400,uniqCombo),rep(600,uniqCombo))
#True B is constant
B_true = rep(20,length(theta_true))


##Vector of cutoffs- these can be modified, the cutoffs used should be approximately this size
#or be to the same order of magnitude
some_cutoffs <-vector(mode = "numeric",length(theta_true))
for(i in 1:length(theta_true)){
  some_cutoffs[i]<-genABCcutoff(B_true[i],M_fix[i],level = 0.15,theta_true = theta_true[i])
}

#Set up Matrix that holds statistics: Bias, MSE, CI width, Coverage etc.
PosteriorInfoCOV <- matrix(0,length(theta_true),9
                           ,dimnames = list(Trial = 1:length(theta_true),Stat = c("Cov_trial","b_trial","theta_trial","M", "COVMeanEst","Bias","MSE","CI_width","Coverage")))

nsim = 100

PosteriorInfoCOV[,"b_trial"]<-B_true
PosteriorInfoCOV[,"M"]<- M_fix
PosteriorInfoCOV[,"theta_trial"]<-theta_true
PosteriorInfoCOV[,"Cov_trial"]<-CoVariationBeta(B_true,theta_true)

#For benchmarking
Ttest <- Sys.time()


#Iterate through all the possible thetas, Bs and M
for( i in 1:length(theta_true)){
  nsim = 100
  theta_trial <-PosteriorInfoCOV[i,"theta_trial"]; B_trial <- PosteriorInfoCOV[i,"b_trial"]
  M <- PosteriorInfoCOV[i,"M"]
  
  COV_trial<-PosteriorInfoCOV[i,"Cov_trial"]
  
  #Repeat estimation of Cov 100 times
  PosteriorIntervals <- replicate(nsim,CovEstFun(B_trial,theta_trial,M,BThetaPrior2,some_cutoffs[i]))
  
  #remove NAs
  if(sum(is.na(PosteriorIntervals[3,]))>0){
    PosteriorIntervals <- PosteriorIntervals[,-which(is.na(PosteriorIntervals[1,]))]
    nsim = nsim -1
  }
  #Calculate statistic: Mean estimate, bias, MSE, CI width, Coverage
  PosteriorInfoCOV[i,"COVMeanEst"]<- mean(PosteriorIntervals[1,])
  PosteriorInfoCOV[i,"Bias"] <-(mean(PosteriorIntervals[1,])-COV_trial)/COV_trial
  PosteriorInfoCOV[i,"MSE"] <-sum((PosteriorIntervals[1,]-COV_trial)^2)/nsim/COV_trial^2
  PosteriorInfoCOV[i,"CI_width"] <- mean(PosteriorIntervals[3,] - PosteriorIntervals[2,])
  PosteriorInfoCOV[i,"Coverage"] <- sum(PosteriorIntervals[2,] <= COV_trial & COV_trial <= PosteriorIntervals[3,])/nsim
  
  
  
}

#Returns time after last set of computation
Sys.time()-Ttest

PosteriorCOV.data<- as.data.frame(PosteriorInfoCOV)

##The Graphs

BiasCov <-ggplot(data = PosteriorCOV.data)+
  geom_line(aes(M,Bias,group=factor(round(Cov_trial,3)),colour=factor(round(Cov_trial,3))))+
  geom_hline(aes(yintercept = 0))+
  ylab("Relative Bias")+
  xlab("Mosquitoes")+
  labs(title="A",colour="C. of Variation")+
  theme(legend.position="none")

MSECov <-ggplot(data = PosteriorCOV.data)+
  geom_line(aes(M,MSE,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  geom_hline(aes(yintercept=0))+
  ylab("Relative MSE")+
  xlab("Mosquitoes")+
  labs(title="B", colour = "C. of Variation")+
  theme(legend.position="none")

CovCov <-ggplot(data = PosteriorCOV.data)+
  geom_line(aes(M,Coverage,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  geom_hline(aes(yintercept=0.95),linetype = "dashed")+
  ylim(0,1)+
  labs(title="C", colour = "C. of Variation")+
  xlab("Mosquitoes")+
  theme(legend.position="none")

CICov <-ggplot(data = PosteriorCOV.data)+
  geom_line(aes(M,CI_width,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  labs(title="D",colour = "C. of Variation")+
  xlab("Mosquitoes")+
  ylab("CI width")+
  theme(legend.position="none")

MeanEstCov <-ggplot(data=PosteriorCOV.data)+
  geom_line(aes(M,COVMeanEst,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  labs(title="E",colour = "C. of Variation")+
  xlab("Mosquitoes")+
  ylab("Mean Estimate")+
  theme(legend.position="none")


MSE2 <-ggplot(data=PosteriorCOV.data)+
  geom_line(aes(M,COVMeanEst,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  labs(title="E",colour = "Coefficient of Variation")+
  theme(legend.position = "right")

hetlegend<-g_legend(MSE2)

Hetgrid <-grid.arrange(BiasCov,MSECov,CovCov,CICov,MeanEstCov,hetlegend,ncol = 2)

#ggsave("COVgrid.png",Hetgrid,dpi= 600,units ="in")
