#Posterior Coefficient of Variation - Heterogeneity Definition v2 with moments- 2x2 +1 grid
# large scale test - shows Bias, MSE,CI width, Coverage, Mean estimate from estimation with several combinations
#                    of B, theta and M

source("hetmosqfuns.R")
source("ABCfun.R")

library(emdbook)
library(coda)
library(ggplot2)
theme_set(theme_bw())
library(gridExtra)

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



#PRIORS
BThetaPrior2 <- GenPriorBTheta(1000,Bsize = 20,Bp=0.95,TShape = 5,TScale = 2)
#hist(BThetaPrior2[,"theta"])
#hist(BThetaPrior2[,"B"])
BThetaPrior2.data = as.data.frame(BThetaPrior2)


#######Trial Keeping B relatively constant, with a thin prior, aim of estimating Coefficient of Variation with moments



#Setting values theta and M to repeat ABC with
uniqueTheta =c(2,5,10,15,20)
uniqCombo =length(uniqueTheta)
theta_true = rep(uniqueTheta,6)
M_fix = c(rep(25,uniqCombo),rep(50,uniqCombo),rep(100,uniqCombo),rep(200,uniqCombo),rep(400,uniqCombo),rep(600,uniqCombo))
#B is kept constant
B_true = rep(20,length(theta_true))



##Vector of cutoffs- these can be modified, the cutoffs used should be approximately this size
#or be to the same order of magnitude
some_cutoffs <-vector(mode = "numeric",length(theta_true))
for(i in 1:length(theta_true)){
  some_cutoffs[i]<-genMomentCutoff(B_true[i],M_fix[i],level = 0.15,theta_true = theta_true[i])
}

#Set up matrix to collect statistics: Bias, MSE, Mean estimate, CI width, Coverage, etc.
PosteriorInfoMomCOV <- matrix(0,length(theta_true),10
                           ,dimnames = list(Trial = 1:length(theta_true),Stat = c("Cov_trial","b_trial","theta_trial","M", "COVMeanEst","Bias","MSE","CI_width","Coverage","Mean_Acceptance")))

nsim = 100

PosteriorInfoMomCOV[,"b_trial"]<-B_true
PosteriorInfoMomCOV[,"M"]<- M_fix
PosteriorInfoMomCOV[,"theta_trial"]<-theta_true
PosteriorInfoMomCOV[,"Cov_trial"]<-CoVariationBeta(B_true,theta_true)




##Vector of cutoffs- these can be modified, the cutoffs used should be approximately this size
some_cutoffs <-vector(mode = "numeric",length(theta_true))
for (i in 1:length(theta_true)){
  some_cutoffs[i]<- genMomentCutoff(B_true[i],M_fix[i],level = 0.1,theta=theta_true[i])
}

#For benchmark
Ttest <- Sys.time()

#Iterate through values of theta and B and M to repeat the MSE with
for( i in 1:length(theta_true)){
  nsim = 100
  theta_trial <-PosteriorInfoMomCOV[i,"theta_trial"]; B_trial <- PosteriorInfoMomCOV[i,"b_trial"]
  M <- PosteriorInfoMomCOV[i,"M"]
  
  COV_trial<-PosteriorInfoMomCOV[i,"Cov_trial"]
  
  #Replicate Estimation function 100 times
  PosteriorIntervals <- replicate(nsim,CovMomEstFun(B_trial,theta_trial,M,BThetaPrior2,some_cutoffs[i]))
  
  #remove NAs
  if(sum(is.na(PosteriorIntervals[3,]))>0){
    allrejects <- sum(is.na(PosteriorIntervals[3,]))
    PosteriorIntervals <- PosteriorIntervals[,-which(is.na(PosteriorIntervals[1,]))]
    nsim = nsim -allrejects
  }
  
  #Collecting statistics: Mean estimate, Bias, MSE, CI width, Coverage and Acceptance
  PosteriorInfoMomCOV[i,"COVMeanEst"]<- mean(PosteriorIntervals[1,])
  PosteriorInfoMomCOV[i,"Bias"] <-(mean(PosteriorIntervals[1,])-COV_trial)/COV_trial
  PosteriorInfoMomCOV[i,"MSE"] <-sum((PosteriorIntervals[1,]-COV_trial)^2)/nsim/COV_trial^2
  PosteriorInfoMomCOV[i,"CI_width"] <- mean(PosteriorIntervals[3,] - PosteriorIntervals[2,])
  PosteriorInfoMomCOV[i,"Coverage"] <- sum(PosteriorIntervals[2,] <= COV_trial & COV_trial <= PosteriorIntervals[3,])/nsim
  PosteriorInfoMomCOV[i,"Mean_Acceptance"] <- mean(PosteriorIntervals[4,])
  
  
}

Sys.time()-Ttest

PosteriorMomCOV.data<-as.data.frame(PosteriorInfoMomCOV)

##The Graphs

BiasMomCov <-ggplot(data = PosteriorMomCOV.data)+
  geom_line(aes(M,Bias,group=factor(round(Cov_trial,3)),colour=factor(round(Cov_trial,3))))+
  geom_hline(aes(yintercept = 0))+
  ylab("Relative Bias")+
  xlab("Mosquitoes")+
  labs(title="A",colour="C. of Variation")+
  theme(legend.position="none")

MSEMomCov <-ggplot(data = PosteriorMomCOV.data)+
  geom_line(aes(M,MSE,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  geom_hline(aes(yintercept=0))+
  ylab("Relative MSE")+
  xlab("Mosquitoes")+
  labs(title="B", colour = "C. of Variation")+
  theme(legend.position="none")

CovMomCov <-ggplot(data = PosteriorMomCOV.data)+
  geom_line(aes(M,Coverage,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  geom_hline(aes(yintercept=0.95),linetype = "dashed")+
  ylim(0,1)+
  labs(title="C", colour = "C. of Variation")+
  xlab("Mosquitoes")+
  theme(legend.position="none")

CIMomCov <-ggplot(data = PosteriorMomCOV.data)+
  geom_line(aes(M,CI_width,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  labs(title="D",colour = "C. of Variation")+
  xlab("Mosquitoes")+
  ylab("CI width")+
  theme(legend.position="none")

MeanEstMomCov<-ggplot(data=PosteriorMomCOV.data)+
  geom_line(aes(M,COVMeanEst,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  labs(title="E",colour = "C. of Variation")+
  xlab("Mosquitoes")+
  ylab("Mean Estimate")+
  theme(legend.position="none")


#extract legend grob from this graph
MSE2 <-ggplot(data=PosteriorMomCOV.data)+
  geom_line(aes(M,COVMeanEst,group=factor(round(Cov_trial,3)),colour =factor(round(Cov_trial,3))))+
  labs(title="E",colour = "Coefficient of Variation")+
  theme(legend.position = "right")

hetlegend<-g_legend(MSE2)

Hetgrid <-grid.arrange(BiasMomCov,MSEMomCov,CovMomCov,CIMomCov,MeanEstMomCov,hetlegend,ncol = 2)

ggsave("MOMCOVgrid.png",Hetgrid,dpi= 600,units ="in",width=5)



