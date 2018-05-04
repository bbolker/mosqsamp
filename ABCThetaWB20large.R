#ABC for Theta and B joint distribution- Grid Plot for estimating Theta


source("hetmosqfuns.R")
source("ABCfun.R")
library(emdbook)
library(coda)
library(ggplot2)
library(gridExtra)

#Prior
BThetaPrior2 <- GenPriorBTheta(1000,Bsize = 20,Bp=0.95,TShape = 5,TScale = 2)




#######Trial: Keeping B relatively constant, with a thin prior, aim of estimating theta

#Setting values theta and M to repeat with
uniqueTheta =c(2,5,10,15,20)
uniqCombo =length(uniqueTheta)
theta_true = rep(uniqueTheta,6)
#Fixed value
M_fix = c(rep(25,uniqCombo),rep(50,uniqCombo),rep(100,uniqCombo),rep(200,uniqCombo)
          ,rep(400,uniqCombo),rep(600,uniqCombo))
B_true = rep(22,length(theta_true))

#STEP 1: Set Threshold

##Vector of cutoffs- these can be modified, the cutoffs used should be approximately this size
some_cutoffs2 <-vector(mode = "numeric",length(theta_true))
for(i in 1:length(theta_true)){
  some_cutoffs2[i]<-genABCcutoff(B_true[i],M_fix[i],level = 0.1,theta_true = theta_true[i])
}



#Set up Matrix
PosteriorInfoTheta <- matrix(0,length(theta_true),9,
                             dimnames = list(Trial = 1:length(theta_true), 
                                             Stat =c("b_trial","theta_trial"
                                                     ,"M","Mean_Estimate_theta"
                                                     ,"Bias","MSE", "CI_width"
                                                     ,"Coverage","Mean_Acceptance")))

nsim =100

PosteriorInfoTheta[,"b_trial"] <-B_true
PosteriorInfoTheta[,"M"] <- M_fix
PosteriorInfoTheta[,"theta_trial"] <- theta_true

Ttest <- Sys.time()


#Iterate through theta, B and M values
for(i in 1:length(theta_true)){
  nsim =100
  theta_trial <- PosteriorInfoTheta[i,"theta_trial"]; B_trial <-PosteriorInfoTheta[i,"b_trial"]
  M <- PosteriorInfoTheta[i,"M"]
  
  #Generates B and theta point estimates and CIs
  PosteriorIntervals <-replicate(nsim,thetaBEstFun(B_trial,theta_trial,M,BThetaPrior2,some_cutoffs2[i]))
  
  #remove NAs
  if(sum(is.na(PosteriorIntervals[3,]))>0){
    allrejects <- sum(is.na(PosteriorIntervals[3,]))
    PosteriorIntervals <- PosteriorIntervals[,-which(is.na(PosteriorIntervals[1,]))]
    nsim = nsim -allrejects
  }
  
  #Collect statistics- mean estimates, Bias, MSE, CI width, Coverage and also here acceptance
  PosteriorInfoTheta[i,"Mean_Estimate_theta"] <- mean(PosteriorIntervals[4,])
  PosteriorInfoTheta[i,"Bias"] <- (mean(PosteriorIntervals[4,])-theta_trial)/theta_trial
  PosteriorInfoTheta[i,"MSE"] <- sum((PosteriorIntervals[4,]-theta_trial)^2)/nsim/theta_trial^2
  PosteriorInfoTheta[i,"CI_width"] = mean(PosteriorIntervals[6,] - PosteriorIntervals[5,])
  PosteriorInfoTheta[i,"Coverage"] <- sum(PosteriorIntervals[5,] <= theta_trial & theta_trial <= PosteriorIntervals[6,])/nsim
  PosteriorInfoTheta[i,"Mean_Acceptance"] <- mean(PosteriorIntervals[7,])
  
}

#Benchmark
Sys.time()-Ttest

#Convert to data frame
PosteriorInfoTheta.data <- as.data.frame(PosteriorInfoTheta)

##The Graphs

#<<ThetaEstimation,echo=FALSE,warning=FALSE>>=
ThetaCov <-ggplot(data= PosteriorInfoTheta.data)+
  geom_line(aes(M,Coverage,group = factor(theta_trial),colour = factor(theta_trial)))+
  geom_hline(aes(yintercept = 0.95),linetype = "dashed")+
  ylim(0,1)+
  theme(legend.position="none")+
  xlab("Mosquitoes")

ThetaCI <-ggplot(data= PosteriorInfoTheta.data)+
  geom_line(aes(M,CI_width,group = factor(theta_trial),colour = factor(theta_trial)))+
  ylab("CI width")+
  labs(colour = "theta")+theme(legend.position="none")+
  xlab("Mosquitoes")

ThetaBias <- ggplot(data= PosteriorInfoTheta.data)+
  geom_line(aes(M,Bias,group = factor(theta_trial),colour = factor(theta_trial)))+
  geom_hline(aes(yintercept = 0))+
  ylab("Relative Bias")+
  theme(legend.position="none")+
  xlab("Mosquitoes")

ThetaMSE <- ggplot(data= PosteriorInfoTheta.data)+
  geom_line(aes(M,MSE,group = factor(theta_trial),colour = factor(theta_trial)))+
  geom_hline(aes(yintercept = 0))+
  ylab("Relative MSE")+
  labs(colour = "theta")+theme(legend.position="none")+
  xlab("Mosquitoes")

ThetaMeanEst <- ggplot(data = PosteriorInfoTheta.data)+
  geom_line(aes(M,Mean_Estimate_theta,group = factor(theta_trial),colour = factor(theta_trial)))+
  ylab("Mean Estimate")+
  labs(colour = "theta")+theme(legend.position="none")+
  xlab("Mosquitoes")




#get legend grob
ThetaMSE2 <- ggplot(data= PosteriorInfoTheta.data)+
  geom_line(aes(M,MSE,group = factor(theta_trial),colour = factor(theta_trial)))+
  geom_hline(aes(yintercept = 0))+
  ylab("Relative MSE")+
  labs(colour = "theta")+theme(legend.position="top")+
  xlab("Mosquitoes")

thetalegend <- g_legend(ThetaMSE2)



grid.arrange(ThetaBias,ThetaMSE,ThetaCov,ThetaCI,ThetaMeanEst,thetalegend,ncol=2)

