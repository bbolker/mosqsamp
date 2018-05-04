#ABC for Theta and B joint distribution- Grid Plot for estimating Theta


source("hetmosqfuns.R")
source("ABCfun.R")
library(emdbook)
library(coda)
library(ggplot2)
library(gridExtra)



#Prior
BThetaPrior3 <- GenPriorBTheta(1000,Bsize = 6, Bp =0.2,TShape = 40,TScale = 0.5)
hist(BThetaPrior3[,"theta"])

#######Trial: Keeping Theta relatively constant, with a thin prior, aim of estimating B



uniqCombo = length(seq(15,60, by =10))

B_parms = rep(seq(15,60, by =10),6)
M_fix = c(rep(25,uniqCombo),rep(50,uniqCombo),rep(100,uniqCombo),rep(200,uniqCombo)
          ,rep(400,uniqCombo),rep(600,uniqCombo))

theta_true = 10


##Vector of cutoffs- these can be modified, the cutoffs used should be approximately this size

some_cutoffs2 <-vector(mode = "numeric",length(B_parms))
for(i in 1:length(B_parms)){
  some_cutoffs2[i]<-genABCcutoff(B_parms[i],M_fix[i],level = 0.1,theta_true = theta_true)
}


PosteriorInfoBHet <- matrix(0,length(B_parms),9,
                            dimnames = list(Trial = 1:length(B_parms)
                                            ,Stat =c("b_trial","theta_trial"
                                                    ,"M","Mean_Estimate"
                                                    ,"Bias","MSE", "CI_width"
                                                    ,"Coverage","Mean_Acceptance")))

nsim =100

PosteriorInfoBHet[,"b_trial"] <-B_parms
PosteriorInfoBHet[,"M"] <- M_fix
PosteriorInfoBHet[,"theta_trial"] <- theta_true

Ttest <- Sys.time()

for(i in 1:length(B_parms)){
  nsim=100
  theta_trial <- PosteriorInfoBHet[i,"theta_trial"]; B_trial <-PosteriorInfoBHet[i,"b_trial"]
  M <- PosteriorInfoBHet[i,"M"]
  
  PosteriorIntervals <-replicate(nsim,thetaBEstFun(B_trial,theta_trial,M,BThetaPrior3,some_cutoffs2[i]))
  
  #remove NAs
  if(sum(is.na(PosteriorIntervals[3,]))>0){
    allrejects <- sum(is.na(PosteriorIntervals[3,]))
    PosteriorIntervals <- PosteriorIntervals[,-which(is.na(PosteriorIntervals[1,]))]
    nsim = nsim -allrejects
  }
  
  #Saving results to matrix
  PosteriorInfoBHet[i,"Mean_Estimate"] <- mean(PosteriorIntervals[1,])
  PosteriorInfoBHet[i,"Bias"] <- (mean(PosteriorIntervals[1,])-B_trial)/B_trial
  PosteriorInfoBHet[i,"MSE"] <- sum((PosteriorIntervals[1,]-B_trial)^2)/nsim/B_trial^2
  PosteriorInfoBHet[i,"CI_width"] = mean(PosteriorIntervals[3,] - PosteriorIntervals[2,])
  PosteriorInfoBHet[i,"Coverage"] <- sum(PosteriorIntervals[2,] <= B_trial 
                                         & B_trial <= PosteriorIntervals[3,])/nsim
  PosteriorInfoBHet[i,"Mean_Acceptance"] <-mean(PosteriorIntervals[7,])
  
  
}

Sys.time()-Ttest

PosteriorInfoBHet.data <- as.data.frame(PosteriorInfoBHet)



Cov <- ggplot(data= PosteriorInfoBHet.data)+
  geom_line(aes(M,Coverage,group=factor(b_trial),colour=factor(b_trial)))+
  geom_hline(aes(yintercept = 0.95),linetype = "dashed")+
  theme(legend.position="none")+
  ylim(0,1)+
  xlab("Mosquitoes")
CI <-ggplot(data= PosteriorInfoBHet.data)+
  geom_line(aes(M,CI_width,group=factor(b_trial),colour=factor(b_trial)))+
  labs(colour = "B")+theme(legend.position="none")+
  ylab("CI width")+
  xlab("Mosquitoes")
Bias <-ggplot(data= PosteriorInfoBHet.data)+
  geom_line(aes(M,Bias, group=factor(b_trial),colour=factor(b_trial)))+
  theme(legend.position="none")+
  geom_hline(aes(yintercept = 0))+theme(legend.position="none")+
  ylab("Relative Bias")+
  xlab("Mosquitoes")
MeanEst <- ggplot(data = PosteriorInfoBHet.data)+
  geom_line(aes(M,Mean_Estimate,group=factor(b_trial),colour=factor(b_trial)))+
  labs(colour = "B")+theme(legend.position="none")+
  ylab("Mean Estimate")+
  xlab("Mosquitoes")
MSE <- ggplot(data= PosteriorInfoBHet.data)+
  geom_line(aes(M,MSE,group=factor(b_trial),colour=factor(b_trial)))+
  geom_hline(aes(yintercept = 0))+
  ylab("Relative MSE")+
  labs(colour="B")+theme(legend.position="none")+
  xlab("Mosquitoes")


MSE2<- ggplot(data= PosteriorInfoBHet.data)+
  geom_line(aes(M,MSE,group=factor(b_trial),colour=factor(b_trial)))+
  labs(colour="B")+theme(legend.position="left")

hetlegend <-g_legend(MSE2)


grid.arrange(Bias,MSE,Cov,CI,MeanEst,hetlegend,ncol=2)
