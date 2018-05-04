##ABC case: Assuming Homogeneous probability

source("ABCfun.R")
source("hetmosqfuns.R")
library(ggplot2)

#When probability is homogeneous and assumed to be:
# -	The only parameter being estimated is B
# - Consider different sizes of B and M




#Generate Posterior
# Given values drawn from the prior, and the true bird population, true theta, a cutoff and some M # of mosquitoes
# outputs the posterior of B
EstHomogPosterior <- function(B,theta=0,M,Prior,ABCcutoff){
  
  if(theta == 0){
    #True occupancy simulated from equal probability
    S_true <- sHetfun2(M,rep(1,B))
  }
  if(theta > 0){
    #True occupancy simulated from Beta distributed probability
    mu = 1/B
    Shape1 = mu*theta
    Shape2 = theta - Shape1
    prob_true = rbeta(B,Shape1,Shape2)
    S_true <- sHetfun2(M,prob_true)
  }
  
  # Occupancy spectrum is simulated given the number of birds drawn from prior
  All <- SimulateOccupancy(M,Prior,HomogAssump = TRUE)
  
  #get distances from all occupancy spectrums to S_true
  Dist <-AllDistance(S_true,All)
  
  #Rejects all values > epsilon
  posterior <- GetPosterior(Prior,Dist,ABCcutoff)
  
  #output
  posterior
}


B_true = 30

#===============================================================================
# Homogeneous Assumption is true

# 1. Draw B_i from prior of B
BPrior <-GenPriorBnBinom(1000,6,0.2)

# Cutoff value is predetermined using 10th percentile of distances
cutoff10<-genABCcutoff(mean(BPrior),500,level = 0.1)

# Returns Posterior
BPosterior <- EstHomogPosterior(B_true,theta=0,M=500,BPrior,cutoff10)

# Convert to data frame for ggplot
BPrior.data <- as.data.frame(BPrior)
BPosterior.data <- as.data.frame(BPosterior)

#===============================================================================
#Homogeneous Assumption is false

# Cutoff value is predetermined using 10th percentile of distances with theta
cutoff10theta <- genABCcutoff(mean(BPrior),500,level = 0.1,theta_true = 5)

# Returns Posterior, generates S_true (the true occupancy) when theta =5
BHetPosterior<- EstHomogPosterior(B_true,theta = 5,M = 500,BPrior,cutoff10theta)

# Convert to data frame for ggplot
BHetPosterior.data = as.data.frame(BHetPosterior)



#Plots
Heterogeneity <- ggplot()+
  geom_density(data=BPrior.data,aes(BPrior,fill="Prior"),alpha=0.7)+
  geom_density(data=BHetPosterior.data,aes(BHetPosterior,fill ="Posterior"),alpha=0.7)+
  geom_vline(xintercept = B_true,alpha=0.7)+
  xlab("B")+
  labs(title ="B")+theme(legend.position="none")

Homogeneous <- ggplot()+
  geom_density(data=BPrior.data,aes(BPrior,fill="Prior"),alpha=0.7)+
  geom_density(data=BPosterior.data,aes(BPosterior,fill="Posterior"),alpha=0.7)+
  geom_vline(xintercept = B_true,alpha=0.7)+
  xlab("B")+ 
  theme(legend.position="none")+
  labs(title = "A")



grid.arrange(Homogeneous, Heterogeneity,ncol= 2)