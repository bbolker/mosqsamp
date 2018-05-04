#Posterior Coefficient of Variation - Heterogeneity Definition v2
#Outputs Graph of Posterior from ABC after with one combination of B=25 and theta=5

source("hetmosqfuns.R")

source("ABCfun.R")
library(emdbook)
library(coda)
library(ggplot2)


#true values
B_true = 25
theta_true = 5
#Fixed value
M_fix = 500

mu = 1/B_true
Shape1 = theta_true / B_true
Shape2 = theta_true - Shape1

#Generate beta probabilities
prob_true <- rbeta(B_true,Shape1,Shape2)

#Get true Occupancy spectrum
S_true <- sHetfun2(M_fix,prob_true)

#Prior and data frame
BThetaPrior2 <- GenPriorBTheta(1000,Bsize = 20,Bp=0.9,TShape = 5,TScale = 2)
BThetaPrior2.data = as.data.frame(BThetaPrior2)

#Generate ABC curoffs
cutoff12<-genABCcutoff(B_true,M_fix,level=0.1,theta_true = 5)

#Simulate occupancy
AllOccup <- SimulateOccupancy(M_fix,BThetaPrior2,HomogAssump = FALSE)

#Compute Distances
Dist <- AllDistance(S_true,AllOccup)
#Get Posterior
Posterior2 <- GetPosterior(BThetaPrior2,Dist,cutoff12)

#HPDregion
HPDregionplot(mcmc(BThetaPrior2),col="black",main ="")
HPDregionplot(mcmc(Posterior2),h=30, col = "red",add=TRUE)
points(B_true,theta_true, pch=20,col="blue")
points(mean(BThetaPrior2[,1]),mean(BThetaPrior2[,2]),pch=23,col="black")
points(mean(Posterior2[,1]),mean(Posterior2[,2]),pch=17,col="red")


## Coefficient of Variation in a beta distribution

CoVariationBeta <- function(b,theta){
  ((b-1)/(theta + 1))^0.5
}


# ~CoVariation of the beta distributed capture probabilities, the priors, posteriors, and the true
CoVariationPrior <- CoVariationBeta(BThetaPrior2[,"B"],BThetaPrior2[,"theta"])
CoVariationPosterior <- CoVariationBeta(Posterior2[,"B"],Posterior2[,"theta"])
CoVariation_true <-CoVariationBeta(B_true,theta_true)


HetPrior.data <-as.data.frame(CoVariationPrior)
HetPosterior.data <-as.data.frame(CoVariationPosterior)

#Acceptance rate of simulations in ABC
Acceptance = length(HetPosterior.data$CoVariationPosterior)/length(HetPrior.data$CoVariationPrior)*100

hetTitle= paste("Density estimation of COV using Occupancy, Acceptance = ",Acceptance, "%")

#GoodPlot for results

ggplot()+
  geom_density(data=HetPrior.data,aes(CoVariationPrior,fill="Prior"),alpha=0.7)+
  geom_density(data=HetPosterior.data,aes(CoVariationPosterior,fill="Posterior"),alpha=0.7)+
  geom_vline(xintercept = CoVariation_true,alpha=0.7)+
  xlab("Heteorogeneity")+
  labs(title = hetTitle)







