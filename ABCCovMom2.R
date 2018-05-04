#Posterior Coefficient of Variation - Heterogeneity Definition v2 using moments
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
prob_true <- rbeta(B_true,Shape1,Shape2)

S_true <- sHetfun2(M_fix,prob_true)
true_Mom <-getMoments(S_true)


#PRIOR distribution
BThetaPrior2 <- GenPriorBTheta(1000,Bsize = 20,Bp=0.9,TShape = 5,TScale = 2)
#hist(BThetaPrior2[,"theta"])
#hist(BThetaPrior2[,"B"])

Het2Prior <-CoVariationBeta(BThetaPrior2[,"B"],BThetaPrior2[,"theta"])
Het2Prior.data <- as.data.frame(Het2Prior)


#Predetermined cutoff value
cutoff10<-genMomentCutoff(B_true,M_fix,level=0.1,theta_true = mean(BThetaPrior2[,"theta"]))

#Simulate Occupancy spectrum Moments from Prior
AllOccup <- SimulateOccupancyMom(M_fix,BThetaPrior2,HomogAssump = FALSE)
#Compute distances between simulated and the true occupancy moments
Dist <- AllMomDistance(true_Mom,AllOccup)
#Reject all distances greater than the cutoff
MOMPosterior <- GetPosterior(BThetaPrior2,Dist,cutoff10)


#HPDregionPlot
HPDregionplot(mcmc(BThetaPrior2),col="black",main ="HPDRegions of B and Theta")
HPDregionplot(mcmc(MOMPosterior), col = "red",add=TRUE)
points(B_true,theta_true, pch=20)




# CoVariation in prior, posteror and the true value
CoVariationPrior <- CoVariationBeta(BThetaPrior2[,"B"],BThetaPrior2[,"theta"])
CoVariationPosterior <- CoVariationBeta(MOMPosterior[,"B"],MOMPosterior[,"theta"])
CoVariation_true <-CoVariationBeta(B_true,theta_true)

#Convert to dataframe
HetPrior.data <-as.data.frame(CoVariationPrior)
HetPosterior.data <-as.data.frame(CoVariationPosterior)

#How many simulations were accepted? add to title of graph
Acceptance = length(HetPosterior.data$CoVariationPosterior)/length(HetPrior.data$CoVariationPrior)*100
hi =paste("Posterior Density Estimation using Moments, Acceptance: ",Acceptance, "%")

#GoodPlot for results
ggplot()+
  geom_density(data=HetPrior.data,aes(CoVariationPrior,fill="Prior"),alpha=0.7)+
  geom_density(data=HetPosterior.data,aes(CoVariationPosterior,fill="Posterior"),alpha=0.7)+
  geom_vline(xintercept = CoVariation_true,alpha=0.7)+
  xlab("Heteorogeneity")+
  labs(title = hi)




