##version January 2/2018
source("mosqfuns2.R")

##==============================================================================

##Table of Contents:
#1. define heterogeneous probability vectors using geometric form where:
#form: p^0 p p^2 p^3 p^4 p^5 ... p^(B-1)
##values of p: 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0

#2. Key Estimator functions
#       - Chao Estimator - Homogeneous, Heterogeneous
#       - Doublet Estimator - Homogeneous, Heterogeneous

#3. Simulate a test of the chao estimator - generate array and data frame ChaoTest

#will produce content for plotting: hetmosqbatch.WKPriorsMatrices.Rdata
#       -

#4. Test the 4 proposed estimators
#       a. chao  - Doublet estimator
#       b. (chao  - Doublet estimator)/chao1
#       c. (chao  - Doublet estimator)/doublet estimator
#       d. |chao  - Doublet estimator|
#===============================================================================
#Changes to original:


#===============================================================================
#               Preliminaries
#===============================================================================
library(fossil)
library(reshape2)
library(ggplot2)


#===============================================================================
#               1. Heterogeneous Probability vectors
#===============================================================================

## fill a matrix with the probability vectors in each column
## need to set B beforehand
B<- 40; ps <- seq(1,.1,by=-0.1) 
##empty matrix with B rows and length p columns
geoProbMat <- matrix(nrow = B,ncol=length(ps))

##create cooresponding metric of probability heterogeneity
probGeoHeterog <- rep(0,length(ps))

##alternatively can use p
probGeoHeterog2 <- ps

##populate each of the p probability columns
for (j in 1:length(ps)){
  ##column j will be p^0 p p^2 p^3 p^4 p^5 ... p^(B-1)
  geoProbMat[,j] = ps[j]^(0:(B-1))
  ##column j is then divided by its sum to normalize it
  geoProbMat[,j] = geoProbMat[,j]/sum(geoProbMat[,j])
  
  ##fill the corresponding heterogeneity vector with sd
  probGeoHeterog[j] = sd(geoProbMat[,j])
  
}
#note: the last column is all 1's - the homogeneous probability case

#plot(probGeoHeterog)

##### 0.14261134 0.12826921 0.11474180 0.10172505 0.08891559 0.07595545
#[7] 0.06232068 0.04699718 0.02736777 0.00000000
#decreasing and linear


#plot(probGeoHeterog2)

# 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
#increasing and linear


#===============================================================================
#               2. Key Estimator functions
#===============================================================================


#===============================================================================
#VHetfunction for drawing a multinomial sample with heterogeneous
#probabilities

vHetfun <- function(M,prob) { 
  c(rmultinom(1,size=M,prob=prob))
}

## One realization of Chao using Homogenous probability
#===============================================================================


Chaofun <- function (B,M){
  v <- vfun(B,M)
  chao1(v)
}

## One realization of Chao using custom probability
#===============================================================================

ChaoHetFun <- function(M,prob){
  v <-vHetfun(M,prob)
  chao1(v)
}

#===============================================================================
# Doublet estimator with W and M inputs, returns Bhat

Bhat_doublet_fun <- function (W,M){
  M*(M-1)/W
}

## Testing Doublet using custom Probability
#===============================================================================

WEstfun<-function(B,M){
  v <- vfun(B,M)
  W <-sum(v*(v-1))
  Bhat_doublet_fun(W,M)
}

## Testing Doublet using custom Probability
#===============================================================================

WHetFun <- function (M, prob) {
  v <- vHetfun(M,prob)
  W <-sum(v*(v-1))
  Bhat_doublet_fun(W,M)
}


#===============================================================================
#               3. Test Chao
#===============================================================================

###
#definitions

#relative bias: (Estimate - Btrue)/Btrue

# MSE of estimate over n realizations: sum ((estimate - Btrue) ^2 )/ n

# relative MSE : MSE of estimate / B^2

#jan 3/18
#Jan 9/18

##Define all of the parameters
nsim<-1000;
ps <- seq(1,.1,by=-0.1) 
Bs <- c(32,56,100,178,316,562,1000)

Ms = seq(40,1000,by=50)

#ps<- c(0.5,1)
#Bs <- c(20,40)
#Ms <- c(10,30)


#produce empty array with B x M x p x Estimator x Stat
Chaotest <- array(0, c(length(Bs),length(Ms),length(ps),2,2),dimnames = list(B=Bs,M=Ms,p=ps,Estimator=c("Chao", "Doublet"),Stat=c("Bias","MSE")))

#Chaotest <- array(0, c(length(Bs),length(Ms),length(ps),2),dimnames = list(B=Bs,M=Ms,p=ps,Stat=c("Bias","MSE")))




for (b in 1:length(Bs)){
  #choose value of B
  B = Bs[b]
  #define the geo probability matrix according to said value
  geoProbMat <- matrix(nrow = B,ncol=length(ps))
  
  for (j in 1:length(ps)){
    #each column takes on p^0 p^1 p^2 .... p^(B-1)
    geoProbMat[,j]= ps[j]^(0:(B-1))
    geoProbMat[,j] = geoProbMat[,j]/sum(geoProbMat[,j])
    
    #associated heterogeneity vector
    probGeoHeterog[j] = sd(geoProbMat[,j])
    ##add another dimension for heterogeneity?
  }
  
  for (p in 1:length(ps)) {
    #iterate through probability vectors
    for (m in 1:length(Ms)){
      #iterate through M values
      
      ##choose M value
      M = Ms[m]
      
      #Simulate Chao estimate
      EstChao <- replicate(nsim,ChaoHetFun(M,geoProbMat[,p]))
      
      #Calculte Bias and MSE for Chao Estimates
      EstChaoBias <- mean((EstChao -B)/B)
      EstChaoRelMSE <- sum((EstChao - B)^2)/nsim
      
      #Assignment
      Chaotest[b,m,p,"Chao","Bias"]<- EstChaoBias
      Chaotest[b,m,p,"Chao","MSE"]<-EstChaoRelMSE
      
      #Chaotest[b,m,p,"Bias"]<- EstChaoBias
      #Chaotest[b,m,p,"MSE"]<-EstChaoRelMSE
      
      #Simulate Doublet estimates
      EstDoublet <- replicate(nsim,WHetFun(M,geoProbMat[,p]))
      
      #Calculate Bias and MSE for Chao Estimates
      EstDoubletBias <- mean((EstDoublet-B)/B)
      EstDoubletRelMSE <- sum((EstDoublet - B)^2)/nsim
      
      Chaotest[b,m,p,"Doublet","Bias"] <- EstDoubletBias
      Chaotest[b,m,p,"Doublet", "MSE"] <- EstDoubletRelMSE
      
    }
  }
  }


ChaoTest.data <- na.omit(melt(Chaotest))  

##Calculating the corresponding Heterogeneity - funtion of B and p  
geoHetMat <- matrix(nrow = length(Bs),ncol=length(ps),dimname=list(B=Bs,p=ps))

for (i in 1:length(Bs)){
  B = Bs[i]
  geoProbMat <- matrix(nrow = B,ncol=length(ps))
  for (j in 1:length(ps)){
    #each column takes on p^0 p^1 p^2 .... p^(B-1)
    geoProbMat[,j]= ps[j]^(0:(B-1))
    geoProbMat[,j] = geoProbMat[,j]/sum(geoProbMat[,j])
    
    #associated heterogeneity vector
    geoHetMat[i,j] = sd(geoProbMat[,j])
    ##add another dimension for heterogeneity?
  }
}

B = ChaoTest.data$B
p = ChaoTest.data$p

char_B <- as.character(B)
char_p <- as.character(p)
het <- geoHetMat[char_B[1],char_p[1]]
for(i in 2:length(B)){
  het = c(het, geoHetMat[char_B[i],char_p[i]])
  
}

ChaoTest.data$Heterogeneity <- het


##Visualization of Chao Test



#plotting all together
ggplot(data= ChaoTest.data) + 
  geom_point(mapping = aes(x=Heterogeneity,y=value,group = M, color=M),subset(ChaoTest.data, Estimator == "Chao") ) +
  geom_line(mapping= aes(x=Heterogeneity,y=value,group = M, color = M),subset(ChaoTest.data,Estimator == "Doublet")) +
  facet_grid(Stat ~ B)

#Only a couple values of p, but comparing Chao and Doublet side by side
BiasChao.data <- subset(ChaoTest.data, Stat != "MSE")
ChaoVrsDoublet.data <- subset(BiasChao.data, p == c(1.0,0.9))



#Create Data column that only says Heterogeneous vrs Homogeneous
temp<-ChaoVrsDoublet.data$p
compare <- temp == 1.0
Heterogeneity <- rep("Heterogeneous",length(temp))
for (i in 1:length(Heterogeneity)){
  if (compare[i]){
    Heterogeneity[i]="Homogeneous"
  }
}
ChaoVrsDoublet.data$HeterogeneityType<-Heterogeneity


ggplot(data=ChaoVrsDoublet.data)+
  geom_line(mapping=aes(x=M,y=value,color=Estimator))+
  scale_x_log10()+
  labs(y="")+
  facet_grid(HeterogeneityText~B)



##Only the bias
BiasChao.data <- subset(ChaoTest.data, Stat != "MSE")




ggplot(data= BiasChao.data) + 
  geom_line(linetype = "dotdash",mapping = aes(x=M,y=value,group = Heterogeneity, color=Heterogeneity),subset(BiasChao.data, Estimator == "Chao") ) +
  geom_line(mapping= aes(x=M,y=value,group = Heterogeneity, color = Heterogeneity),subset(BiasChao.data,Estimator == "Doublet")) +
  facet_grid(Stat ~ B)

#Separated by alternative metric
ggplot(data= BiasChao.data) + 
  geom_point(mapping = aes(x=M,y=value,group = p, color=p),subset(BiasChao.data, Estimator == "Chao") ) +
  geom_line(mapping= aes(x=M,y=value,group = p, color = p),subset(BiasChao.data,Estimator == "Doublet")) +
  facet_grid(Stat ~ B)

#Separated by alternative metric  ##Only Chao
ggplot(data= BiasChao.data) + 
  geom_line(mapping = aes(x=M,y=value,group = p, color=p),subset(BiasChao.data, Estimator == "Chao") ) +
  facet_grid(Stat ~ B)

#Separated by alternative metric
ggplot(data= BiasChao.data) + 
  geom_point(mapping = aes(x=B,y=value,group = p, color=p),subset(BiasChao.data, Estimator == "Chao") ) +
  geom_line(mapping= aes(x=B,y=value,group = p, color = p),subset(BiasChao.data,Estimator == "Doublet")) +
  facet_grid(Stat ~ M)

##Only the MSE
MSEChao.data <- subset(ChaoTest.data, Stat == "MSE")
ggplot(data= MSEChao.data) + 
  geom_point(mapping = aes(x=M,y=value,group = Heterogeneity, color=Heterogeneity),subset(MSEChao.data, Estimator == "Chao") ) +
  geom_line(mapping= aes(x=M,y=value,group = Heterogeneity, color = Heterogeneity),subset(MSEChao.data,Estimator == "Doublet")) +
  facet_grid(Stat ~ B)


#===============================================================================
#               4. Testing the Proposed Estimators
#===============================================================================


#               4A. Preliminary Functions
#===============================================================================


SimProEstimate1 <- function (M,prob){
  v <-vHetfun(M,prob)
  W <- sum(v*(v-1))
  chao1(v) - Bhat_doublet_fun(W,M)
}

SimProEstimate2 <- function (M,prob){
  v <-vHetfun(M,prob)
  W <- sum(v*(v-1))
  (chao1(v) - Bhat_doublet_fun(W,M))/chao1(v)
}

SimProEstimate3 <- function (M,prob){
  v <-vHetfun(M,prob)
  W <- sum(v*(v-1))
  (chao1(v) - Bhat_doublet_fun(W,M))/Bhat_doublet_fun(W,M)
}

SimProEstimate4 <- function (M,prob){
  v <-vHetfun(M,prob)
  W <- sum(v*(v-1))
  abs(chao1(v)-Bhat_doublet_fun(W,M))
}



#               4B. Creating Test Dataset
#===============================================================================

##Define all of the parameters
nsim<-1000;
ps <- seq(1,.1,by=-0.1) 
Bs <- c(32,56,100,178,316,562,1000)

Ms = seq(40,1000,by=50)

#ps<- c(0.5,1)
#Bs <- c(20,40)
#Ms <- c(10,30)


#produce empty array with B x M x p x Proposed1-4 x Stat
ProposedTest <- array(0, c(length(Bs),length(Ms),length(ps),4,2),dimnames = list(B=Bs,M=Ms,p=ps,Estimator=c("No1", "No2", "No3", "No4"),Stat=c("Mean Estimate","Variation in Estimate")))


for (b in 1:length(Bs)){
  #choose value of B
  B = Bs[b]
  #define the geo probability matrix according to said value
  geoProbMat <- matrix(nrow = B,ncol=length(ps))
  
  for (j in 1:length(ps)){
    #each column takes on p^0 p^1 p^2 .... p^(B-1)
    geoProbMat[,j]= ps[j]^(0:(B-1))
    geoProbMat[,j] = geoProbMat[,j]/sum(geoProbMat[,j])
    
    #associated heterogeneity vector
    probGeoHeterog[j] = sd(geoProbMat[,j])
    ##add another dimension for heterogeneity?
  }
  
  for (p in 1:length(ps)) {
    #iterate through different probability vectors
    for (m in 1:length(Ms)){
      #iterate through M values
      
      ##choose M value
      M = Ms[m]
      
      #Simulate Proposed Estimates
      EstProposed1 <- replicate(nsim,SimProEstimate1(M,geoProbMat[,p]))
      EstProposed2 <- replicate(nsim,SimProEstimate2(M,geoProbMat[,p]))
      EstProposed3 <- replicate(nsim,SimProEstimate3(M,geoProbMat[,p]))
      EstProposed4 <- replicate(nsim,SimProEstimate4(M,geoProbMat[,p]))
      
      #Take Mean
      ProposedTest[b,m,p,"No1","Mean Estimate"] <- mean(EstProposed1)
      ProposedTest[b,m,p,"No2","Mean Estimate"] <- mean(EstProposed2)
      ProposedTest[b,m,p,"No3","Mean Estimate"] <- mean(EstProposed3)
      ProposedTest[b,m,p,"No4","Mean Estimate"] <- mean(EstProposed4)
      
      #Take Variance
      ProposedTest[b,m,p,"No1","Variation in Estimate"] <- var(EstProposed1)
      ProposedTest[b,m,p,"No2","Variation in Estimate"] <- var(EstProposed2)
      ProposedTest[b,m,p,"No3","Variation in Estimate"] <- var(EstProposed3)
      ProposedTest[b,m,p,"No4","Variation in Estimate"] <- var(EstProposed4)
      
      

    }
  }
}

#Melt into Dataset
EstimateTest.data <- na.omit(melt(ProposedTest)) 



##Calculating the corresponding Heterogeneity - only depends on B and p N
geoHetMat <- matrix(nrow = length(Bs),ncol=length(ps),dimname=list(B=Bs,p=ps))


for (i in 1:length(Bs)){
  B = Bs[i]
  geoProbMat <- matrix(nrow = B,ncol=length(ps))
  for (j in 1:length(ps)){
    #each column takes on p^0 p^1 p^2 .... p^(B-1)
    geoProbMat[,j]= ps[j]^(0:(B-1))
    geoProbMat[,j] = geoProbMat[,j]/sum(geoProbMat[,j])

    
    #associated heterogeneity vector
    geoHetMat[i,j] = sd(geoProbMat[,j])
    ##add another dimension for heterogeneity?
  }
}
geoProbMat

B = EstimateTest.data$B
p = EstimateTest.data$p

char_B <- as.character(B)
char_p <- as.character(p)
het <- geoHetMat[char_B[1],char_p[1]]
for(i in 2:length(B)){
  het = c(het, geoHetMat[char_B[i],char_p[i]])
  
}

EstimateTest.data$Heterogeneity <- het

Estimate1Test.data <- subset(EstimateTest.data,Stat == "Mean Estimate")
EstimateAllTest.data <- subset(Estimate1Test.data, B %in% c(56,100))
EstimateAllTestB.data <- subset(Estimate1Test.data, M %in% c(90,340,640,890))

