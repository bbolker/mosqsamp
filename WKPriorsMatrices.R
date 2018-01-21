##version January 9-18/2018
source("mosqfuns2.R")

##==============================================================================================

##Contents:
#1. B given K Distributions
#2. The Number of Doublets

#will produce content for plotting: hetmosqbatch.WKPriorsMatrices.Rdata

#===============================================================================
#Changes to original:


#===============================================================================
#               Preliminaries
#===============================================================================
library(ggplot2)




#===============================================================================
#               1. The Number of Distinct Birds (K)
#===============================================================================

sfun <- function(B,M,collapse=TRUE) {
  tt <- table(factor(vfun(B,M),levels=0:M))
  if (!collapse) tt else paste(tt,collapse=".")
}
Kfun <- function(B,M) {
  unname(B-sfun(B,M,collapse=FALSE)[1])
}


nsim <- 1000; M <- 20; B <- 40
Kvec <- replicate(nsim,Kfun(B,M))

#Distribution of K for B = 40
Ktab <- table(Kvec)/nsim
#plot(Ktab)


#Creating the P(B|K) matrix where each column represents a different value of K,
# and each row a different B

#Setting parameters: Potential Bird Populations, Fix Number of Mosquitoes, and nsim
BPossibles <- 10:50
M <- 50
nsim <- 1000
KPossibles <- 1:M

PofB<-rep(1/length(BPossibles),length(BPossibles)) 
#assume: P(B) vector is discrete uniform
##Next step: Try another assumption

# Empty P(K|B) matrix
KgivenB <- matrix(0,length(BPossibles),M,dimnames = list(B = BPossibles,K = KPossibles))

# populate P(K|B)
for (b in 1:length(BPossibles)){
  B <- BPossibles[b]
  Kvec <- replicate(nsim,Kfun(B,M))
  Ktab <- table(Kvec)/nsim
  k.char <- names(Ktab)
  k.integer <- as.integer(k.char)
  KgivenB[b,k.integer] = Ktab[k.char]
}

# Calculate P(K) = P(B) X P(K|B) //works out to matrix product
PofK <- PofB %*% KgivenB


#Populate P(B|K)
BgivenK = matrix(0,length(BPossibles),M,dimnames=list(B=BPossibles,K=1:M))
for (j in 1:M) {
  BgivenK[,j]<-PofB*KgivenB[,j]/PofK[,j]
}

K.array <- array(0,c(length(BPossibles),length(KPossibles),1),dimnames = list(B=BPossibles,K=KPossibles,Probability="B | K"))

Bnames<- as.character(BPossibles)
Knames<- as.character(KPossibles)

for (k in Knames){
  K.array[Bnames,k,1]<-BgivenK[,k]
}

#Convert to data frame
K.data <- na.omit(melt(K.array))

##Next Step: work on visualizing


K.data30 <- subset(K.data,K==30)

##Calculation of Mean, Variance for K==30 Graph

#mean
meanK <-sum(K.data30$value*K.data30$B)
meanK
#Variance
varK <-sum(K.data30$value*(K.data30$B)^2) - meanK^2
varK



M = 50
K = 30
## method of moments estimator for K
Bhat_approx <- function(K,M,type=c("frac","exp"),upper=1e7) {
  type <- match.arg(type) ## "frac" by default
  if (K>=M) stop("K must be < M")
  tmpf <- if (type=="exp") {
    function(B) B*(1-exp(-M/B))-K
  } else {
    function(B) B*(1-(1-(1/B))^M)-K
  }
  r <- try(uniroot(tmpf,interval=c(K+0.001,upper))$root)
  if (inherits(r,"try-error")) {
    r <- NA
  }
  return(r)
}

Bhat_mm <- sapply(K,Bhat_approx,M=M)
#Give estimate by method of moments
#Bhat_mm






#===============================================================================
#               2. The Number of Doublets
#===============================================================================

##Simulate one realization of W with Homogeneous probability
Wfun <- function (B,M){
  v <- vfun(B,M)
  sum(v*(v-1)) ###  
}

##Estimate B from W and M
Bhat_doublet_fun <- function (W,M){
  M*(M-1)/W
}




#simsetup for Matrix
nsim <- 1000; M <- 20; B <- 40
Wvec <- replicate(nsim,Wfun(B,M))
Wtab <- table(Wvec)/nsim

Wvec2 <- as.numeric(names(Wtab))


Bhat_doublet <- Bhat_doublet_fun(Wvec2,M)

#Doublet statistic and Bhat
plot(Wvec2,Bhat_doublet,xlab="doublet statistic",ylab="estimated B", main="Doublet Statistic and Bhat for M=20")

Wtab1 <- Wtab
names(Wtab1)=sprintf("%6.3f", Bhat_doublet)
plot(Wtab1, main = "Bhat distribution when actual B = 40")

Wtab2<-Wtab1
names(Wtab2)=sprintf("%6.f", round(Bhat_doublet))
plot(Wtab2, main = "Rounded Bhat distribution when actual B = 40")



#CREATING THE MATRIX OF P(B|W)
#Setting parameters: Potential Bird Populations, Fix Number of Mosquitoes, and nsim
BPossibles <- 10:50
M <- 50
nsim <- 1000
#P(B) vector assumed to be Discrete uniform
PofB<-rep(1/length(BPossibles),length(BPossibles)) 

#W has to be an even number
WPossibles<-seq(from=0, to = M*(M-1), by = 2)

#Create empty matrix
WgivenB = matrix(0,length(BPossibles),length(WPossibles),dimnames=list(B=BPossibles,W=WPossibles))

#fill
for (i in 1:length(BPossibles)){
  B<- BPossibles[i]
  Wvec <- replicate(nsim,Wfun(B,M))
  Wtab <- table(Wvec)/nsim
  w <- names(Wtab)
  WgivenB[i,w]<-Wtab[w]
}

##Remove all the unnecessary columns - these have W value which are impossible or unlikely
indices <-0
for (i in 1:length(WPossibles)){
  
  ##Will be a vector of only TRUE if all members of a column are zero
  allZeros <- WgivenB[,i]==rep(0,length(BPossibles)) 
  #indices are collected of non-zero columns
  if(FALSE %in% allZeros){
    indices <- c(indices,i)
  }
}
#All non-empty indices replace the old matrix
WgivenB<- WgivenB[,indices]
#Remove all W values which didn't occur
WPossibles<- WPossibles[indices]

BgivenW = matrix(0,length(BPossibles),length(WPossibles),dimnames=list(B=BPossibles,W=WPossibles))

PofWtimesB = PofB %*%WgivenB

for (k in 1:length(WPossibles)){
  BgivenW[,k]<-PofB*WgivenB[,k]/PofWtimesB[,k]
}

#BgivenW

Doublet.array <- array(0,c(length(BPossibles),length(WPossibles),1),dimnames = list(B=BPossibles,W=WPossibles,Probability="P(B | W)"))

Bnames<- as.character(BPossibles)
Wnames<- as.character(WPossibles)

for (w in Wnames){
  Doublet.array[Bnames,w,1]<-BgivenW[,w]
}

#Convert to data frame
Doublet.data <- na.omit(melt(Doublet.array))

##Next Step: work on visualizing


Bhat_doublet_fun(122,50)

Doublet.data122 <- subset(Doublet.data,W==122)


##Calculation of Mean, Variance for K==30 Graph
#mean
meanW <- sum(Doublet.data122$value*Doublet.data122$B)
meanW
#variance
varW<- sum(Doublet.data122$value*(Doublet.data122$B)^2) - meanW^2
varW


