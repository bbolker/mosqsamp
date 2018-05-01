#New functions Heterogeneous Mosquito sampling

source("mosqfuns2.R")


library(fossil)

#===============================================================================
#VHetfunction for drawing a multinomial sample with heterogeneous
#probabilities

vHetfun <- function(M,prob) { 
  c(rmultinom(1,size=M,prob=prob))
}

sHetfun <- function(M,prob){
  table(vHetfun(M,prob),dnn = "Occupancy spectrum including S_0")
}

#===============================================================================
#SHetfunction2 for drawing a multinomial sample with heterogeneous
#probabilities, corresponding function to vHetfun, special format- creating a matrix
#of values, 
#made for use with ABCfun functions

sHetfun2 <- function(M,prob){
  s<- matrix(data=0,nrow = 1, ncol=length(0:M),dimnames=list(TrueOccup= 1,Occupancy = 0:M))
  s_fill <- sHetfun(M,prob)
  s[1,names(s_fill)]<-s_fill
  s
}


## One realization of Chao using Homogenous probability
#===============================================================================


Chaofun <- function (B,M){
  v <- vfun(B,M)
  chao1(v)
}

## One realization of Chao using inputed probability
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