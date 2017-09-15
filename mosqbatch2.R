source("mosqfuns2.R")
library(plyr)
## library(doMC)
## registerDoMC(2)
library(foreach)
nsim <- 1000

estres0 <- matrix(NA,nrow=4,ncol=3,
                 dimnames=list(method=c("MM","MLE","doublets","collision"),
                 stat=c("bias","var","MSE")))

getstats <- function(M,B,nsim) {
    estres <- estres0
    Kvec <- replicate(nsim,Kfun(B,M))
    Wvec <- replicate(nsim,Wfun(B,M))
    Cvec <- replicate(nsim,Cfun(B))
    ##
    Ktab <- table(Kvec)/nsim
    Kvec2 <- as.numeric(names(Ktab))

    library(plyr)
    rfit <- laply(Kvec2,fitB2,M=M)

    Bhat_mm <- sapply(Kvec2[Kvec2<M],Bhat_approx,M=M)
    B_mean_mm <- sum(Bhat_mm*Ktab[Kvec2<M])

    Wtab <- table(Wvec)/nsim
    Wvec2 <- as.numeric(names(Wtab))
    Bhat_doublet <- 1 + sqrt(1+4*M*(M-1)/Wvec2)/2

    Ctab <- table(Cvec)/nsim
    Cvec2 <- as.numeric(names(Ctab))
    Bhat_collision <- sapply(Cvec2,Bhat_T)

    estres["MM","bias"] <- bias_mm <- B_mean_mm -B
    estres["MM","var"] <- sum((Bhat_mm-B_mean_mm)^2*Ktab[Kvec2<M])
    estres["MM","MSE"] <- sum((Bhat_mm-B)^2*Ktab[Kvec2<M])

    Bhat_mle <- exp(na.omit(rfit[,1]))
    B_mean_mle <- sum(Bhat_mle*Ktab[Kvec2<M])
    estres["MLE","bias"] <- bias_mle <- B_mean_mle -B
    estres["MLE","var"] <- sum((Bhat_mle-B_mean_mle)^2*Ktab[Kvec2<M])
    estres["MLE","MSE"] <- sum((Bhat_mle-B)^2*Ktab[Kvec2<M])
    B_mean_doublet <- sum(Bhat_doublet[-1]*Wtab[-1])
    estres["doublets","bias"] <- bias_doublet <- B_mean_doublet -B
    estres["doublets","var"] <- sum((Bhat_doublet[-1]-B_mean_doublet)^2*Wtab[-1])
    estres["doublets","MSE"] <- sum((Bhat_doublet[-1]-B)^2*Wtab[-1])
    ##
    B_mean_collision <- sum(Bhat_collision*Ctab)
    estres["collision","bias"] <- bias_collision <- B_mean_collision -B
    estres["collision","var"] <- sum((Bhat_collision-B_mean_collision)^2*Ctab)
    estres["collision","MSE"] <- sum((Bhat_collision-B)^2*Ctab)

    estres
}

Mfrac <- seq(0.1,1,by=0.1)
Bvec <-  round(10^seq(1.5,3,length=7))
rr <- expand.grid(Mfrac=Mfrac,B=Bvec)
rr <- transform(rr,M=round(Mfrac*B))
rr <- subset(rr,M>3 & M<600)
with(rr,plot(M,B))
## checkpointing??
a0 <- aaply(rr,1,
      function(x) {getstats(M=x[["M"]],B=x[["B"]],nsim=nsim)},
            .progress="text") ## ,.parallel=TRUE)

save("a0",file="mosqbatch1.RData")
