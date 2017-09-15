## new mosquito-sampling functions

## simulate single sample: number of times each bird is bitten
rfun <- function(B,M) { c(rmultinom(1,size=M,prob=rep(1,B))) }

## simulate sample + tabulate
sfun <- function(B,M,collapse=TRUE) {
    tt <- table(factor(rfun(B,M),levels=0:M))
    if (!collapse) tt else paste(tt,collapse=".")
}

## multinomial coefficient
mchoose <- function(n,log=FALSE) {
    m <- lfactorial(sum(n))-sum(lfactorial(n))
    if (log) m else exp(m)
}

## Maxwell-Boltzmann occupancy probability
occprob <- function(n,B,M,log=FALSE,add.zero=TRUE) {
    if (add.zero) nx <- c(B-sum(n),n) else {
        nx <- n; n <- n[-1]
    }
    r <- -M*log(B)+mchoose(nx,log=TRUE)+
        lfactorial(M)-sum(n*lfactorial(seq_along(n)))
    if (log) r else exp(r)
}

## NLL(K)
nllfun <- function(logB,K,M) {
    if (any(logB-log(K) > 25)) warning("results may be numerically unstable")
    M*logB-lfactorial(exp(logB))+lfactorial(exp(logB)-K)
}

## MLE fit (see fitB2 for improved version)
fitB <- function(K,M) {
    require(bbmle)
    m1 <- mle2(nllfun,optimizer="optimize",
               lower=log(M),upper=25,
               start=list(logB=10),
               data=list(K=K,M=M))
    cc <- confint(m1,quietly=TRUE)
    list(fit=coef(m1),confint=cc)
}

## wrapper for fitB (see fitB2 for improved version)
ffun <- function(K,M) {
    if (K<M) {
        f <- try(fitB(K,M),silent=TRUE)
        if (inherits(f,"try-error")) rep(NA,3) else unlist(f)
    } else {
        c(NA,lboundfun(M),NA)
    }
}

Kfun <- function(B,M) {
    ## simulate K
    unname(B-sfun(B,M,collapse=FALSE)[1])
}

## same as occprob(B-M,B,M ... ?)
collision_prob <- function(M,B,inverse=FALSE,log=FALSE) {
    r <- sum(log(B-seq(0,M-1)))-M*log(B)
    if (inverse) {
        if (log) r else exp(r)
    } else {
        r2 <- 1-exp(r)
        if (log) log(r2) else r2
    }
}

## method of moments estimator for K
Bhat_approx <- function(K,M,type=c("frac","exp"),upper=1e7) {
    type <- match.arg(type)
    if (K>=M) stop("K must be < M")
    tmpf <- if (type=="exp") function(B) B*(1-exp(-M/B))-K else
    function(B) B*(1-(1-(1/B))^M)-K
    r <- try(uniroot(tmpf,interval=c(K+0.001,upper))$root)
    if (inherits(r,"try-error")) NA else r
}

## expected time to first collision
Texpfun <- function(B,maxn=4*B) {
    nvec <- 2:maxn
    sum(exp(-(nvec*(nvec-1))/(2*B)))
}

## root-finding solution of B
Bhat_T <- function(tau) {
    r <- try(uniroot(function(x) { tt <- Texpfun(x)
                                   ## cat(tau,tt,"\n")
                                   tt-tau },
                     interval=c(tau,1e6))$root,silent=TRUE)
    if (inherits(r,"try-error")) NA else r
}


## simulate doublets
Wfun <- function(B,M) {
    r <- rfun(B,M)
    mean(r*(r-1))
}

## simulate time to first collision
Cfun <- function(B,ssize=2*B) {
    ss <- sample(B,size=ssize,replace=TRUE)
    min(which(duplicated(ss)))
}

## derivative of lfactorial(exp(x)-c)
grad_exp_lfactorial <- function(x,c=0) {
    digamma(exp(x)-c+1)*exp(x)
}

## derivative of NLL(K)
dnllfun <- function(logB,K,M) {
    M-grad_exp_lfactorial(logB)+grad_exp_lfactorial(logB,K)
}

## revised MLE: uses root-finding on d(NLL) and NLL directly, instead of
##  going through mle2 [should be faster & more robust]
fitB2 <- function(K,M,maxval=15,alpha=0.05) {
    nllfun.off <- function(x,K,M,val) {
        nllfun(x,K,M)-(val+qchisq(1-alpha,1)/2)
    }
    u1 <- try(uniroot(dnllfun,K=K,M=M,interval=c(log(K)+0.001,maxval))$root,silent=TRUE)
    if (inherits(u1,"try-error")) return(rep(NA,3))
    val <- nllfun(u1,K,M)
    u.lo <- try(uniroot(nllfun.off,K=K,M=M,val=val,interval=c(log(K)+0.0001,u1))$root,silent=TRUE)
    if (inherits(u.lo,"try-error")) u.lo <- NA
    u.hi <- try(uniroot(nllfun.off,K=K,M=M,val=val,interval=c(u1,maxval))$root,silent=TRUE)
    if (inherits(u.hi,"try-error")) u.hi <- NA
    c(u1,u.lo,u.hi)
}
