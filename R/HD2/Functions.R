library(Rcpp)
library(tidyverse)
library(purrr)
library(INLA)
library(inlabru)
library(sp)
library(pryr)
library(fields)
library(nloptr)
source("circularHarmonics.R")
source("osc_precision.R")
source("hd_precision.R")
source("objective.R")
source("priorbetaXW_osc_HD.R")
source("priortheta_osc_HD.R")              
source("llik.R")
source("marginalposterior.R")

ldmvnorm <- function(x, m, Q){
    ## ----------------------------------------------------------------
    ## input: x (evaluation point), 
    ##        m (mean vector), 
    ##        Q (precision matrix/inverse covariance Sigma^{-1})
    ## output: log density of multivariate normal distribution
    ## ----------------------------------------------------------------
    k      <- length(x)
    xminmu <- Matrix(x-m,ncol=1)
    L      <- (Matrix::Cholesky(Q, perm=TRUE) %>% (Matrix::expand))$L
    halflogdetQ   <- sum(log(diag(L)))     #!!!
    ## print(paste("halflogdetQ is:", halflogdetQ))
    out    <- halflogdetQ - (k/2)*log(2*pi) - .5*sum((xminmu)* (Q %*% xminmu)) # remove matrix multiplication
    out <- as.numeric(out)
    return(out)
}

## theta <- xpar.theta  <- c(3, -1/2, -3.8, log(7), log(10))
## theta.X <- theta[1:3]
## theta.Z <- theta[4:5]    
## QX     <- osc.precision(theta=theta.X, mesh=mesh)
## QZ     <- temp.precision(theta=theta.Z, mesh=mesh1d) 

## for(i in 1:1000){
##     Sys.sleep(1)
##     ## par(mfrow=c(1,2))
##     ## plot(i1)
##     ## plot(i2)
##     ## i1 <- image(QZ[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2)
##     ## i2 <- image(QX[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2)    
##     plot(image(QZ[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2))
##     Sys.sleep(2)
##     plot(image(QX[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2), add=TRUE)
##     Sys.sleep(4)
## }

## x<-Z
## m<-matrix(rep(0, length(Z)), ncol=1)
## Q<-QZ

split.lines <- function(mesh, sp, ep, filter.zero.length = TRUE, tol= 1e-8, return.filter.index=TRUE) {
    ## locations for splitting
    loc = as.matrix(rbind(sp,ep))
    idx = 1:dim(sp)[1]
    ## Filter out segments not on the mesh
    t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
    t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
    ## if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
    sp = matrix(sp[!((t1==0) | (t2==0)),], ncol=2)
    ep = matrix(ep[!((t1==0) | (t2==0)),], ncol=2)
    idx = idx[!((t1==0) | (t2==0))]
    loc = as.matrix(rbind(sp,ep))
    ## Split them segments into parts
    if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
    np = dim(sp)[1]
    sp.idx = t(rbind(1:np,np+1:np))
    splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
    ##plot(data$mesh)
    ##points(loc)
    ##points(splt$split.loc,col="blue)
    sp = matrix(splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]], ncol=2) # Start point of new segments
    ep = matrix(splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]], ncol=2) # End points of new segments
    idx = idx[splt$split.idx[,1]]
    origin = splt$split.origin
    sl = apply((ep-sp)^2,MARGIN=1,sum)
    filter.index <- !(sl < tol^2)
    ## Filter out zero length segments
    if ( filter.zero.length ) {
        ## sl = apply((ep-sp)^2,MARGIN=1,sum)
        sp = sp[!(sl<tol^2),]
        ep = ep[!(sl<tol^2),]
        origin = origin[!(sl<tol^2)]
        idx = idx[!(sl<tol^2)]
        filter.index = filter.index[!(sl<tol^2)]
        ## drop filter.index
        return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
    }
    return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc, filter.index=filter.index))
}


split.lines <- function(mesh, sp, ep, filter.zero.length = TRUE, tol= 1e-8, return.filter.index=TRUE) {
    ## locations for splitting
    loc = as.matrix(rbind(sp,ep))
    idx = 1:dim(sp)[1]
    ## Filter out segments not on the mesh
    t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
    t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
    ## if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
    sp = matrix(sp[!((t1==0) | (t2==0)),], ncol=2)
    ep = matrix(ep[!((t1==0) | (t2==0)),], ncol=2)
    idx = idx[!((t1==0) | (t2==0))]
    loc = as.matrix(rbind(sp,ep))
    ## Split them segments into parts
    if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
    np = dim(sp)[1]
    sp.idx = t(rbind(1:np,np+1:np))
    splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
    ##plot(data$mesh)
    ##points(loc)
    ##points(splt$split.loc,col="blue)
    sp = matrix(splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]], ncol=2) # Start point of new segments
    ep = matrix(splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]], ncol=2) # End points of new segments
    idx = idx[splt$split.idx[,1]]
    origin = splt$split.origin
    sl = apply((ep-sp)^2,MARGIN=1,sum)
    filter.index <- !(sl < tol^2)
    ## Filter out zero length segments
    if ( filter.zero.length ) {
        ## sl = apply((ep-sp)^2,MARGIN=1,sum)
        sp = sp[!(sl<tol^2),]
        ep = ep[!(sl<tol^2),]
        origin = origin[!(sl<tol^2)]
        idx = idx[!(sl<tol^2)]
        filter.index = filter.index[!(sl<tol^2)]
        ## drop filter.index
        return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
    }
    return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc, filter.index=filter.index))
}



do.rbind <- function(x) do.call("rbind", x)
do.cbind <- function(x) do.call("cbind", x)



## Correlation function of head directional process
cor_HD <- function(theta, kappa, phi)
{
    K <- 1000
    helper1 <-function(theta, mu=0){
        k<- 1:K
        sapply(theta, function(x) sum((1/(kappa^4+2*phi*(kappa^2)*k^2+(k^4))) * (cos(k*x)-mu)))
    }
    helper2 <-function(){
        k<- 1:K
        sum((1/(kappa^4+2*phi*(kappa^2)*k^2+(k^4))))
    }
    BminC <- helper1(theta, mu=1)
    ApluC <- (1/(2*kappa^4)) + helper2()
    return(1+(BminC/ApluC))
}


if(FALSE){
    logkappa<-log(seq(0.2, 2, len=10))
    tau <- 10
    theta0<-0
    theta<-seq(0,2*pi,len=100)
    ## 
    par(mfrow=c(1,1))
    plot(theta, cor_HD(theta, kappa=0.2, phi=-.8), ylim=c(-1,1), type="l")
    ## 
    for(i in 1:length(kappa))
        lines(theta, cor_HD(theta, kappa=kappa[i], phi=-.8), ylim=c(-1,1), lty=1)
    ## 
    par(mfrow=c(1,1))
    plot(theta, cor_HD(theta, kappa=1, phi=-.8), ylim=c(-1,1))
    cor_HD(theta, kappa, phi=-0.001)
}
