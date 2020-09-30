options(warn=-1)                        #suppress warnings
set.seed(111086)
counter <- 0
## quilt.plot
library(Rcpp)
library(tidyverse)
library(purrr)
##install.packages("INLA", repos=c(getOption("repos"),
##INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)   
library(inlabru)
library(sp)
library(pryr)
library(fields)
library(nloptr)
## require(rgdal, quietly=TRUE)
## scp -r /home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/* ipapasta@xserver2.maths.ed.ac.uk:/home/ipapasta/Software/R/grid_fields/R/
## 
## source R and cpp functions.
## 
source("load_data.R")
## ---------------------------
source("Functions.R")
source("circularHarmonics.R")
## ---------------------------
source("osc_precision.R")
source("hd_precision.R")
## ---------------------------
source("objective.R")
source("priorbetaXW_osc_HD.R")
source("priortheta_osc_HD.R")              
source("gradient_osc_HD.R")
source("hessian_osc_HD.R")
source("llik.R")
source("marginalposterior.R")



##
## load mesh object:
##
## if(!sim)
load("mesh.RData")                     # ! rgdal missing from
p <- mesh$n                         


if(FALSE){
    k <- 1.6
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
}


Ypos.tmp <- data.frame(
    hd=X$hd, speed=X$speed, time=X$synced_time,
    coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(coords.lead = lead(coords)) %>%
    mutate(time.lead = lead(time)) %>%
    mutate(hd.lead = lead(hd)) %>%
    head(-1) #can we do better than this?

## head(Ypos.tmp)
##         hd     speed       time       coords  coords.lead  time.lead  hd.lead
## 1 1.560514  7.602041 0.02207147 94.32308.... 93.99159.... 0.05244587 1.534571
## 2 1.534571 11.343946 0.05244587 93.99159.... 93.83494.... 0.10063787 1.628249
## 3 1.628249  3.304126 0.10063787 93.83494.... 93.57528.... 0.13195947 1.310288
## 4 1.310288  8.673476 0.13195947 93.57528.... 93.29851.... 0.16445867 1.089065
## 5 1.089065 10.162832 0.16445867 93.29851.... 93.32308.... 0.19626667 1.110728
## 6 1.110728  2.668065 0.19626667 93.32308.... 93.48895.... 0.22825387 1.151682

##
## ## split lines and remove lines with length less than a given tolerance.
##

line.segments <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords), filter.zero.length=FALSE,
                             ep=do.call("rbind",Ypos.tmp$coords.lead), tol=.0)

## line.segments.all.lines <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords), filter.zero.length=FALSE,
##                              ep=do.call("rbind",Ypos.tmp$coords.lead), tol=.0)

## line.segments <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords), filter.zero.length = FALSE,
##                              ep=do.call("rbind",Ypos.tmp$coords.lead), tol=0.1)


## > cbind(tail(do.call("rbind",Ypos.tmp$coords)), tail(do.call("rbind",Ypos.tmp$coords.lead)))
##              [,1]     [,2]     [,3]     [,4]
## [47751,] 36.39626 70.55233 35.71713 71.02992
## [47752,] 35.71713 71.02992 34.93738 71.48559
## [47753,] 34.93738 71.48559 34.13147 71.87465
## [47754,] 34.13147 71.87465 33.27284 72.13771
## [47755,] 33.27284 72.13771 32.40027 72.32761
## [47756,] 32.40027 72.32761 31.51980 72.50708


## names(line.segments)
##  "sp"           "ep"           "split.origin" "idx"          "split.loc"
## > cbind(head(line.segments$sp), head(line.segments$ep))
##          [,1]     [,2]     [,3]     [,4]
## [1,] 94.32308 14.49926 93.99159 14.40524
## [2,] 93.57528 14.51369 93.29851 14.69392
## [3,] 93.13883 14.79136 93.35087 14.51937
## [4,] 93.92352 13.83616 94.28017 13.31189
## [5,] 94.40114 13.14841 94.69219 12.75735
## [6,] 94.70665 12.73793 95.03859 12.31383



## df <- data.frame(origin=line.segments$split.origin,
##                  sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
##                  ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
##     group_by(origin) %>%
##     summarize(sp=list(sp), ep=list(ep)) %>%
##     mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
##     mutate(ep = lapply(ep, function(x) do.call("rbind", x)))


## 
## filter.index contains the # of lines used in the integration
## 
df <- data.frame(origin=line.segments$split.origin,
                           filter.index=line.segments$filter.index,
                           sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                           ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
    group_by(origin) %>%
    summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
    mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
    mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 
    ## mutate(filter.index = lapply(filter.index, function(x) do.call("c", filter.index))) 


## ## 
## Ypos <- inner_join(Ypos.tmp %>% mutate(origin=1:nrow(Ypos.tmp)), df) %>%
##     ## length of line segments
##     mutate(Li = map2(sp, ep, function(x, y) apply(y-x, 1, function(z) as.numeric(sqrt(sum(z^2)))))) %>% #!!! 
##     ## times at the nodes of 1d mesh
##     mutate(Ti = pmap(list(time, time.lead, Li), function(x, y, z){
##         o <- as.numeric(x) + (cumsum(as.numeric(z)/sum(as.numeric(z))))*(as.numeric(y) - as.numeric(x))
##         return(matrix(c(x, o), ncol=1))})) %>%
##     mutate(Delta.T = lapply(Ti, function(x) diff(x))) %>% 
##     mutate(s.midpoints = map2(sp, ep, function(x, y) (x+y)/2 )) %>%
##     mutate(t.midpoints = lapply(Ti, function(x){
##         as.vector((na.omit(lag(x)) + na.omit(lead(x)))/2)
##     })) %>%
##     ## circular interpolation for broken lines.
##     ## ----------------------------------------
##     ## mutate(hdi = pmap(list(hd, t.midpoints), function (x, y) rep(x, length(y)) )) %>%
##     mutate(hdi = pmap(list(hd, lead(hd), t.midpoints), function (x, y, z) {        
##         circ.int <- x + (0:(length(z)-1))*(y-x)/length(z)
##         if(is.na(circ.int)) {
##             circ.int <- x
##         }
##         return(circ.int)
##     })) %>%
##     mutate(hdi.diff = pmap(list(hdi), function (x) c(x[-1] - x[-c(length(x))], 0) ))


interpolate <- function(x, y, z){
    interp <- (x) + (cumsum((z)/sum((z))))*((y) - (x))
    Delta <- diff(c(x, interp))
    attr(Delta, "data") <- c(x, interp[-length(interp)]) # note that
                                                         # this way
                                                         # the last
                                                         # value from
                                                         # the
                                                         # combined
                                                         # data vector
                                                         # will be
                                                         # missing so
                                                         # needs to be
                                                         # appended later
    return(Delta)
}

## length of line segments
## duration of time intervals + data as attribute
## change in head direction + data as attribute
Ypos <- inner_join(Ypos.tmp %>%
                   mutate(origin=1:nrow(Ypos.tmp)), df) %>%
    mutate(Li = map2(sp, ep,
                     function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>% #!!! 
    mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
    mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 

filter.index  <- do.call("c", Ypos$filter.index)

## multiplicities for trapezoidal rule
intermediate.knots <- cbind(do.call("rbind",Ypos$sp)[filter.index,1][-1] -
                            do.call("rbind",Ypos$ep)[filter.index,1][-sum(filter.index)])
functions.multiplicity <- intermediate.knots

for(i in 1:length(functions.multiplicity)) {
    if(i==1){
        if(functions.multiplicity[i]!=0) functions.multiplicity[i] <- 1
        else functions.multiplicity[i] <- 2
    }else{
        if(functions.multiplicity[i]!=0) {
            functions.multiplicity[i] <- 1
            functions.multiplicity[i-1] <- 1} else{
                                                functions.multiplicity[i] <- 2
                                            }
    }
}

functions.multiplicity <- c(1, functions.multiplicity, 1)
##  %>%
## mutate(Delta.HD = lapply(hdi, function(x) {
##     o <- diff(x)
##     o[is.na(o)] <- 1e-4 #why does NA occur here?
##     return(o)
## }))
## circ.int <- x + (0:(length(z)-1))*(y-x)/length(z)
## if(is.na(circ.interp)) {
##     circ.interp <- x
## } 

## o <- Ypos.all.lines[do.call("c", Ypos.all.lines$filter.index),]



## nn <- length(o1)
## o1 <- do.call("c", Ypos.all.lines$hdi)[do.call("c", Ypos.all.lines$filter.index)]
## o2 <- do.call("c", Ypos.all.lines$Delta.HD)[do.call("c", Ypos.all.lines$filter.index)]/ do.call("c", Ypos.all.lines$Delta.T)[do.call("c", Ypos.all.lines$filter.index)]
## par(mfrow=c(2,1))
## plot(o1, type="l")
## plot(o2, type="l")

## plot(rank(o1)/(1+nn), rank(o2)/(1+nn))



## dim(do.rbind(Ypos$Li))
## length(do.call("c", Ypos$Li))


## for HD.data and T.data, the last element is missing and needs to be appended
## coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep)[filter.index,],1))
## HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD)[filter.index], )
## T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T)[filter.index]))

## no filter
coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))

## dim(coords.trap)
## length(functions.multiplicity)
      

## A matrices
proj.s  <- inla.mesh.projector(mesh, loc=coords.trap)
s.psi   <- data.frame(proj.s$proj$bary)
names(s.psi) <- c("psi.1", "psi.2", "psi.3")
data <- list(Ypos=Ypos, Y=Y)


## 
## matrices of basis functions
##
## ---------------------------------------------------------------
order.HD <- 5
Aosc.obs  <- inla.spde.make.A(mesh=mesh, loc=as.matrix(data$Y %>% dplyr:: select(position_x, position_y)))
Aosc      <- proj.s$proj$A
Ahd.obs   <- circular.make.A(data$Y$hd, order=order.HD)
Ahd       <- circular.make.A(HD.data, order=order.HD)
## ---------------------------------------------------------------
Aobs <- inla.row.kron(Ahd.obs, Aosc.obs)
A    <- inla.row.kron(Ahd, Aosc)
## L     <- Matrix(do.call("c", data$Ypos$Li))*Matrix(do.call("c", data$Ypos$hdi)) #!! revisit what Li means here and whether Gamma(t) was implemented

## integration 
Delta.T  <- diff(T.data)
Delta.HD <- diff(HD.data)
Delta.s  <- c(diff(do.call("c", Ypos$Li)), 0) #append a 0 value (fix this with lead)         

Delta.time <- c((diff(T.data)[1]/2), diff(T.data, 2)/tail(diff(T.data,1), -1), (tail(diff(T.data),1)/2))
Delta.sHD <- c(sqrt(Delta.s^2 + Delta.HD^2), 0) #LVCF

W     <- Delta.time*Delta.sHD
tW    <- t(W)
DW    <- Diagonal(length(W), W)
## W.old <- Matrix(t(t(L.old)%*%A))



## 
## 
## ## hyperparameters for priors
##
logrhoL    <- log(10)
logsigmaU  <- log(1) #
alpha1     <- .001
alpha2     <- .4
logrho.HDL    <- log(10)
logsigma.HDU  <- log(1) #
alphaHD1      <- .001
alphaHD2      <- 1e-6
hyperpar      <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU), alphaHD1=alphaHD1, alphaHD2=alphaHD2,
                   rho.HDL=exp(logrho.HDL), sigma.HDU=exp(logsigma.HDU))        
par.theta  <- c(log(22), log(1.75), -3, log(pi/1.5), log(.1))
Xest       <-   Xinit    <- Matrix(rep(0, dim(A)[2]), ncol=1)


betaest       <-   betainit <- nrow(Y)/
    sum((Ypos %>% mutate(speed.lead = lead(speed),
                         dt=c(diff(time)[1], diff(time) ))  %>% head(-1) %>%
         mutate(val=dt*((speed + speed.lead)/2)))$val)
gradest       <- NULL
Hessianest    <- NULL


par=par.theta; fn = pthetapc.prop.marg.post_osc_HD; hyperpar=hyperpar; data=data; type="biased"
X=Xinit; beta=betaest; mesh=mesh; A=A; W=W; DW=DW; Aobs=Aobs; tAobs=t(Aobs); tA=t(A); order.HD=order.HD; acc=1e-6
print.verbose=TRUE; control=list(maxit=500)

opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_HD, hyperpar=hyperpar, data=data, type="biased",
                   X=Xinit, beta=betaest, mesh=mesh, A=A, tA=tA, W=W, tW=tW, DW=DW, Aobs=Aobs, tAobs=tAobs, order.HD=order.HD, acc=1e-6,
                   print.verbose=TRUE, control=list(maxit=10000, trace=TRUE))

print(paste("SAVING output"))

save(list=ls(), file="mle_oscillating_HD.RData")


## opt.theta <- cobyla(x0=par.theta, fn = pthetapc.prop.marg.post_osc_HD, hyperpar=hyperpar, data=data, type="biased",
##                    X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, tW=tW, Aobs=Aobs, tAobs=tAobs, order.HD=order.HD, acc=1e-6,
##                    print.verbose=TRUE)


## opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_temp, hyperpar=hyperpar, data=data, 
##                    X=Xinit, Z=Zinit, beta=betaest, mesh=mesh, mesh1d=mesh1d, A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                    W=W, W.input.to.sparseMatrix=W.input.to.sparse.matrix,
##                    acc=1e-6, print.verbose=TRUE, control=list(maxit=500))
