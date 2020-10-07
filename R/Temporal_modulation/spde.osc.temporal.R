## 
## seed for reproducibility
##
## set.seed(111086) 
## !! quilt.plot
## 
## load packages
##
counter <- 0
## sim     <- FALSE
library(Rcpp)
library(tidyverse)
library(purrr)
library(INLA)   #install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
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
source("Functions.R")
source("osc_precision.R")
source("objective.R")
source("temp_precision.R")
source("priorbetaXZ_osc_temp.R")
source("priortheta_osc_temp.R")              
source("gradient_osc_temp.R")
source("hessian_osc_temp.R")
source("llik.R")
source("marginalposterior.R")

##
## load mesh object:
##
## if(!sim)
load("mesh.RData")                     # ! rgdal missing from
                                           # maths computing server so
                                           # currently saved locally
                                           # and uploaded for
                                           # practical
                                           # reasons. Request to be
                                        # installed.

if(FALSE){
    k <- 1.6
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
}

## if(sim)
    ## load("simmesh.RData");     load("/home/ipapasta/Desktop/simmesh.RData")

p <- mesh$n
options(warn=-1)                        #suppress warnings
##
## create a regular mesh of points between all line segments defined
## by the trajectory of the mouse (positional data)
##
Ypos.tmp <- data.frame(speed=X$speed, time=X$synced_time, coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(coords.lead = lead(coords)) %>%
    mutate(time.lead = lead(time)) %>% head(-1)
## need to alter code to split lines with inlabru:::split.lines
## line.segments <- inlabru:::split_lines(mesh, sp=do.call("rbind", Ypos.tmp$coords), ep=do.call("rbind",Ypos.tmp$coords.lead))
line.segments <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords), ep=do.call("rbind",Ypos.tmp$coords.lead), tol=.5)
## line.segments.full <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords), ep=do.call("rbind",Ypos.tmp$coords.lead), tol=1e-5)

## ----------------------
## fix issue on tolerance
## ----------------------
df <- data.frame(origin=line.segments$split.origin,
                 sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                 ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
    group_by(origin) %>%
    summarize(sp=list(sp), ep=list(ep)) %>%
    mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
    mutate(ep = lapply(ep, function(x) do.call("rbind", x)))

## df.full <- data.frame(origin=line.segments.full$split.origin,
##                  sp=I(lapply(as.list(apply(line.segments.full$sp, 1, as.list)), unlist)),
##                  ep=I(lapply(as.list(apply(line.segments.full$ep, 1, as.list)), unlist))) %>%
##     group_by(origin) %>%
##     summarize(sp=list(sp), ep=list(ep)) %>%
##     mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
##     mutate(ep = lapply(ep, function(x) do.call("rbind", x)))



o <- inner_join(Ypos.tmp %>% mutate(origin=1:nrow(Ypos.tmp)), df)
## o.full <- inner_join(Ypos.tmp %>% mutate(origin=1:nrow(Ypos.tmp)), df.full)

## Ypos.tmp %>% mutate(origin=1:nrow(Ypos.tmp), df)

Ypos <- o %>%
    ## length of line segments
    mutate(Li = map2(sp, ep, function(x, y) apply(y-x, 1, function(z) as.numeric(sqrt(sum(z^2)))))) %>%
    ## times at the nodes of 1d mesh
    mutate(Ti = pmap(list(time, time.lead, Li), function(x, y, z){
        o <- as.numeric(x) + (cumsum(as.numeric(z)/sum(as.numeric(z))))*(as.numeric(y) - as.numeric(x))
        return(matrix(c(x, o), ncol=1))
    })) %>%
    mutate(s.midpoints = map2(sp, ep, function(x, y) (x+y)/2 )) %>%
    mutate(t.midpoints = lapply(Ti, function(x){
        as.vector((na.omit(lag(x)) + na.omit(lead(x)))/2)
    }))

## Ypos.full <- o.full %>%
##     ## length of line segments
##     mutate(Li = map2(sp, ep, function(x, y) apply(y-x, 1, function(z) as.numeric(sqrt(sum(z^2)))))) %>%
##     ## times at the nodes of 1d mesh
##     mutate(Ti = pmap(list(time, time.lead, Li), function(x, y, z){
##         o <- as.numeric(x) + (cumsum(as.numeric(z)/sum(as.numeric(z))))*(as.numeric(y) - as.numeric(x))
##         return(matrix(c(x, o), ncol=1))
##     })) %>%
##     mutate(s.midpoints = map2(sp, ep, function(x, y) (x+y)/2 )) %>%
##     mutate(t.midpoints = lapply(Ti, function(x){
##         as.vector((na.omit(lag(x)) + na.omit(lead(x)))/2)
##     })) 

## Creation of 1D mesh 
t.nodes <- c(Ypos$time[[1]], do.call("c", (Ypos %>% mutate(t.nodes = pmap(list(time, time.lead, Li), function(x, y, z){
    o <- as.numeric(x) + (cumsum(as.numeric(z)/sum(as.numeric(z))))*(as.numeric(y) - as.numeric(x))
    return(o)
})))$t.nodes))

## t.nodes.full <- c(Ypos.full$time[[1]], do.call("c", (Ypos.full %>%
##                          mutate(t.nodes = pmap(list(time, time.lead, Li), function(x, y, z){
##                              o <- as.numeric(x) + (cumsum(as.numeric(z)/sum(as.numeric(z))))*(as.numeric(y) - as.numeric(x))
##                              return(o)
##                          })))$t.nodes))



## Creation of W matrix
mesh1d  <- inla.mesh.1d(loc=t.nodes, order=2)
proj.t  <- inla.mesh.projector(mesh1d, loc=do.call("c", Ypos$t.midpoints))
proj.s  <- inla.mesh.projector(mesh, loc=do.call("rbind", Ypos$s.midpoints))
s.tv    <- data.frame(mesh$graph$tv[ proj.s$proj$t, , drop=FALSE ])
names(s.tv)  <- c("v1", "v2", "v3")
s.psi  <- data.frame(proj.s$proj$bary)
names(s.psi) <- c("psi.1", "psi.2", "psi.3")

df.W <- cbind(data.frame(line.id=1:nrow(s.psi)),
              data.frame(origin=line.segments$split.origin),
              data.frame(Li=do.call("c", Ypos$Li)),
              s.tv,
              data.frame(knot1=1:length(do.call("c", Ypos$Li)),
                         knot2=2:(1+length(do.call("c", Ypos$Li)))),
              s.psi,
              data.frame(psi.tilde.1 = rep(1/2, length(do.call("c", Ypos$Li)))),
              data.frame(psi.tilde.2 = rep(1/2, length(do.call("c", Ypos$Li))))) %>%
    mutate(jk = pmap(list(knot1, knot2, v1, v2, v3), function(x, y, z, w, k) expand.grid(c(x, y), c(z, w, k)))) %>%
    mutate(x  = pmap(list(psi.tilde.1, psi.tilde.2, psi.1, psi.2, psi.3, Li), function(x, y, z, w, k, L){
        matrix(apply(expand.grid(c(x, y), c(z, w, k)), 1, prod)*L, ncol=1) 
    })) %>%
    mutate(jkx = map2(jk, x, function(x1, x2) {
        df <- cbind(x1, x2)
        names(df) <- c("knot", "v", "value")
        df
    })) %>%
    as_tibble

zero.entries <- data.frame(expand.grid(c(1,2), setdiff(1:p, unique(do.call("rbind", df.W$jkx)$v))),
                          rep(0, dim(expand.grid(c(1,2), setdiff(1:p, unique(do.call("rbind", df.W$jkx)$v))))[1]))
names(zero.entries) <- c("knot", "v", "value")
W.input.to.sparse.matrix <- rbind(do.call("rbind", df.W$jkx), zero.entries) %>%arrange(v)
W  <- sparseMatrix(i=W.input.to.sparse.matrix$knot, j=W.input.to.sparse.matrix$v, x=W.input.to.sparse.matrix$value)
object_size(W) #needs pryr


## W <- W[c(line.segments$split.origin[1], line.segments$split.origin), ]
## rownames(W) <- c(line.segments$split.origin[1], line.segments$split.origin)


## W.input.to.sparse.matrix <- W.input.to.sparse.matrix %>% mutate(bool.filter.rows = knot %in% line.segments$split.origin)
## W.input.to.sparse.matrix.CHECK <- W.input.to.sparse.matrix %>% mutate(bool.filter.rows = knot %in% unique(c(line.segments$split.origin[1], line.segments$split.origin)))

## W <- sparseMatrix(i=W.input.to.sparse.matrix[W.input.to.sparse.matrix.CHECK$bool.filter.rows==TRUE,]$knot,
##                   j=W.input.to.sparse.matrix[W.input.to.sparse.matrix.CHECK$bool.filter.rows==TRUE,]$v,
##                   x=W.input.to.sparse.matrix[W.input.to.sparse.matrix.CHECK$bool.filter.rows==TRUE,]$value)

## mutate(v.ordered = plyr::mapvalues(W.input.to.sparse.matrix$v, from=unique(W.input.to.sparse.matrix$v), to=1:length(unique( W.input.to.sparse.matrix$v)))) %>%
    ## o <- SPL.xyz((mesh$loc[mesh$graph$tv[15,],-3])[1,], (mesh$loc[mesh$graph$tv[15,],-3])[2,], (mesh$loc[mesh$graph$tv[15,],-3])[3,], ID="1")
## plot(mesh, asp=1, xlim=c(0,100), ylim=c(0,100))
## plot(SpatialPolygons(list(o)) , add=TRUE, col=2, lwd=2)

dim(W) 
data <- list(Ypos=Ypos, Y=Y)

## df.W %>% dplyr::select(v1, v2)
## o <- df.W %>% gather(v1, v2, v3, knot1, knot2, -origin)
## o <- df.W %>% dplyr::select(v1, v2, v3, knot1, knot2, psi.1, psi.2, psi.3, psi.tilde.1, psi.tilde.2, Li) %>%
##     gather(value, psi, -v1,-v2,-v3,-knot1,-knot2, -Li) %>% as_tibble

if(FALSE){
    ## points(mycoords, col=2, pch=16, add=TRUE)
    plot(line.segments.xp$sp[,1], line.segments.xp$sp[,2], col=4, pch=16, asp=1)
    points(line.segments.xp$ep[,1], line.segments.xp$ep[,2], col=4, pch=16)
    plot(mesh, asp=1, add=TRUE)
    arrows(do.call("rbind",Ypos$coords)[10:20,1],
           do.call("rbind",Ypos$coords)[10:20,2],
           do.call("rbind",Ypos$lead)[10:20,1],
           do.call("rbind",Ypos$lead)[10:20,2], add=TRUE, lty=3)   
}

##
Aobs  <- inla.spde.make.A(mesh=mesh, loc=as.matrix(data$Y %>% dplyr:: select(position_x, position_y)))
A     <- proj.s$proj$A
Atildeobs  <- inla.spde.make.A(mesh=mesh1d, loc=as.matrix(data$Y$firing_times))
Atilde     <- proj.t$proj$A


## 
## ## hyperparameters for priors
##
logrhoL    <- log(10)
logsigmaU  <- log(1) #
alpha1     <- .001
alpha2     <- 1e-6
logrho.tempL    <- log(10)
logsigma.tempU  <- log(1) #
alphaT1     <- .001
alphaT2     <- 1e-6
hyperpar   <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU),
                   alphaT1=alphaT1, alphaT2=alphaT2, rho.tempL=exp(logrho.tempL), sigma.tempU=exp(logsigma.tempU))        
par.theta  <- c(3, -1/1.5, -3.8, log(.7), log(.3))
Xest       <-   Xinit    <- Matrix(rep(0, mesh$n), ncol=1)
Zest       <-   Zinit    <- Matrix(rep(0, mesh1d$n), ncol=1)

betaest       <-   betainit <- nrow(Y)/sum((Ypos %>% mutate(speed.lead = lead(speed), dt=c(diff(time)[1], diff(time) ))  %>%
                                            head(-1) %>%
                                            mutate(val=dt*((speed + speed.lead)/2)))$val)
                                        #(number of firing events)/
                                        #(\int_\Gamma ds)
gradest       <- NULL
Hessianest    <- NULL



opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_temp, hyperpar=hyperpar, data=data, 
                   X=Xinit, Z=Zinit, beta=betaest, mesh=mesh, mesh1d=mesh1d, A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
                   W=W, W.input.to.sparseMatrix=W.input.to.sparse.matrix,
                   acc=1e-6, print.verbose=TRUE, control=list(maxit=500))

save(list=ls(), file="mle_oscillating_temporal_modulation.RData")


## opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_fixedphirho, hyperpar=hyperpar, data=data, type="biased",
##                        X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, acc=1e-6, print.verbose=FALSE, control=list(maxit=500))


## cobyla
if(FALSE){
    opt.theta <- cobyla(x0=par.theta, fn = pthetapc.prop.marg.post_osc, hyperpar=hyperpar, data=data, type="origin",
                        X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, acc=1e-6, lower=c(-Inf, -Inf, -Inf), upper=c(Inf, Inf, Inf))
}












