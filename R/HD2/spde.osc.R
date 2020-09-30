## 
## seed for reproducibility
##
## set.seed(111086) 
## !! quilt.plot
## 
## load packages
##
counter <- 0
sim     <- FALSE
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
source("priorXbeta_osc.R")
source("priortheta_osc.R")              
source("gradient_osc.R")
source("hessian_osc.R")
source("llik.R")
source("pthetamargipost_osc_stable.R")
sourceCpp("hessian_functions.cpp")

##
## load mesh object:
##
if(!sim)
    load("mesh.RData")                     # ! rgdal missing from
                                           # maths computing server so
                                           # currently saved locally
                                           # and uploaded for
                                           # practical
                                           # reasons. Request to be
                                           # installed.

if(sim)
    load("simmesh.RData");     load("/home/ipapasta/Desktop/simmesh.RData")

p <- mesh$n
options(warn=-1)                        #suppress warnings
##
## create a regular mesh of points between all line segments defined
## by the trajectory of the mouse (positional data)
##
Ypos <- data.frame(speed=X$speed, synced_time=X$synced_time, coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(lead = lead(coords)) %>% head(-1) %>%
    ## compute the unit vector in the direction of the next point
    mutate(dir    = map2(coords, lead, function(x, y) (y-x))) %>% ##/sqrt(sum((y-x)^2))
    ## compute the length of the lines joining the points
    mutate(Li     = map2(coords, lead, function(x, y) as.numeric(sqrt(sum((y-x)^2)))))

##
## set a threshold value to be used as a rule for the segmentation (k = ceiling of |Li|/h)
##
h <- quantile(unlist(Ypos$Li), 0) 


Ypos <- Ypos %>%
    ## the final number of line segments that compose each line in s (if 1 then line is not segmented)
    mutate(k   = ceiling(as.numeric(Li)/h)) %>%
    ## and the length of the line segments is given by sLi
    mutate(sLi = map2(Li, k, function(Li, k) Li/k)) %>%
    mutate(endpoints= pmap(list(coords, dir, k), function(coords, dir, k){
        ## -------------------------------------------------------------
        ## compute the coordinates of the endpoints of all line segments
        ## if k==1 (1 line segment),
        ## Test: code must return the actual coordinates of the
        ## initial line segment. DONE
        ## -------------------------------------------------------------
        endpoints <- coords # matrix(nrow=k+1, ncol=2)
        for(i in 1:k) endpoints <- rbind(endpoints, coords + i*(dir/k))
        return(unname(endpoints))
    })) %>% as_tibble %>%
    mutate(midpoints = lapply(endpoints, function(endpoints){
        out <- matrix(apply(endpoints, 2,
                            function(col.vec) {
                                out <- NULL
                                for(i in 1:(length(col.vec)-1)) out[i] <- (col.vec[i]+col.vec[i+1])/2
                                return(as.numeric(out))
                            }), nrow=(nrow(endpoints)-1))
        return(out)
    })) %>%
    mutate(midpoints.weights = map2(midpoints, sLi, function(midpoints, sLi){
        out <- matrix(rep(sLi, nrow(midpoints)), ncol=1)
        return(out)
    }))


if(!sim){
    data <- list(Ypos=Ypos, Y=Y)
}

if(sim){
    load("Ysim.RData")
    data <- list(Ypos=Ypos, Y=Ysim)
}

## -------------------------------------------------------
## Mesh creation:
## We need some reference about the study region, which can
## be provided by the location points or just a domain.
## The location, supplied on the loc argument, are used
## as initial triangulation nodes. A single polygon can
## be supplied to determine the domain extent on the loc.domain
## argument.  If we supply the point locations,
## or the domain is supplied using the loc.domain
## argument, the algorithm find a convex hull mesh.
## A non convex hull mesh can be made when we provide a
## (list of) set of polygons on the boundary argument,
## where each element of this list is of inla.mesh.segment()
## class.  So, one of these three options is mandatory.
## The other mandatory argument is the max.edge.
## This argument specifies the maximum allowed
## triangle edge lengths in the inner domain and in the
## outer extension.  So, it is a scalar or length two vector.
## This argument is numeric on the SAME SCALE UNIT
## as the coordinates.
## A good mesh needs to have triangles as regular
## as possible in size and shape. To help this requirement
## in addition to max.edge , we have the min.angle argument,
## which can  be scalar or length two vector, to specify
## the minimum internal angles of the triangles on
## the inner domain and on the outer extension.
if(FALSE)
{
    k <- 1.
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
    ##         k <- 1.5
    ## mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03,120), cutoff=k/2)
    ## 
    k <- 1.
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
    ## 
    par(mfrow=c(1,2))
    plot(mesh, asp=1)
    plot(mycoords, add=TRUE, col=2, pch=16)
    ## 
    plot(trajectory)
    points(mycoords, pch=16, cex=0.7, col=2)
    ## plot(trajectory, add=TRUE)
    save(mesh, file="mesh.RData")
}
## We can also build boundaries using the
## inla.nonconvex.hull() function.
## args(inla.nonconvex.hull)
## function (points, convex = -0.15,
##           concave = convex,
##           resolution = 40,
##           eps = NULL, crs = NULL)
## NULL
## In this function we provide the points and set
## some constraint.  We can control the shape of
## the boundary including its convexity, concavity
## and resolution. Here, we make
## some boundaries and build a mesh with each one
## to better understand it.
## Advice from Finn Lindgren: instead of setting
## k to an integer for the creation of the mesh,
## use an ACTUAL PHYSICAL PARAMETER. Note from
## the comments above, the max.edge parameter
## is numeric on the SAME SCALE UNIT as the
## the coordinates.
## 
## ## 
##
Aobs  <- inla.spde.make.A(mesh=mesh, loc=data$Y %>% dplyr:: select(position_x, position_y) %>% as.matrix)
A     <- inla.spde.make.A(mesh=mesh, loc=as.matrix(do.call("rbind",data$Ypos$midpoints)))
L     <- Matrix(do.call("rbind", data$Ypos$midpoints.weights))
W     <- Matrix(t(t(L)%*%A))

## 
## 
## ## hyperparameters for priors
##
logrhoL    <- log(10)
logsigmaU  <- log(1) #
alpha1     <- .001
alpha2     <- 1e-6
hyperpar   <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU))        
par.theta  <- c(3, -1/2, -3.8)
Xest       <-   Xinit    <- Matrix(rep(0, mesh$n), ncol=1)

betaest       <-   betainit <- nrow(Y)/sum((Ypos %>% mutate(speed.lead = lead(speed), dt=c(diff(synced_time)[1], diff(synced_time) ))  %>%
                                            head(-1) %>% mutate(val=dt*((speed + speed.lead)/2)))$val)
                                        #(number of firing events)/
                                        #(\int_\Gamma ds)
gradest       <- NULL
Hessianest    <- NULL

type <- "biased"
if(type=="biased"){
    opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc, hyperpar=hyperpar, data=data, type="biased",
                       X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, acc=1e-6, print.verbose=FALSE, control=list(maxit=500))
    save(list=ls(), file="mle_oscillating_biased.RData")
}else{
    opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc, hyperpar=hyperpar, data=data, type="origin",
                       X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, acc=1e-6, print.verbose=FALSE, control=list(maxit=500))
    save(list=ls(), file="mle_oscillating_origin.RData")
}

## opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_fixedphirho, hyperpar=hyperpar, data=data, type="biased",
##                        X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, acc=1e-6, print.verbose=FALSE, control=list(maxit=500))


## cobyla
if(FALSE){
    opt.theta <- cobyla(x0=par.theta, fn = pthetapc.prop.marg.post_osc, hyperpar=hyperpar, data=data, type="origin",
                        X=Xinit, beta=betaest, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, acc=1e-6, lower=c(-Inf, -Inf, -Inf), upper=c(Inf, Inf, Inf))
}












