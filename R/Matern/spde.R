## 
## set.seed for reproducibility
##
set.seed(111086)

## 
## load packages
##
counter <- 0
library(Rcpp)
library(tidyverse)
library(purrr)
library(INLA)   #install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(inlabru)
library(sp)
library(pryr)
library(fields)
## require(rgdal, quietly=TRUE)

## 
## source R and cpp functions.
## 
source("load_data.R")
source("Functions.R")
source("priorXbeta.R")
source("priortheta.R")
source("gradient.R")
source("hessian.R")
source("llik.R")
source("pthetamargipost.R")
source("pthetapcmargipost.R")
sourceCpp("hessian_functions.cpp")

##
## load mesh object:
## 
load("mesh.RData")                      # ! rgdal missing from maths
                                        # computing server so
                                        # currently saved locally and
                                        # uploaded for practical
                                        # reasons. Request to be
                                        # installed.
p <- mesh$n
options(warn=-1)                        #suppress warnings



## 
## INITIAL PARAMETERS
##
## 2

## args <- commandArgs(trailingOnly = TRUE) !! Rscript
## if((args %>% length)!= 4)
## {
##     paste("Incorrect number of parameters passed by user: check alpha (smoothness), max.edge size (triangulation), range (effective correlation) and scale (variance of Gaussian process) are appropriately given.")    
## }

## 
## Matern shape parameter: range from 3/2 to 2
## 
alpha <- 2   # smooth spatial function
## alpha <- args[1] %>% as.numeric      !! Rscript
## ------------------------------------------------------
## no of parameters in log-linear temporal predictor
## ------------------------------------------------------
## p     <- 1
## p  <- args[2]     !! Rscript

##
## Size of mesh
##
k     <- 10
## k     <- args[2] %>% as.numeric      !! Rscript

##
## Reading and preparing data.
##


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
h <- quantile(unlist(Ypos$Li), 0.4)

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

data <- list(Ypos=Ypos, Y=Y)



## !! revisit selection of hyperparameters.
## 
## log.par.theta <- log(c(1.,.34))         
##


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

if(FALSE)
{
    ## the mesh loaded in line 36 of
    ## this script is identical to
    ## the mesh generated by the
    ## below code. Due to missing
    ## libraries in maths server,
    ## mesh is loaded from an .RData
    ## object.
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 10*k), cutoff=k/2)
    plot(mesh, asp=1)
    plot(mycoords, add=TRUE, col=2, pch=16)
    ## plot(trajectory, add=TRUE)
    save(mesh, file="mesh.RData")
}

Aobs  <- inla.spde.make.A(mesh=mesh, loc=Y %>% dplyr:: select(position_x, position_y) %>% as.matrix)
A     <- inla.spde.make.A(mesh=mesh, loc=as.matrix(do.call("rbind",Ypos$midpoints)))
L     <- Matrix(do.call("rbind", Ypos$midpoints.weights))

if(FALSE)
{
    ## -------------------------------
    ## Oscillating covariance function
    ## -------------------------------
    ## inla.mesh.fem
    
}

## PARAMETER INITIALIZATION AND ACCURACY
log.par.theta <-   c(log(5), log(.1))
Xest          <-   Xinit    <- Matrix(rep(0, mesh$n), ncol=1)
betaest       <-   betainit <-
    nrow(Y)/sum((Ypos %>%
    mutate(speed.lead = lead(speed), dt=c(diff(synced_time)[1], diff(synced_time) ))  %>% head(-1) %>%
    mutate(val=dt*((speed + speed.lead)/2)))$val) #(number of firing events)/ (\int_\Gamma ds)

## load("initialXbeta.RData")
## Xinit       <- Matrix(o$Xest, ncol=1)
## betainit    <- o$betaest 
gradest     <- NULL
Hessianest  <- NULL
## acc     <- 1e-6

## --------------------------------------------------

## 
## prior parameters
##
logrhoL    <- log(20)
logsigmaU  <- log(5)
alpha1     <- .4
alpha2     <- .01
hyperpar   <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU))        


spde  <- inla.spde2.pcmatern(mesh, alpha=alpha, prior.range=c(exp(logrhoL), alpha1), prior.sigma=c(exp(logsigmaU), alpha2)) 


if(FALSE)
{    
    spde = inla.spde2.pcmatern(mesh, prior.range = c(exp(logrhoL), alpha1), prior.sigma = c(exp(logsigmaU), alpha2))
    ## 
    QQ = inla.spde2.precision(spde, theta = c(log(10),log(26)))
    ## - theta: log range and log sigma (standard deviation parameter)
    sdd = sqrt(diag(inla.qinv(QQ)))
    local.plot.field(sdd, mesh)
    points(df$locx, df$locy)
}


## 
## penalised complexity prior
## 
opt.theta <- optim(par=log.par.theta, fn = pthetapc.prop.marg.post, spde=spde, data=data, X=Xinit, beta=betaest, hyperpar=hyperpar,
                   mesh=mesh, L=L, A=A, Aobs=Aobs, acc=1e-7, control=list(trace=TRUE, maxit=1000), method="Nelder-Mead")

opt.theta <- optim(par=c(2.3, 3.25), fn = pthetapc.prop.marg.post, spde=spde, data=data, X=Xinit, beta=betaest, hyperpar=hyperpar,
                   mesh=mesh, L=L, A=A, Aobs=Aobs, acc=1e-7, control=list(trace=TRUE, maxit=5), method="Nelder-Mead")


if(FALSE){
    opt.theta <- optim(par=log.par.theta, fn = ptheta.prop.marg.post, spde=spde, data=data, X=Xinit, beta=betainit,
                       mesh=mesh, L=L, A=A, Aobs=Aobs, acc=1e-7, control=list(trace=TRUE, maxit=5), method="BFGS")
    opt.theta <- optim(par=log.par.theta, fn = ptheta.prop.marg.post2, spde=spde, data=data, X=Xinit, beta=betainit,
                       mesh=mesh, L=L, A=A, Aobs=Aobs, acc=1e-7, control=list(trace=TRUE, maxit=2), method="Nelder-Mead")
    opt.theta <- optim(par=log.par.theta, fn = ptheta.prop.marg.post2, spde=spde, data=data, X=Xinit, beta=betainit,
                       mesh=mesh, L=L, A=A, Aobs=Aobs, acc=1e-7, control=list(trace=TRUE), method="BFGS")
    opt.theta <- optim(par=opt.theta$par, fn = ptheta.prop.marg.post, spde=spde, data=data, X=Xinit, beta=betainit,
                       mesh=mesh, L=L, A=A, Aobs=Aobs, acc=1e-7, control=list(trace=TRUE, maxit=5), method="Nelder-Mead")
}


if(FALSE)
    {
        save(list=ls(), file="current_intensity.RData")
    }

if(FALSE)
    {
        outXbeta  <- ptheta.prop.marg.post.out(opt.theta$par, spde=spde, data=data, X=Xinit, beta=betainit,
                                               mesh=mesh, L=L, A=A, Aobs=Aobs, acc=acc)
    }

## obj.to.save <- ls()## marg.post.optimizer=opt.theta, A.Aobs.L = list(A, Aobs, L))


## To run from terminal use
## 
## nohup R CMD BATCH spde.R &
## 
## TODO: use Rscript

initstring <- log.par.theta %>% round(2) %>% as.character %>% paste(collapse='_')

file.to.save <- paste("lgcp_spde_fit_", "shape_", alpha, "_mesh_edge_size", round(k,2), "_initial_values_", initstring, ".RData", sep="")

file.to.save <- paste("lgcp_spde_fit_", "shape_", alpha, "_mesh_edge_size_", round(k,2),
                      "_alpha1_", alpha1,"_alpha2_", alpha2, "_rhoL_", exp(logrhoL),"_sigmaU_", exp(logsigmaU),
                      "_initial_values_", initstring, ".RData", sep="")




save(list=ls(), file=file.to.save)

##  par.theta <- opt$par
##  theta.par <- exp(c(-10.257475,   5.253098))



