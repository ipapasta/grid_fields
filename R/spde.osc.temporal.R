## 
## seed for reproducibility
##
## set.seed(111086) 
## !! quilt.plot
## sim     <- FALSE
## library(Rcpp)
library(tidyverse)
library(purrr)
library(INLA)   #install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
inla.setOption(pardiso.license = "/Users/ipapasta/pardiso.lic")
library(inlabru)
library(sp)
library(fields)
library(nloptr)
library(pals)
## require(rgdal, quietly=TRUE)
## scp -r /home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/* ipapasta@xserver2.maths.ed.ac.uk:/home/ipapasta/Software/R/grid_fields/R/
## 
## source R and cpp functions.
## 
source("load_data.R")
source("Functions.R")
source("osc_precision.R")
source("hd_precision.R")
## source("objective.R")
source("temp_precision.R")
## source("priorbetaXZ_osc_temp.R")
## source("priortheta_osc_temp.R")              
## source("gradient_osc_temp.R")
## source("hessian_osc_temp.R")
## source("llik.R")
## source("marginalposterior.R")

k    <- 5
mesh      <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
## 
fem.mesh  <- inla.mesh.fem(mesh, order = 2)
## Bspatial.phi0 = matrix(c(0,1,0,0), nrow=1)
## Bspatial.phi1 = matrix(c(0,0,1,0), nrow=1)
## Bspatial.phi2 = matrix(c(0,0,0,1), nrow=1)
## M0.spatial = fem.mesh$c0 # C
## M1.spatial = fem.mesh$g1
## M2.spatial = fem.mesh$g2
## 
p         <- mesh$n                         
theta.nodes <- seq(0, 2*pi, len=30)
mesh.hd     <- inla.mesh.1d(theta.nodes, boundary="cyclic", degree=1)
fem.mesh.hd <- inla.mesh.fem(mesh.hd, order = 2)
##
## Bhd.phi0 = matrix(c(0,1,0), nrow=1)
## Bhd.phi1 = matrix(c(0,0,1), nrow=1)
## M0.hd = fem.mesh.hd$c0
## M1.hd = fem.mesh.hd$g1
## M2.hd = fem.mesh.hd$g2
## 
nodes       <- c(mesh.hd$loc, 2*pi)
intervals   <- head(cbind(nodes, lead(nodes)), -1)

Ypos.tmp <- data.frame(
    hd=X$hd, time=X$synced_time,
    coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(coords.lead = lead(coords)) %>%
    mutate(time.lead = lead(X$synced_time)) %>%
    mutate(hd.lead = lead(X$hd)) %>%
    head(-1)

## options(warn=2)


## o <- lapply(1:nrow(Ypos.tmp), function(k) Ypos.tmp$hd == Ypos.tmp$hd.lead)


Ypos.tmp <- Ypos.tmp %>% mutate(HD.split = map2(hd, hd.lead, split.arcs),
                                L.arcs = lapply(HD.split,
                                                function(x) apply(x, 1,
                                                                  function(y) abs(y[2]-y[1]))),
                                time.split = pmap(list(time, time.lead, L.arcs), function(x,y,z){
                                    o <- interpolate(x,y,z)
                                    oo <- c(attr(o, "data"), y)
                                    ooo <- head(cbind(oo, lead(oo)), -1)
                                    colnames(ooo) <- NULL
                                    return(ooo)
                                }),
                                coords.split=pmap(list(coords, coords.lead, L.arcs), function(x,y,z){
                                    interpolate2(x, y, z)
                                }),
                                new.time = lapply(time.split, function(x) x[,1, drop=FALSE]),
                                new.time.lead= lapply(time.split, function(x) x[,2, drop=FALSE]),
                                new.hd = lapply(HD.split, function(x) x[,1, drop=FALSE]),
                                new.hd.lead = lapply(HD.split, function(x) x[,2, drop=FALSE]),
                                new.coords = lapply(coords.split, function(x) x[,1:2, drop=FALSE]),
                                new.coords.lead = lapply(coords.split, function(x) x[,3:4, drop=FALSE])
                                )

Ypos.tmp <- Ypos.tmp %>% dplyr::select(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
    unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))

names(Ypos.tmp) <- c("time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")
## Ypos.tmp$Time.split[3][[1]]
    

line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
                             filter.zero.length=FALSE,
                             ep=Ypos.tmp$coords.lead, tol=.0)

## line.segments.filtered <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords),
##                              filter.zero.length=TRUE,
##                              ep=do.call("rbind",Ypos.tmp$coords.lead), tol=1e-4)
## dim(line.segments.filtered$sp)
## filter.index contains the # of lines used in the integration

df <- data.frame(origin=line.segments$split.origin,
                 filter.index=line.segments$filter.index,
                 sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                 ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
    group_by(origin) %>%
    summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
    mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
    mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 




## length of line segments
## duration of time intervals + data as attribute
## change in head direction + data as attribute


## oo <- inner_join(Ypos.tmp %>%
##                    mutate(origin=1:nrow(Ypos.tmp)), df)


Ypos <- inner_join(Ypos.tmp %>%
                   mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
    mutate(Li = map2(sp, ep,
                     function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
    mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
    mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 



filter.index  <- do.call("c", Ypos$filter.index)



## ------------------------------
## integration weights are T.data
## ------------------------------
coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))



## ---------------------------------------
## SpatialPointsDataFrame and SpatialLines
## ---------------------------------------
Y.spdf    <- SpatialPointsDataFrame(SpatialPoints(cbind(Y$position_x, Y$position_y)), as.data.frame(Y%>%dplyr::select(-c(position_x, position_y))))
Ypos.sldf <- SpatialLinesDataFrame(SpatialLines(lapply(as.list(1:nrow(Ypos)),function(k) Lines(list(Line(cbind(c(Ypos$coords[k,1],
                                                                                    Ypos$coords.lead[k,1]),
                                                                                  c(Ypos$coords[k,2],
                                                                                    Ypos$coords.lead[k,2])))),
                                                                                  ID=k))),
                                   Ypos %>% dplyr::select(-c(coords, coords.lead)))

data <- list(Ypos=Ypos, Y=Y, Yspdf=Y.spdf, Ypos.sldf = Ypos.sldf)




## SpatialPointsDataFrame and SpatialLines



## dim(coords.trap)
## length(functions.multiplicity)
      



## circular domain and temporal domain meshes
mesh1d  <- inla.mesh.1d(loc=c(T.data[seq(1, length(T.data), by = 300)], T.data[length(T.data)]), order=2)
print(paste("mesh1d:", head(diff(mesh1d$loc))))

## POSITIONAL - Used for INTEGRAL
Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
A      <- inla.row.kron(Ahd, Aosc)


## OBSERVED
Aosc.obs  <- inla.spde.make.A(mesh=mesh,
                              loc=as.matrix(data$Y %>%
                                            dplyr:: select(position_x, position_y)))
Ahd.obs   <- inla.spde.make.A(mesh=mesh.hd, data$Y$hd)

Aobs      <- inla.row.kron(Ahd.obs, Aosc.obs)
Atildeobs    <- inla.spde.make.A(mesh=mesh1d, data$Y$firing_times)

## c(3.1305657, -0.8904750, -3.3000715,  1.2177724,  1.1621642,  3.0267431, -0.2246174)




## nrow(A)==nrow(At); 92264 each row of above matrices contains
## non-zero values at knots wrapping a distinct line segment.
## Ck <- sapply(dGamma, function(x) rep(x, 6))

## (coords.trap, HD.data, T.data)

dGamma <- c(do.call("c", Ypos$Li))
dT  <- diff(T.data)
Atmp <- as(A, "dgTMatrix")
A.indices <- cbind(cbind(Atmp@i+1, Atmp@j+1), Atmp@x)
A.indices <- A.indices[order(A.indices[,1]),] %>% as.data.frame
Attmp <- as(Atilde, "dgTMatrix")
At.indices <- cbind(cbind(Attmp@i+1, Attmp@j+1), Attmp@x)
At.indices <- At.indices[order(At.indices[,1]),] %>% as.data.frame
names(A.indices) <- c("tk", "i", "psi.ot") #ot: omega x theta
names(At.indices) <- c("tk", "l", "psi.t")


df.unnest <- full_join(At.indices %>% group_by(tk) %>% nest(),
                A.indices %>% group_by(tk) %>% nest(), by="tk") %>%
    arrange(tk) %>%
    ungroup %>% 
    mutate(
        time = T.data,
        time.lag = c(0, time[-length(time)]), #issue with NAs, use 0 (will be discarded later)
        direction     = HD.data,
        direction.lag = c(0, HD.data[-length(direction)]),
        coords = I(coords.trap),
        coords.lag = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
        dGamma=c(dGamma,0),
        dGamma.lead = lead(dGamma),
        dGamma.lag = lag(dGamma),
        val = pmap(list(data.x, data.y, dGamma), function(x, y, z){
            oo <- expand.grid(1:nrow(x), 1:nrow(y))
            ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
            oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val=ooo)
            oooo
        }))

df <- df.unnest %>% dplyr::select(-c("data.x", "data.y")) %>%
    unnest(val)


df.W <- rbind(df %>% mutate(group=tk,
                            dGamma.lag=0) %>%
              dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val),
              df %>% 
              filter(tk!=1) %>%
              mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                     group=tk-1,
                     dGamma=0) %>%
              dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val)) %>%
    arrange(group) %>%
    mutate(dGamma.trap = dGamma + dGamma.lag) ## %>%
    ## select(-c(dGamma, dGamma.lead, dGamma.lag))

## print(df.W %>% select(tk,time, group, dGamma, dGamma.lag, dGamma.trap), n=60)

tol <- 0
df.dGamma.sum.k.kplus1 <- df.W %>% group_by(group, l, i) %>%
    summarize(val = sum(max(dGamma.trap*val, tol)))

df.dGamma.sum.k.kplus1 <- df.W %>% group_by(group, l, i) %>%
    summarize(val = sum(max(dGamma.trap*val, tol)),
              time = unique(time),
              direction=unique(direction),
              coords=unique(coords))

## df.dGamma.sum.k.kplus1$val[df.dGamma.sum.k.kplus1$val==0] <- 1e-16
## the integration scheme is stable, which can be achieved by ensuring
## positive weights on all basis functions that interact with the
## line/curve of integration.



## Include arguments in functions

df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)), space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
mapindex2space.direction_index <- function(index){    
    f<-function(index.single){
        as.numeric(df.indices[which(df.indices$cross==index.single),c("dir","space")])
    }
    t((Vectorize(f, vectorize.args="index.single"))(index))
}

mapindex2space.direction_basis <- function(index){    
    f<-function(index.single){
        o <- as.numeric(df.indices[which(df.indices$cross==index.single),c("dir","space")])
        return(c(mesh.hd$loc[o[1]], mesh$loc[o[2],-3]))
    }
    t((Vectorize(f, vectorize.args="index.single"))(index))
}

W <- sparseMatrix(i=df.dGamma.sum.k.kplus1$l,
                  j=df.dGamma.sum.k.kplus1$i,
                  x=df.dGamma.sum.k.kplus1$val/2)
W         <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(A)-ncol(W)))
W.ipoints.M2 <- as(W, "dgTMatrix")
W.ipoints.M2 <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2@i+1], hd=mapindex2space.direction_basis(W.ipoints.M2@j+1)[,1],
                        coords.x1 =mapindex2space.direction_basis(W.ipoints.M2@j+1)[,2],
                        coords.x2 =mapindex2space.direction_basis(W.ipoints.M2@j+1)[,3],
                        weight=W.ipoints.M2@x) %>% arrange(firing_times)




## the following B matrices are intended to be used with inla.spde2.generic - see spde2_implementation.pdf
## for definition of matrices and quantities below. spde2_implementation presents methods for general non-stationary models
## see Lindgren et al paper for information on GMRF precision matrices when the latent field is non-stationary
## there are two ways for defining models: one with inla.spde2.generic and the other with inla.rgeneric.define
## inla.spde2.generic provides support for matern models (this includes oscillating models too)
## inla.rgeneric.define allows user to build the model from scratch - this includes priors of hyperparameters
## Ideally, we need all models below to be fit with inla.rgeneric.define.
B.phi0.matern = matrix(c(0,1,0), nrow=1)
B.phi1.matern = matrix(c(0,0,1), nrow=1)
B.phi0.oscillating = matrix(c(0,1,0,0), nrow=1)
B.phi1.oscillating = matrix(c(0,0,1,0), nrow=1)
B.phi2.oscillating = matrix(c(0,0,0,1), nrow=1)

## the following commands implement the finite element method and can
## be used to obtain useful quantities such as the M matrices which
## are also defined in spde2_implementation but can be used both in
## inla.spde2.generic and inla.rgeneric.define
fem.mesh    <- inla.mesh.fem(mesh, order = 2)
fem.mesh.hd <- inla.mesh.fem(mesh.hd, order = 2)
fem.mesh.temporal <- inla.mesh.fem(mesh1d, order = 2)
## M matrices for spatial oscillating model
M0 = fem.mesh$c0
M1 = fem.mesh$g1
M2 = fem.mesh$g2
## M matrices for temporal model
M0.temporal = fem.mesh.temporal$c0
M1.temporal = fem.mesh.temporal$g1
M2.temporal = fem.mesh.temporal$g2
## M matrices for circular/directional model
M0.hd = fem.mesh.hd$c0
M1.hd = fem.mesh.hd$g1
M2.hd = fem.mesh.hd$g2


## for syntax on how to write new models see "git-books" or even better
## vignette("rgeneric", package="INLA")
## if you can't open the vignette then you probably have an older version of INLA.
## Install most recent _development_ version
## Finn suggested we amend the lgcp code to include inla.mode="experimental" (recent feature of INLA)
## 

## source all custom-made built models for inla.rgeneric.define
source("rgeneric_models.R")
## define models
## oscilalting.rgeneric is used for M0
## space.direction.rgeneric is used for M1 and M2
## temporal.rgeneric is used for M1 and M2
oscillating.rgeneric     <- inla.rgeneric.define(oscillating.model, M = list(M0=M0, M1=M1, M2=M2))
circular.rgeneric        <- inla.rgeneric.define(circular1D.model,  M = list(M0=M0.hd, M1=M1.hd, M2=M2.hd))
temporal.rgeneric        <- inla.rgeneric.define(temporal.model,    M=list(M0.temporal=M0.temporal, M1.temporal=M1.temporal, M2.temporal=M2.temporal))
space.direction.rgeneric <- inla.rgeneric.define(space.direction.model, M=list(M0.space=M0, M1.space=M1, M2.space=M2,
                                                                               M0.direction=M0.hd, M1.direction=M1.hd, M2.direction=M2.hd))

## ----------------
## Fitting M0 model
## ----------------
## current implementation below is correct
cmp.oscillating.rgeneric <- coordinates ~ 
    spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
    Intercept
cmp.space <- firing_times ~ spde2(space(firing_times), model=space.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) 

fit.oscillating.rgeneric <- lgcp(cmp.oscillating.rgeneric, data = Y.spdf, samplers = Ypos.sldf,
                                 domain = list(coordinates = mesh), options=list(verbose = TRUE))

## plot estimated intensity
pr.int <- predict(fit.oscillating.rgeneric, pxl, ~ spde2)


ggplot() + gg(pr.int) + gg(mycoords, color="red", size=0.2)+
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar")+
    xlim(0,100)+ ylim(0,100)+   
    coord_equal() + theme_classic()

## ----------------
## Fitting M1 model
## ----------------
## NOTES: the integration points that are supplied are incorrect here but were used to verify that the code runs.
## A correct implementation would need to compute the W.ipoints for M1 correctly-Work in progress

cmp.space.direction <- firing_times ~
    spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
          mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE))))

fit.space.direction <- lgcp(cmp.space.direction, data = as.data.frame(Y.spdf),
                            ips=W.ipoints.M2,
                            domain = list(firing_times = mesh1d),
                            options=list(
                                num.threads=8,
                                verbose = TRUE, bru_max_iter=1))



## ----------------
## Fitting M2 model
## ----------------
## current implementation below is correct
cmp.space.direction.time <- firing_times ~
    spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
          mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE)))) +
    time(firing_times, model=temporal.rgeneric, mapper=bru_mapper(mesh1d, indexed=TRUE)) + Intercept

fit.space.direction.time <- lgcp(cmp.space.direction.time, data = as.data.frame(Y.spdf),
                                 ips=W.ipoints.M2,
                                 domain = list(firing_times = mesh1d),
                                 options=list(
                                     num.threads=8,
                                     verbose = TRUE, bru_max_iter=1))



## -------------------------------------------------------------------------
## inlabru implementation
## -------------------------------------------------------------------------

if(FALSE){
    save(fit.space.direction.time, file="sdt_inlabru.RData")
    ## 
    plot(fit.space.direction.time$marginals.hyperpar[[1]], type="l")
    sigmaLN  <- 1; murho    <- 30
    lines(seq(0,100,len=100), dlnorm(seq(0,100,len=100), sigmaLN, log=FALSE), lty=2)
    ##
    plot(fit.space.direction.time$marginals.hyperpar[[2]], type="l")
    lines(seq(0,10,len=100), dexp(seq(0,10,len=100), 1/2, log = FALSE), lty=2)
    ## 
    plot(fit.space.direction.time$marginals.hyperpar[[3]], type="l")
    lines(seq(-1,1,len=100), prior.phi_osc(seq(-1,1,len=100), a=1, b=20, lg=FALSE), lty=2)
    ## 
    plot(fit.space.direction.time$marginals.hyperpar[[4]], type="l")
    lines(seq(0,10,len=100), dexp(seq(0,10,len=100), 1/(2*pi), log = FALSE), lty=2)
    ## 
    plot(fit.space.direction.time$marginals.hyperpar[[5]], lty=2, type="l")
    lines(seq(0,10,len=100), dexp(seq(0,10,len=100), 1, log = FALSE), lty=2)
    ## 
    plot(fit.space.direction.time$marginals.hyperpar[[6]], type="l")
    lines(seq(0,200,len=100), dexp(seq(0,200,len=100), 1/100, log=TRUE), lty=2)
    ## 
    plot(fit.space.direction.time$marginals.hyperpar[[7]], type="l")
    lines(seq(0,10,len=100), dexp(seq(0,10,len=100), 1/3, log = FALSE), lty=2)


    ##
    slack <- 5
    maxx  <- max(pos.coords[,"x"]+slack)
    maxy  <- max(pos.coords[,"y"]+slack)
    minx  <- min(pos.coords[,"x"]-slack)
    miny  <- min(pos.coords[,"y"]-slack)
    N      <- 100
    coords  <- expand.grid(seq(minx, maxx, len=N), seq(miny, maxy, len=N))
    coords  <- matrix(rep(c(40,45),N*N), ncol=2,byrow=TRUE)
    times   <- rep(10,N*N)
    dir     <- seq(0,2*pi, length=N*N)

    dat.2.predict <- data.frame(firing_times=times, hd=dir, coords.x1=coords[,1], coords.x2=coords[,2])
    pr.int.full <- predict(fit.space.direction.time, dat.2.predict, ~ Intercept + spde2 + time)
    pr.int.full <- predict(fit.space.direction.time, NULL, ~ spde2_eval(coords, dir) + time_eval(times))
    pr.int.full <- predict(fit.space.direction.time, NULL, ~ spde2_eval(coords, dir) + time_eval(times))

    p1 <- ggplot(pr.int.full) + 
        geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
        scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                           labels=paste0(0:5,expression("pi"),"/",3)) +
        scale_y_continuous(breaks=seq(-3.5,3.5, by=0.5), limits=c(-3.5, 3.5))+
        geom_line(aes(x=hd, y=mean)) +
        geom_hline(yintercept=0, colour="grey")+
        coord_polar(start = pi, direction=1) + theme_minimal()

    p2<- ggplot(pr.int.full) + 
        geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
        geom_line(aes(x=hd, y=mean)) +
        geom_hline(yintercept=0, colour="grey")+
        coord_fixed()+
        theme_minimal()

    grid.arrange(p1,p2, nrow=1)
    
    pr.int.full <- predict(fit.space.direction.time, NULL, ~ c(Intercept = Intercept_latent))
    pr.int <- predict(fit.space.direction.time, NULL, ~ spde2_eval(coords, dir) + time_eval(times))
    pr.int <- predict(fit.space.direction.time,NULL, ~ Intercept + spde2(coords))
    p3  <- ggplot(pr.int.full, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                             limits=c(min(pr.int.full$mean),max(pr.int.full$mean)))+
        coord_fixed()+ 
        theme_classic() + theme(legend.text=element_text(size=11))
    ## 
    
    pxl <- pixels(mesh, nx=500, ny=500)    
    pr.int <- predict(fit.space.rgeneric, pxl, ~ spde2)

    pr.int <- predict(fit.space.direction.time, pxl, ~ spde2) 
    
    library(pals)
    ggplot() +
        gg(pr.int) +
        gg(mycoords, color="red", size=0.2)+
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar")+
        xlim(0,100)+
        ylim(0,100)+   
        coord_equal() + theme_classic()


    ## 
    ## rgeneric
    oscillating.rgeneric <- inla.rgeneric.define(oscillating.model, M = list(M0=M0, M1=M1, M2=M2))
    circular.rgeneric    <- inla.rgeneric.define(circular1D.model, M=list(M0=M0.hd, M1=M1.hd, M2=M2.hd))
    ## cmp.oscillating.rgeneric <- hd2 + coordinates ~ circular(hd2, model=circular.rgeneric) +
    ##     mySPDE(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
    cmp.oscillating.rgeneric <- coordinates ~ 
        spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
        Intercept

    ## cmp.oscillating.rgeneric <- firing_times ~ spde2(space(firing_times), model=oscillating.rgeneric, ..., mapper=...) + temporal(firing_times, model="OU??")

    space.rgeneric <- inla.rgeneric.define(oscillating.model,
                                            M=list(M0.space=M0, M1.space=M1, M2.space=M2))
    cmp.space <- firing_times ~ spde2(space(firing_times), model=space.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) ## +
        ## time(firing_times, model="ou")
    ## samplers.space <- cbind(Ypos.sldf@data$time, Ypos.sldf@data$time.lead)
    ## ips  <- ipoints(samplers.space.direction, mesh1d)
    ## 
    fit.space <- lgcp(cmp.space, data = Y.spdf@data, domain = list(firing_times = mesh1d),                      
                      options=list(verbose = TRUE, bru_max_iter=1))

    
    ## space takes a vector of firing times and outputs a two-column matrix corresponding to those times.
    ## space(firing_times, ...) check ... or fetch data from global environment.
    
    ## spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
    ## Intercept
    
    ## +        direction(Y.spdf$hd)
    ## hd 
    ## form <- coordinates  ~ spde2 + hd + Intercept 
    fit.oscillating.rgeneric <- lgcp(cmp.oscillating.rgeneric, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh),
                                     ## formula=form,
                                     options=list(verbose = TRUE))
    ## ------------------------------------------------------- if you
    ## give ips then you don't need samplers and  domain stands
    ## for integration domain - used to construct integration
    ## scheme. First take domain and samplers (subintervals).
    ## By default it would put integral
    ## -------------------------------------------------------
    ## ------------------------------------- temporal
    ## -------------------------------------
    cmp.temporal <- firing_times ~ 
        temporal(firing_times, model = "ar1", mapper=bru_mapper(mesh1d, indexed=TRUE)) +
        + Intercept
    fit.temporal <- lgcp(cmp.temporal, data = Y.spdf, samplers = matrix(c(0, max(Y$firing_times)), nrow=1), domain = list(firing_times = mesh1d),
                                     ## formula=form,
                         options=list(verbose = TRUE))
    ## samplers a matrix for temporal process (each row codes an interval in time )
    ## if formula is coordinates ~ . (tells that predictor is linear)
    ## -------------------------------------
    ## directional
    ## -------------------------------------
    circular.rgeneric <- inla.rgeneric.define(circular1D.model, M=list(M0=M0.hd, M1=M1.hd, M2=M2.hd))
    cmp.circular <- hd ~ 
        directional(hd, model = circular.rgeneric, mapper=bru_mapper(mesh.hd, indexed=TRUE)) +
        + Intercept

    fit.circular <- lgcp(cmp.circular, data = Y.spdf, samplers = Ypos.sldf, domain = list(hd = mesh.hd),
                         options=list(verbose = TRUE))

    
    ## -------------------------------------
    ## oscillating and temporal - incorrect
    ## -------------------------------------
    cmp.oscillating.temporal <- coordinates + firing_times ~
        spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
        f(firing_times, model = "ar1", mapper=bru_mapper(mesh1d, indexed=TRUE)) +
        + Intercept
    fit.oscillating.temporal <- lgcp(cmp.oscillating.temporal, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh, firing_times = mesh1d),
                                     ## formula=form,
                                     options=list(verbose = TRUE))    ## 
    ## 
    
    plot(fit.oscillating.rgeneric$marginals.hyperpar[[1]], type="l", xlim=c(-7,10))
    lines(fit.oscillating.rgeneric$marginals.hyperpar[[2]], lty=2)
    lines(fit.oscillating.rgeneric$marginals.hyperpar[[3]], lty=2)
    
    pxl <- pixels(mesh, nx=500, ny=500)
    pr.int <- predict(fit.oscillating.rgeneric, pxl, ~ spde2)

    library(pals)
    ggplot() +
        gg(pr.int) +
        gg(mycoords, color="red", size=0.2)+
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar")+
        xlim(0,100)+
        ylim(0,100)+   
        coord_equal() + theme_classic()
    ## 
    
    matern.spde2 <- inla.spde2.pcmatern(mesh,
                                        prior.sigma = c(2, 0.01),
                                        prior.range = c(20, 0.01))



    ## fit.oscillating.rgeneric.no.samplers <- lgcp(cmp.oscillating.rgeneric,
    ## data = Y.spdf, samplers = NULL, domain = list(coordinates = mesh), options=list(verbose = TRUE))    
    ## pr.int <- predict(fit.oscillating.rgeneric.no.samplers, pxl, ~ mySPDE + Intercept)
    
    ## 
    ## cmp.matern <- coordinates ~ mySPDE(coordinates, model = matern.spde2, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
    ## fit.matern <- lgcp(cmp.matern, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh), options=list(verbose = TRUE))    

    ## 
    ## spde2
    ##
    oscillating.spde2 = inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2, 
                                           B0 = B.phi0.oscillating, B1 = B.phi1.oscillating, B2 = B.phi2.oscillating, theta.mu = c(0, log(20), -5), 
                                           theta.Q = diag(c(10, 10, 10)), transform = "log")
    cmp.oscillating.spde2 <- coordinates ~ mySPDE(coordinates, model = oscillating.spde2, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
    fit.oscillating.spde2 <- lgcp(cmp.oscillating.spde2, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh), options=list(verbose = TRUE))

    plot(fit.oscillating.spde2$marginals.hyperpar[[1]], type="l", xlim=c(-5,5))
    lines(fit.oscillating.spde2$marginals.hyperpar[[2]], lty=2)
    lines(fit.oscillating.spde2$marginals.hyperpar[[3]], lty=2)
}






## theta.hd <- seq(2.5,7.5,len=5)
## theta.temp <- seq(1,10,len=10)
## theta.var <- seq(0.1,4, len=4)
## if(FALSE){
##     theta.rho <- seq(15,25,len=10)
##     val.prof <- NULL
##     for(i in 1:length(theta.rho)){    
##         par.theta.profile <- c(log(22.9),
##                                log(0.41),
##                                -log((1-(-0.92))/(1+(-0.92))),
##                                log(3.37),
##                                log(3.19),
##                                log(theta.rho[i]),
##                                log(0.79))
##         val.prof[i] <- pthetapc.prop.marg.post_osc_temp(par.theta.profile, data=data, 
##                                                         X=Xinit, Z=Zinit, beta=betaest,
##                                                         mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
##                                                         A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                                                         W=W, 
##                                                         acc=1e-3)    
##     }
## }
## fit.space <- lgcp(cmp.space, data = Y.spdf@data, domain = list(firing_times = mesh1d),                      
##                   options=list(verbose = TRUE, bru_max_iter=1))
## space.rgeneric <- inla.rgeneric.define(oscillating.model,
##                                        M=list(M0.space=M0, M1.space=M1, M2.space=M2))
## time(firing_times, model="ou")
## samplers.space <- cbind(Ypos.sldf@data$time, Ypos.sldf@data$time.lead)
## ips  <- ipoints(samplers.space.direction, mesh1d)
##
## ----------------------------------
## load("modelfit.RData")
## Xest       <-   Xinit    <- Matrix(rep(0, ncol(A)), ncol=1)
## Zest       <-   Zinit    <- Matrix(rep(0, mesh1d$n), ncol=1)
## betaest       <-   betainit <- 0
## nrow(Y)/sum((Ypos %>% mutate(speed.lead = lead(speed), dt=c(diff(time)[1], diff(time) ))  %>%
##            head(-1) %>%
##            mutate(val=dt*((speed + speed.lead)/2)))$val)
                                        #(number of firing events)/
                                        #(\int_\Gamma ds)
## gradest       <- NULL
## Hessianest    <- NULL
## set.seed(111086)
## par.theta <- c(log(21.77),
##                log(0.48),
##                -log((1-(-0.91))/(1+(-0.91))),
##                log(3.23),
##                log(3.01),
##                log(23.36),
##                log(0.79))
## ## old implementation
## if(FALSE){
## opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_temp, data=data, 
##                    X=Xinit, Z=Zinit, beta=betainit,
##                    mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
##                    A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                    W=W, acc=1e-2,
##                    print.verbose=TRUE,
##                    control=list(maxit=5000))
## save(opt.theta, Xest, Zest, betaest, Hessianest, file="fitted_model.RData")
## }
