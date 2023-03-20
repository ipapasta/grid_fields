set.seed(111086) 
library(tidyverse)
library(purrr)
library(INLA)  
library(inlabru)
library(sp)
library(fields)
library(pals)
## inla.setOption(pardiso.license = "/Users/ipapasta/pardiso.lic")
simulation.hpp <- FALSE
source("load_data_analysis.R")
source("Functions.R")
bru_options_set(inla.mode = "experimental")
bru_options_set(control.compute = list(openmp.strategy="huge"))
##}else{
##  bru_options_set(control.compute = list(openmp.strategy="pardiso"))
## }
## ===============================================================================
mycoords.mesh2d <- cbind(runif(65*65, 0, 100), runif(65*65, 0, 100))
k         <- 4.5
mesh      <- inla.mesh.2d(mycoords.mesh2d, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
## k         <- 4
## mycoords.mesh2d <- expand.grid(seq(0, 100, len=50), seq(0,100, len=50))
## mesh      <- inla.mesh.2d(mycoords.CV, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
## mesh      <- inla.mesh.2d(mycoords.mesh2d, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
X         <- X.train
Y         <- Y.train
##
## size of discretized fields 
p           <- mesh$n
## circular mesh
mesh.hd     <- inla.mesh.1d(seq(0, 2*pi, len=21), boundary="cyclic", degree=1)
p.hd        <- mesh.hd$n

if(for_haavard){
    k         <- 10
    mesh      <- inla.mesh.2d(mycoords.mesh2d, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
    mesh.hd   <- inla.mesh.1d(seq(0, 2*pi, len=14), boundary="cyclic", degree=1)
}

## line splits
Ypos.ls         <- split.segments.wrapper.function(X=X, mesh=mesh, mesh.hd =mesh.hd)
Yposraw.ls      <- split.segments.wrapper.function(X=dat$X, mesh=mesh, mesh.hd =mesh.hd)
Ypos            <- Ypos.ls$Ypos
Yposraw         <- Yposraw.ls$Ypos
filter.index    <- Ypos.ls$filter.index
filterraw.index <- Yposraw.ls$filter.index
 
## ------------------
## Integration points
## ------------------
## 
## CHECK always that: 
##
## train.index <- sort(sample(1:K, size = floor(K/2), replace=FALSE))
## dim(coordsraw.trap)[1] == length(HDraw.data); length(HDraw.data) == length(Traw.data); length(Traw.data) == length(rawindex.CV)
coordsraw.trap  <- rbind(do.call("rbind",Yposraw$sp), tail(do.call("rbind",Yposraw$ep),1))
Liraw.trap      <- c(unlist(Yposraw$Li), tail(unlist(Ypos$Li),1))
HDraw.data      <- c(do.call("c", (Yposraw %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Yposraw$hd.lead, 1))
Traw.data       <- c(do.call("c", (Yposraw %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))
diffTraw.data   <- c(unlist(Yposraw$Ti), tail(unlist(Ypos$Ti),1))
rawindex.CV     <- c((Yposraw%>% dplyr::select(Ti, index.CV) %>% unnest(Ti))$index.CV, Yposraw$index.CV[length(Yposraw$index.CV)])

df.int.points       <- data.frame(coordsraw.trap=I(coordsraw.trap),
                                  Liraw.trap       = Liraw.trap,
                                  HDraw.data       = HDraw.data,
                                  Traw.data        = Traw.data,
                                  diffTraw.data    = diffTraw.data,
                                  rawindex.CV      = rawindex.CV)
df.int.points.train <- df.int.points %>% dplyr::filter(rawindex.CV %in% train.index)
df.int.points.test  <- df.int.points %>% dplyr::filter(!(rawindex.CV %in% train.index))

coords.trap <- df.int.points.train$coordsraw.trap           
HD.data     <- df.int.points.train$HDraw.data                     
T.data      <- df.int.points.train$Traw.data                     
index.CV    <- df.int.points.train$rawindex.CV

## 
## test set
## 
## dim(coords.trap.test)[1] == length(HD.data.test); length(HD.data.test) == length(T.data.test); length(T.data.test) == length(index.CV.test)
coords.trap.test <- df.int.points.test$coordsraw.trap           
HD.data.test     <- df.int.points.test$HDraw.data                     
T.data.test      <- df.int.points.test$Traw.data                     
index.CV.test    <- df.int.points.test$rawindex.CV


## mesh for temporal process
## print(paste("mesh1d:", head(diff(mesh1d$loc))))
mesh1d       <- inla.mesh.1d(loc=c(Traw.data[seq(1, length(Traw.data), by = 25)], Traw.data[length(Traw.data)]), order=2)


## Matrix of basis function evaluations for positional data (INTEGRAL term of Poisson likelihood)
Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
A      <- inla.row.kron(Ahd, Aosc)
## rownames(Atilde) <- rownames(Aosc) <- rownames(Ahd) <- rownames(A) <- index.CV


## ## Matrix of basis function evaluations for observed firing events (SUM term of Poisson likelihood)
Atildeobs <- inla.spde.make.A(mesh=mesh1d, Y$firing_times)
Aosc.obs  <- inla.spde.make.A(mesh=mesh, loc=as.matrix(Y %>% dplyr:: select(position_x, position_y)))
Ahd.obs   <- inla.spde.make.A(mesh=mesh.hd, Y$hd)
Aobs      <- inla.row.kron(Ahd.obs, Aosc.obs)


## Line segment lengths and Time interval lengths
## ----------------------------------------------
dGamma        <- df.int.points.train$Li
dT            <- df.int.points.train$Ti
dGamma.test   <- df.int.points.test$Li
dT.test       <- df.int.points.test$Ti
dGamma.raw    <- df.int.points$Li
dT.raw        <- df.int.points$Ti

## spatial basis functions
## training set
Aosctmp <- as(Aosc, "dgTMatrix")
Aosc.indices <- cbind(cbind(Aosctmp@i+1, Aosctmp@j+1), Aosctmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
Aosc.indices <- Aosc.indices[order(Aosc.indices[,1]),] %>% as.data.frame #
Aosc.indices <- Aosc.indices %>% mutate(index.CV = sort(rep(index.CV, 3)))

## spatial-directional basis functions
## training set
Atmp <- as(A, "dgTMatrix")
A.indices <- cbind(cbind(Atmp@i+1, Atmp@j+1), Atmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
A.indices <- A.indices[order(A.indices[,1]),] %>% as.data.frame #
A.indices <- A.indices %>% mutate(index.CV = sort(rep(index.CV, 6)))


## temporal basis functions
## training set
Attmp <- as(Atilde, "dgTMatrix")
At.indices <- cbind(cbind(Attmp@i+1, Attmp@j+1), Attmp@x) # (i,j, A[i,j]) for which Atilde[i,j] is non-zero (Time)
At.indices <- At.indices[order(At.indices[,1]),] %>% as.data.frame
At.indices <- At.indices %>% mutate(index.CV = sort(rep(index.CV, 2)))


## training set
names(Aosc.indices) <- c("tk", "i", "psi.o", "index.CV") #ot: omega 
names(A.indices)    <- c("tk", "i", "psi.ot", "index.CV") #ot: omega x theta
names(At.indices)   <- c("tk", "l", "psi.t", "index.CV")


## training set
df.prism.M0      <- df.prism.M0.wrapper_CV(Aosc.indices = Aosc.indices, dGamma=dGamma, T.data=T.data, HD.data=HD.data,
                                         coords.trap=coords.trap, index.CV=index.CV) %>% unnest(cols=c(val.M0))
df.prism.M1_M2 <- df.prism.M1.M2.wrapper2_CV(At.indices= At.indices, Aosc.indices=Aosc.indices, A.indices=A.indices,
                                             dGamma=dGamma, T.data=T.data, HD.data=HD.data, coords.trap=coords.trap, index.CV=index.CV)

df.prism.M1                      <- df.prism.M1_M2 %>% dplyr::select(-c(val.M2.space.time, val.M2.space.direction.time)) %>% unnest(cols=c(val.M1))
df.prism.M2.space.time           <- df.prism.M1_M2 %>% dplyr::select(-c(val.M1, val.M2.space.direction.time)) %>% unnest(cols=c(val.M2.space.time))
df.prism.M2.space.direction.time <- df.prism.M1_M2 %>% dplyr::select(-c(val.M1, val.M2.space.time)) %>% unnest(cols=c(val.M2.space.direction.time))



## ------------------------------------------------
## M0 model:  Integration weights integration knots 
## ------------------------------------------------
## training set
df.W.M0 <- df.prism.M0 %>% group_by(index.CV) %>% nest() %>% 
    mutate(W.M0.vector = map(data, function(x){
        tol <- 0
        df.dGamma.sum.k.kplus1.M0 <- rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
                                           dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0),
                                           x %>% dplyr::filter(tk!=min(tk)) %>%
                                           mutate(time=time.lag, direction=direction.lag, coords=coords.lag, group=tk-1, dGamma=0) %>%
                                           dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0)) %>%
            arrange(group) %>%
            mutate(dGamma.trap = dGamma + dGamma.lag) %>% 
            group_by(group, i) %>%
            summarize(val = sum(pmax(dGamma.trap*val.M0, tol))/2,
                      time = unique(time),
                      direction=unique(direction),
                      coords=unique(coords))  %>%
            ungroup %>% group_by(i) %>%
            summarize(val = sum(val))
    }))
mat.M0.tmp <- (do.call("rbind",df.W.M0$W.M0.vector)) %>% group_by(i) %>%
    summarize(val=sum(val))
W.M0 <- sparseVector(i=mat.M0.tmp$i, x=mat.M0.tmp$val, length=mesh$n)
W.ipoints.M0 <- as(W.M0, "sparseMatrix")
df.join.M0 <- data.frame(coords.x1 = mesh$loc[,1], coords.x2= mesh$loc[,2])
W.ipoints.M0 <- left_join(df.join.M0, data.frame(coords.x1 = mesh$loc[W.ipoints.M0@i+1,1],
                           coords.x2 = mesh$loc[W.ipoints.M0@i+1,2],
                           weight=W.ipoints.M0@x), by=c("coords.x1", "coords.x2")) %>%
    dplyr::mutate(weight = replace_na(weight, 0))
    



## ---------------
## space-direction
## ---------------
## train set
df.W.M1 <- df.prism.M1 %>% group_by(index.CV) %>% nest() %>% 
    mutate(W.M1.vector = map(data, function(x){
        tol <- 0
        df.dGamma.sum.k.kplus1.M1 <- rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
              dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1),
              x %>% 
              filter(tk!=min(tk)) %>%
              mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                     group=tk-1,
                     dGamma=0) %>%
              dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1)) %>%
            arrange(group) %>%
            mutate(dGamma.trap = dGamma + dGamma.lag) %>%
            group_by(group, i) %>%
            summarize(val = sum(pmax(dGamma.trap*val.M1, tol))/2,
                      time = unique(time),
                      direction=unique(direction),
                      coords=unique(coords))  %>%
            ungroup %>% group_by(i) %>%
            summarize(val = sum(val))
    }))

mat.M1.tmp <- (do.call("rbind",df.W.M1$W.M1.vector)) %>% group_by(i) %>%
    summarize(val=sum(val))
W.M1       <- sparseVector(i=mat.M1.tmp$i, x=mat.M1.tmp$val, length=mesh$n * mesh.hd$n)
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
W.ipoints.M1 <- as(W.M1, "sparseMatrix")
df.join.M1   <- data.frame(hd        = sort(rep(mesh.hd$loc, mesh$n)),
                           coords.x1 = rep(mesh$loc[,1], mesh.hd$n),
                           coords.x2 = rep(mesh$loc[,2], mesh.hd$n)) 
W.ipoints.M1 <- left_join(df.join.M1, data.frame(hd  = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,1],
                                                 coords.x1 = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,2],
                                                 coords.x2 = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,3],
                                                 weight=W.ipoints.M1@x), by=c("hd", "coords.x1", "coords.x2")) %>%
    dplyr::mutate(weight = replace_na(weight, 0))



## space-time
## training set
df.W.M2.space.time <- df.prism.M2.space.time %>% group_by(index.CV) %>% nest() %>%
    mutate(W.M2.space.time = map(data, function(x){
        tol <- 0
        df.dGamma.sum.k.kplus1.M2.space.time <- rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
                                                      dplyr::select(group, time, coords, dGamma, dGamma.lag, l, i, val.M2),
                                                      x %>% 
                                                      filter(tk!=min(tk)) %>%
                                                      mutate(time=time.lag, coords=coords.lag,
                                                             group=tk-1,
                                                             dGamma=0) %>%
                                                      dplyr::select(group, time, coords, dGamma, dGamma.lag, l, i, val.M2)) %>%
            arrange(group) %>%
            mutate(dGamma.trap = dGamma + dGamma.lag) %>%
            group_by(group, l, i) %>%
            summarize(val = sum(pmax(dGamma.trap*val.M2, tol)),
                      time = unique(time),
                      coords=unique(coords))
    }))
## tol <- 0
mat.tmp.space.time <- (do.call("rbind",df.W.M2.space.time$W.M2.space.time))
W <- sparseMatrix(i=mat.tmp.space.time$l,
                  j=mat.tmp.space.time$i,
                  x=mat.tmp.space.time$val/2,
                  dims=c(mesh1d$n, mesh$n))
df.join.M1   <- data.frame(hd        = sort(rep(mesh.hd$loc, mesh$n)),
                           coords.x1 = rep(mesh$loc[,1], mesh.hd$n),
                           coords.x2 = rep(mesh$loc[,2], mesh.hd$n)) 
## W <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(Aosc)-ncol(W)))
## W <- W %>% rbind(Matrix(0, nrow=ncol(Attmp)-nrow(W), ncol=ncol(W)))                
## 
## W.ipoints.M2.space.direction matrix has the correct format to be used in inlabru
## 
W.ipoints.M2.space.time <- as(W, "dgTMatrix")
W.ipoints.M2.space.time <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2.space.time@i+1],
                                      coords.x1 = mesh$loc[W.ipoints.M2.space.time@j+1,1],
                                      ## mapindex2space.direction_basis(W.ipoints.M2.space.time@j+1)[,2]
                                      coords.x2 =mesh$loc[W.ipoints.M2.space.time@j+1,2],
                                      ## mapindex2space.direction_basis(W.ipoints.M2.space.time@j+1)[,3]
                                      weight=W.ipoints.M2.space.time@x) %>% arrange(firing_times)






## space-direction-time
## train set
df.W.M2.space.direction.time <- df.prism.M2.space.direction.time %>% group_by(index.CV) %>% nest() %>%
    mutate(W.M2.vector = map(data, function(x){
        tol <- 0
        rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
              dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val.M2),
              x %>% 
              filter(tk!=min(tk)) %>%
              mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                     group=tk-1,
                     dGamma=0) %>%
              dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val.M2)) %>%
            arrange(group) %>%
            mutate(dGamma.trap = dGamma + dGamma.lag) %>% group_by(group, l, i) %>%
            summarize(val = sum(pmax(dGamma.trap*val.M2, tol)),
                      time = unique(time),
                      direction=unique(direction),
                      coords=unique(coords))
    }))

mat.tmp.space.direction.time <- (do.call("rbind",df.W.M2.space.direction.time$W.M2.vector))##  %>% group_by(l, i) %>%
    ## summarize(val=sum(val))
W.space.direction.time       <- sparseMatrix(i=mat.tmp.space.direction.time$l,
                                             j=mat.tmp.space.direction.time$i,
                                             x=mat.tmp.space.direction.time$val/2,
                                             dims=c(mesh1d$n, mesh$n * mesh.hd$n))
## 
W.ipoints.M2.space.direction.time <- as(W.space.direction.time, "dgTMatrix")
W.ipoints.M2.space.direction.time <- data.frame(firing_times = mesh1d$loc[W.ipoints.M2.space.direction.time@i+1],
                                                hd           = mapindex2space.direction_basis(W.ipoints.M2.space.direction.time@j+1)[,1],
                                                coords.x1    = mapindex2space.direction_basis(W.ipoints.M2.space.direction.time@j+1)[,2],
                                                coords.x2    = mapindex2space.direction_basis(W.ipoints.M2.space.direction.time@j+1)[,3],
                                                weight       = W.ipoints.M2.space.direction.time@x)


B.phi0.matern = matrix(c(0,1,0), nrow=1)
B.phi1.matern = matrix(c(0,0,1), nrow=1)
B.phi0.oscillating = matrix(c(0,1,0,0), nrow=1)
B.phi1.oscillating = matrix(c(0,0,1,0), nrow=1)
B.phi2.oscillating = matrix(c(0,0,0,1), nrow=1)
fem.mesh    <- inla.mesh.fem(mesh, order = 2)
fem.mesh.hd <- inla.mesh.fem(mesh.hd, order = 2)
fem.mesh.temporal <- inla.mesh.fem(mesh1d, order = 2)
M0 = fem.mesh$c0
M1 = fem.mesh$g1
M2 = fem.mesh$g2
M0.temporal = fem.mesh.temporal$c0
M1.temporal = fem.mesh.temporal$g1
M2.temporal = fem.mesh.temporal$g2
M0.hd = fem.mesh.hd$c0
M1.hd = fem.mesh.hd$g1
M2.hd = fem.mesh.hd$g2


Y.spdf    <- SpatialPointsDataFrame(coords = SpatialPoints(cbind(Y$position_x, Y$position_y)),
                                    data   = as.data.frame(Y%>%dplyr::select(-c(position_x, position_y))))
Ypos.sldf <- SpatialLinesDataFrame(sl   = SpatialLines(lapply(as.list(1:nrow(Ypos)),
                                                              function(k) Lines(list(Line(cbind(c(Ypos$coords[k,1],
                                                                                                  Ypos$coords.lead[k,1]),
                                                                                                c(Ypos$coords[k,2],
                                                                                                  Ypos$coords.lead[k,2])))), ID=k))),
                                   data = Ypos %>% dplyr::select(-c(coords, coords.lead)))

data <- list(Ypos=Ypos, Y=Y, Yspdf=Y.spdf, Ypos.sldf = Ypos.sldf)

## ------------------------------------------------------
## specification of prior distribution of hyperparameters
## ------------------------------------------------------
## spatial model
sigma.range.spatial.oscillating <- .4
mu.range.spatial.oscillating    <- 20
sigma.spatial.oscillating       <- 1/2
a.par.phi.prior.spatial.oscillating <- 2
b.par.phi.prior.spatial.oscillating <- 20
## directional model
rho.directional   <- 1/(2*pi)
sigma.directional <- 1
## 
rho.temporal      <- 1/100
sigma.temporal    <- 1/3
rho.prior.rate.time   <- 1/100
sigma.prior.rate.time <- 1/3
initial.space     <- list(theta1=2.6,theta2=0.5, theta3=-1.4)
initial.direction <- list(theta4=log(pi), theta5=log(1))
initial.time      <- list(theta1=log(100), theta2=log(3))
l = -0.98
u = 1
## 
Pl.Omega         <- Polygon(expand.grid(c(min(X$position_x),max(X$position_x)), c(min(X$position_y),max(X$position_y)))[c(1,2,4,3),])
ID.Omega         <- "Omega"
Pls.Omega        <- Polygons(list(Pl.Omega), ID=ID.Omega)
SPls.Omega       <- SpatialPolygons(list(Pls.Omega))
weights.domain   <- ipoints(domain=mesh, samplers=SPls.Omega)
weights.domain1d <- ipoints(domain=mesh1d, samplers=c(0,max(dat$X$synced_time)))
locs             <- weights.domain@coords
rownames(locs)   <- NULL
source("rgeneric_models.R")

## ----------------
## Fitting M0 model
## ---------------- 
space.rgeneric  <- inla.rgeneric.define(oscillating.model,
                                        M = list(M0=M0, M1=M1, M2=M2),
                                        theta.functions = list(theta.2.phi   = theta.2.phi,
                                                               theta.2.sigma = theta.2.sigma,
                                                               theta.2.rho   = theta.2.rho,
                                                               l=l, u=u),
                                        hyperpar = list(
                                            mu.range.spatial.oscillating        = mu.range.spatial.oscillating,
                                            sigma.range.spatial.oscillating     = sigma.range.spatial.oscillating,
                                            sigma.spatial.oscillating           = sigma.spatial.oscillating,
                                            a.par.phi.prior.spatial.oscillating = a.par.phi.prior.spatial.oscillating,
                                            b.par.phi.prior.spatial.oscillating = b.par.phi.prior.spatial.oscillating),
                                        prior.functions = list(prior.phi_osc = prior.phi_osc),
                                        initial.space=initial.space)

A.spatial.field_constr <- inla.spde.make.A(mesh=mesh, loc=locs,
                                           weights=weights.domain@data[,1],
                                           block=rep(1, nrow(weights.domain@coords)))
cmp.space <- firing_times ~
    spde2(cbind(coords.x1, coords.x2), model=space.rgeneric, mapper=bru_mapper(mesh,indexed
                                                                               =TRUE),
          extraconstr=list(A=as.matrix(A.spatial.field_constr,nrow=1), e=0)) + Intercept(1)

fit.space <- lgcp(cmp.space,
                  data = Y.spdf,
                  ips     = W.ipoints.M0 %>% dplyr::filter(weight!=0),
                  domain  = list(firing_times = mesh1d),
                  options = list(num.threads=8, verbose = TRUE, bru_max_iter=1) )

## ----------------
## Fitting M1 model
## ----------------
space.direction.rgeneric <- inla.rgeneric.define(space.direction.model,
                                                 M=list(M0.space=M0, M1.space=M1, M2.space=M2,
                                                        M0.direction=M0.hd, M1.direction=M1.hd, M2.direction=M2.hd),
                                                 theta.functions = list(theta.2.rho   = theta.2.rho,
                                                                        theta.2.sigma = theta.2.sigma,
                                                                        theta.2.phi   = theta.2.phi,           
                                                                        theta.2.rho.direction = theta.2.rho.direction,
                                                                        theta.2.sigma.direction = theta.2.sigma.direction,
                                                                        l=l, u=u),
                                                 hyperpar = list(
                                                     mu.range.spatial.oscillating        = mu.range.spatial.oscillating,
                                                     sigma.range.spatial.oscillating     = sigma.range.spatial.oscillating,
                                                     sigma.spatial.oscillating           = sigma.spatial.oscillating,                               
                                                     a.par.phi.prior.spatial.oscillating = a.par.phi.prior.spatial.oscillating,
                                                     b.par.phi.prior.spatial.oscillating = b.par.phi.prior.spatial.oscillating,
                                                     rho.directional                     = rho.directional,
                                                     sigma.directional                   = sigma.directional),
                                                 prior.functions = list(prior.phi_osc = prior.phi_osc),
                                                 initial.space=list(theta1 = fit.space$summary.hyperpar$mean[1],
                                                                    theta2 = fit.space$summary.hyperpar$mean[2],
                                                                    theta3 = fit.space$summary.hyperpar$mean[3]),
                                                 initial.direction = initial.direction)

A.spatial.field_constr_along.directions     <- as.matrix(kronecker(Diagonal(mesh.hd$n),
                                                                   as.matrix(inla.spde.make.A(mesh=mesh, loc=locs,
                                                                                              weights=weights.domain@data[,1],
                                                                                              block=rep(1, nrow(weights.domain@coords))), nrow=1)))

cmp.space.direction <- firing_times ~
    spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
          mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE))),
          extraconstr=list(A=as.matrix(A.spatial.field_constr_along.directions,nrow=mesh.hd$n), e=rep(0,mesh.hd$n))) +
    Intercept(1)


fit.space.direction <- lgcp(cmp.space.direction,
                            data    = Y.spdf,
                            ips     = W.ipoints.M1 %>% dplyr::filter(weight!=0),
                            domain  = list(firing_times = mesh1d),
                            options = list(num.threads=8, verbose = FALSE, bru_max_iter=1))

## ----------------
## Fitting space-time model
## ----------------
space.rgeneric     <- inla.rgeneric.define(oscillating.model,
                                           M = list(M0=M0, M1=M1, M2=M2),
                                           theta.functions = list(theta.2.phi   = theta.2.phi,
                                                                  theta.2.sigma = theta.2.sigma,
                                                                  theta.2.rho   = theta.2.rho,
                                                                  l=l, u=u),
                                           hyperpar = list(
                                               mu.range.spatial.oscillating        = mu.range.spatial.oscillating,
                                               sigma.range.spatial.oscillating     = sigma.range.spatial.oscillating,
                                               sigma.spatial.oscillating           = sigma.spatial.oscillating,
                                               a.par.phi.prior.spatial.oscillating = a.par.phi.prior.spatial.oscillating,
                                               b.par.phi.prior.spatial.oscillating = b.par.phi.prior.spatial.oscillating),
                                           prior.functions = list(prior.phi_osc = prior.phi_osc),
                                           initial.space=list(theta1 = fit.space$summary.hyperpar$mean[1],
                                                              theta2 = fit.space$summary.hyperpar$mean[2],
                                                              theta3 = fit.space$summary.hyperpar$mean[3]))

time.rgeneric            <- inla.rgeneric.define(temporal.model,
                                                 M=list(M0.temporal=M0.temporal, M1.temporal=M1.temporal, M2.temporal=M2.temporal),
                                                 theta.functions = list(theta.2.rho.time        = theta.2.rho.time,
                                                                        theta.2.sigma.time      = theta.2.sigma.time),
                                                 hyperpar = list(
                                                     sigma.prior.rate.time   = sigma.prior.rate.time,
                                                     rho.prior.rate.time     = rho.prior.rate.time),
                                                 initial.time  = initial.time)

A.temporal.field_constr <- weights.domain1d$weight

cmp.space.time <- firing_times ~
    spde2(cbind(coords.x1, coords.x2), model=space.rgeneric, mapper=bru_mapper(mesh,indexed=TRUE),
          extraconstr=list(A=as.matrix(A.spatial.field_constr,nrow=1), e=0)) + 
    spde1(firing_times, model=time.rgeneric, mapper=bru_mapper(mesh1d, indexed=TRUE),
          extraconstr=list(A=matrix(A.temporal.field_constr,nrow=1), e=0)) + Intercept(1)


fit.space.time <- lgcp(cmp.space.time, data = Y.spdf,
                       ips=W.ipoints.M2.space.time %>% dplyr::filter(weight!=0),
                       domain = list(firing_times = mesh1d),
                       options=list( num.threads=8, verbose = FALSE, bru_max_iter=1)) 




## ----------------------------------
## Fitting space-direction-time model
## ----------------------------------

space.direction.rgeneric <- inla.rgeneric.define(space.direction.model,
                                                 M=list(M0.space=M0, M1.space=M1, M2.space=M2,
                                                        M0.direction=M0.hd, M1.direction=M1.hd, M2.direction=M2.hd),
                                                 theta.functions = list(theta.2.rho   = theta.2.rho,
                                                                        theta.2.sigma = theta.2.sigma,
                                                                        theta.2.phi   = theta.2.phi,           
                                                                        theta.2.rho.direction = theta.2.rho.direction,
                                                                        theta.2.sigma.direction = theta.2.sigma.direction,
                                                                        l=l, u=u),
                                                 hyperpar = list(
                                                     mu.range.spatial.oscillating        = mu.range.spatial.oscillating,
                                                     sigma.range.spatial.oscillating     = sigma.range.spatial.oscillating,
                                                     sigma.spatial.oscillating           = sigma.spatial.oscillating,                               
                                                     a.par.phi.prior.spatial.oscillating = a.par.phi.prior.spatial.oscillating,
                                                     b.par.phi.prior.spatial.oscillating = b.par.phi.prior.spatial.oscillating,
                                                     rho.directional                     = rho.directional,
                                                     sigma.directional                   = sigma.directional),
                                                 prior.functions = list(prior.phi_osc = prior.phi_osc),
                                                 initial.space=list(theta1 = fit.space.direction$summary.hyperpar$mean[1],
                                                                    theta2 = fit.space.direction$summary.hyperpar$mean[2],
                                                                    theta3 = fit.space.direction$summary.hyperpar$mean[3]),
                                                 initial.direction = list(theta4=fit.space.direction$summary.hyperpar$mean[4],
                                                                          theta5=fit.space.direction$summary.hyperpar$mean[5]))

time.rgeneric            <- inla.rgeneric.define(temporal.model,
                                                 M=list(M0.temporal=M0.temporal, M1.temporal=M1.temporal, M2.temporal=M2.temporal),
                                                 theta.functions = list(theta.2.rho.time        = theta.2.rho.time,
                                                                        theta.2.sigma.time      = theta.2.sigma.time),
                                                 hyperpar = list(
                                                     sigma.prior.rate.time   = sigma.prior.rate.time,
                                                     rho.prior.rate.time     = rho.prior.rate.time),
                                                 initial.time  = list(theta1 = fit.space.time$summary.hyperpar$mean[1],
                                                                      theta2 = fit.space.time$summary.hyperpar$mean[2]))

cmp.space.direction.time <- firing_times ~
    spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
          mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE))),
          extraconstr=list(A=as.matrix(A.spatial.field_constr_along.directions,nrow=mesh.hd$n), e=rep(0,mesh.hd$n))) +
    spde1(firing_times, model=time.rgeneric, mapper=bru_mapper(mesh1d, indexed=TRUE),
          extraconstr=list(A=matrix(A.temporal.field_constr,nrow=1), e=0)) + Intercept(1)


fit.space.direction.time <- lgcp(cmp.space.direction.time, data = Y.spdf,
                                 ips=W.ipoints.M2.space.direction.time %>% dplyr::filter(weight!=0),
                                 domain = list(firing_times = mesh1d),
                                 options=list(
                                     num.threads=8,
                                     verbose = FALSE,
                                     bru_max_iter=1))


save(list=ls(), file="/Users/ipapasta/Documents/Research/Software/R/grid_fields/R/analysis_saved_data/fits.RData")

## ----------------------------------------------------------------
## Predictions for the expected number of firing events in test set
## ----------------------------------------------------------------
## calculation of integration weight for connected line segments
## ("clumps"); note that adjacent clumps in test data are joined
## together so that the number of test clumps is less than the
## original training/test split

## unique(Y.test$index.CV) e
## =  1  2  3  4  6  8 11 13 14 15 16 19 20 21 22 23 24 29 33 37 38 43 44 46 47 49 55 56 57 61
## length(list(1:4, c(6), c(8), c(11), 13:16, 
##       19:24, c(29), c(33), 37:38, 43:44, 46:47, c(49), 55:57, c(61))) = 14
## i.e. there are 14 test clumps after grouping adjacent selected clumps
## so for example the first clump in Y.test can be obtained by filtering for index.CV %in% 1:4
## weights_test  <- weights_line_segments_in_train(X.test=X.test, Y.test = Y.test, mesh=mesh, mesh.hd=mesh.hd, mesh1d=mesh1d)
## weights_train <- weights_line_segments_in_train(X.test=X.train, Y.test = Y.train, mesh=mesh, mesh.hd=mesh.hd, mesh1d=mesh1d)
## weights_test[[1]]$index.CV 1 6 8 11 13 19 29 33 37 43 46 49 55 60
## which is of length 14 I think 1:4 from unique(Y.test$index.CV) are
## now just labelled as 1 here etc..

## use split function to take unique(Y.test$index.CV) and returns a list of consecutive integers
clumps             <- setdiff(1:K, train.index)
    ## split(unique(Y.test$index.CV), cumsum(c(1, diff(unique(Y.test$index.CV)) != 1)))
obs.firings        <- sapply(1:length(clumps), function(i) nrow(Y.test %>% filter(index.CV %in% clumps[[i]])))
## clumps.train       <- split(unique(Y.train$index.CV), cumsum(c(1, diff(unique(Y.train$index.CV)) != 1)))
## can get number of firing events in clump 1 (after grouping) of test
## set using nrow(Y.test %>% filter(index.CV %in% clumps[[1]])) = 513
## ===============================================================================
## ===============================================================================
## Function that calculates the mean and variance for the predictive
## distribution of the number of firing events on segments/clumps of
## test set.  we calculate the mean via Monte Carlo; we can calculate
## the conditional predictive mean (conditioning on latent params;
## Intercept + GRF). Then get unconditional mean via iterated
## expectation (tower law).  For the variance of predictive
## distribution use Var(Y) = E(Var(Y | X)) + Var(E(Y | X )) where in
## our case Var(Y | X) = E(Y | X) as conditional distribution is
## Poisson.
## Inputs; weights.mat: weights matrix (will be extracted from
## weights_test object) post.sample: list of posterior samples
## (intercept & GRF); output of inlabru function "generate"
## observed.total.no.events           <- nrow(as.data.frame(Y.spdf))
## seq.N                              <- seq(round(observed.total.no.events-1000), round(observed.total.no.events+1000), by=1)
## 
##===============================================================================
##===============================================================================
## repeat test on all segments/clumps 
## posterior samples;
samp.M0.space                <- generate(object = fit.space, n.samples = 10000, num.threads=ncores)
samp.M2.space.time           <- generate(object = fit.space.time,  n.samples = 10000, num.threads=ncores)
samp.M1.space.direction      <- generate(object = fit.space.direction,  n.samples = 10000, num.threads=ncores)
samp.M2.space.direction.time <- generate(object = fit.space.direction.time,  n.samples = 10000, num.threads=ncores)
## samp.M1.space.direction2<- generate(object = fit.space.direction2,  n.samples = 5000, num.threads=8)
##
## space
clumps.mean.var.M0            <- lapply(1:length(clumps),
                                        function(i) {                                            
                                            pred.mean.var(weights.mat = df.W.M0.test$W.M0.vector[[i]]$weight,
                                                          post.sample = samp.M0.space)
                                        })
## space-direction
clumps.mean.var.M1.space.direction <- lapply(1:length(clumps),
                                             function(i) {
                                                 pred.mean.var(weights.mat = df.W.M1.test$W.M1.vector[[i]]$weight,
                                                               post.sample = samp.M1.space.direction)
                                             })
## space-time
clumps.mean.var.M2.space.time <- lapply(1:length(clumps),
                                        function(i){
                                            pred.mean.var.M2(weights.mat = df.W.M2.space.time.test$W.M2.matrix[[i]],
                                                             post.sample = samp.M2.space.time)
                                        })
## space-direction-time
clumps.mean.var.M2.space.direction.time <- lapply(1:length(clumps),
                                                  function(i){
                                                      pred.mean.var.M2(weights.mat = df.W.M2.space.direction.time.test$W.M2.matrix[[i]],
                                                                       post.sample = samp.M2.space.direction.time)
                                                  })


# predictive means on each segment/clump
pred.means.M0                       <- sapply(clumps.mean.var.M0, function(i) i[[1]]) 
pred.means.M1.space.direction       <- sapply(clumps.mean.var.M1.space.direction, function(i) i[[1]])
pred.means.M2.space.time            <- sapply(clumps.mean.var.M2.space.time, function(i) i[[1]]) 
pred.means.M2.space.direction.time  <- sapply(clumps.mean.var.M2.space.direction.time, function(i) i[[1]])

##
# predictive variance on each segment/clump
pred.vars.M0                        <- sapply(clumps.mean.var.M0, function(i) i[[2]])
pred.vars.M1.space.direction        <- sapply(clumps.mean.var.M1.space.direction, function(i) i[[2]])
pred.vars.M2.space.time             <- sapply(clumps.mean.var.M2.space.time, function(i) i[[2]])
pred.vars.M2.space.direction.time   <- sapply(clumps.mean.var.M2.space.direction.time, function(i) i[[2]])


out.list <- list(obs.firings   = obs.firings,
                 pred.means.M0 = pred.means.M0,
                 pred.vars.M0  = pred.vars.M0,
                 pred.means.M1.space.direction      = pred.means.M1.space.direction,
                 pred.vars.M1.space.direction       = pred.vars.M1.space.direction,
                 pred.means.M2.space.time           = pred.means.M2.space.time,
                 pred.vars.M2.space.time            = pred.vars.M2.space.time,
                 pred.means.M2.space.direction.time = pred.means.M2.space.direction.time,
                 pred.vars.M2.space.direction.time  = pred.vars.M2.space.direction.time,
                 mesh    = mesh,
                 mesh.hd = mesh.hd,
                 mesh1d  = mesh1d,
                 samp.M0.space                 =  samp.M0.space,               
                 samp.M2.space.time            =  samp.M2.space.time,          
                 samp.M1.space.direction       =  samp.M1.space.direction,     
                 samp.M2.space.direction.time  =  samp.M2.space.direction.time)





save(out.list,  file=paste0("/exports/eddie/scratch/ipapasta/CV_output_M0_M1_M2_", char.to.save, "_seconds_split.RData"))
save(fit.space, file=paste0("/exports/eddie/scratch/ipapasta/CV_fit_space_", char.to.save, "_seconds_split.RData"))
save(fit.space.time, file=paste0("/exports/eddie/scratch/ipapasta/CV_fit_space_time_", char.to.save, "_seconds_split.RData"))
save(fit.space.direction, file=paste0("/exports/eddie/scratch/ipapasta/CV_fit_space_direction_", char.to.save,"_seconds_split.RData"))
save(fit.space.direction.time, file=paste0("/exports/eddie/scratch/ipapasta/CV_fit_space_direction_time_", char.to.save,"_seconds_split.RData"))
## 
if(FALSE){
## space - direction, oscillating directional field
space.direction2.rgeneric <- inla.rgeneric.define(space.direction2.model,
                                                  M=list(M0.space=M0, M1.space=M1, M2.space=M2,
                                                         M0.direction=M0.hd, M1.direction=M1.hd, M2.direction=M2.hd),
                                                  theta.functions = list(theta.2.rho             = theta.2.rho,
                                                                         theta.2.sigma           = theta.2.sigma,
                                                                         theta.2.phi             = theta.2.phi,           
                                                                         theta.2.rho.direction   = theta.2.rho.direction,
                                                                         theta.2.sigma.direction = theta.2.sigma.direction,
                                                                         l=l, u=u),
                                                  hyperpar = list(
                                                      mu.range.spatial.oscillating        = mu.range.spatial.oscillating,
                                                      sigma.range.spatial.oscillating     = sigma.range.spatial.oscillating,
                                                      sigma.spatial.oscillating           = sigma.spatial.oscillating,
                                                      a.par.phi.prior.spatial.oscillating = a.par.phi.prior.spatial.oscillating,
                                                      b.par.phi.prior.spatial.oscillating = b.par.phi.prior.spatial.oscillating,
                                                      rho.directional                     = rho.directional,
                                                      sigma.directional                   = sigma.directional),
                                                  prior.functions = list(prior.phi_osc = prior.phi_osc),
                                                  initial.space     = list(theta1 = fit.space.direction$summary.hyperpar$mean[1],
                                                                           theta2 = fit.space.direction$summary.hyperpar$mean[2],
                                                                           theta3 = fit.space.direction$summary.hyperpar$mean[3]),
                                                  initial.direction = list(theta4 = fit.space.direction$summary.hyperpar$mean[4],
                                                                           theta5 = fit.space.direction$summary.hyperpar$mean[5],
                                                                           theta6 = -2))


    cmp.space.direction2 <- firing_times ~
        spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction2.rgeneric,
              mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE))),
              extraconstr=list(A=as.matrix(A.spatial.field_constr_along.directions,nrow=mesh.hd$n), e=rep(0,mesh.hd$n))) +
        Intercept

    fit.space.direction2 <- lgcp(cmp.space.direction2,
                                 data    = Y.spdf,
                                 ips     = W.ipoints.M1,
                                 domain  = list(firing_times = mesh1d),
                                 options = list( num.threads=ncores, verbose = FALSE, bru_max_iter=1) ) %>%
        bru_rerun()
    bru_options_set(inla.mode = "classic")
    fit.space.direction2 <- fit.space.direction2 %>% bru_rerun()
}


if(FALSE){
    predictive.M0.space.total.no.events      <- predictive.M0(seq=seq.N, weights.mat = W.M0.vector, post.sample = samp.space)
    predictive.M2.space.time.total.no.events <- predictive.M2(seq=seq.N, weights.mat = W.M2.space.time, post.sample = samp.space.time)
    predictive.M1.space.total.no.events      <- predictive.M1(seq=seq.N, weights.mat = W.M1.vector, post.sample = samp.space)
    ## 
    ## pred.vars.M1.space.direction2  = pred.vars.M1.space.direction2,
    ## pred.means.M1.space.direction2 = pred.means.M1.space.direction2,
    ##
    ## pred.means.M1.space.direction2 <- sapply(clumps.mean.var.M1.space.direction2, function(i) i[[1]]) 
    ## pred.vars.M1.space.direction2       <- sapply(clumps.mean.var.M1.space.direction2, function(i) i[[2]])
}

