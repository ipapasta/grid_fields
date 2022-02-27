set.seed(111086) 
library(tidyverse)
library(purrr)
library(INLA)  
library(inlabru)
library(sp)
library(fields)
library(nloptr)
library(pals)
inla.setOption(pardiso.license = "/Users/ipapasta/pardiso.lic")
source("load_data_CV.R")
source("Functions.R")
# source("osc_precision.R")
# source("hd_precision.R")
# source("temp_precision.R")
bru_options_set(inla.mode = "experimental")
## option "experimental" implements Variational Bayes correction this
## is useful for Poisson point process likelihoods.
## This option is implemented in latest INLA testing version
if(inla.pardiso.check() != "SUCCESS: PARDISO IS INSTALLED AND WORKING"){
  ## enable openmp for parallel computing using strategy "huge".
  ## this is going to engage all RAM and core resources of the
  ## computer needs some care when setting this up on Eddie.
  bru_options_set(control.compute = list(openmp.strategy="huge"))
}else{
  bru_options_set(control.compute = list(openmp.strategy="pardiso"))
}
k         <- 5
mesh      <- inla.mesh.2d(mycoords.CV, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
X <- X.train
Y <- Y.train
p           <- mesh$n
## circular mesh
mesh.hd <- inla.mesh.1d(seq(0, 2*pi, len=30), boundary="cyclic", degree=1)
p.hd    <- mesh.hd$n
## line splits
Ypos.ls      <- split.segments.wrapper.function(X=X, mesh=mesh, mesh.hd =mesh.hd)
Ypos  <- Ypos.ls$Ypos 
filter.index <- Ypos.ls$filter.index
df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)),
                         space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
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
Y.spdf    <- SpatialPointsDataFrame(coords = SpatialPoints(cbind(Y$position_x, Y$position_y)),
                                    data   = as.data.frame(Y%>%dplyr::select(-c(position_x, position_y))))
Ypos.sldf <- SpatialLinesDataFrame(sl   = SpatialLines(lapply(as.list(1:nrow(Ypos)),
                                                              function(k) Lines(list(Line(cbind(c(Ypos$coords[k,1],
                                                                                                  Ypos$coords.lead[k,1]),
                                                                                                c(Ypos$coords[k,2],
                                                                                                  Ypos$coords.lead[k,2])))), ID=k))),
                                   data = Ypos %>% dplyr::select(-c(coords, coords.lead)))
data <- list(Ypos=Ypos, Y=Y, Yspdf=Y.spdf, Ypos.sldf = Ypos.sldf)
coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))
mesh1d       <- inla.mesh.1d(loc=c(T.data[seq(1, length(T.data), by = 300)], T.data[length(T.data)]), order=2)
data         <- list(Ypos=Ypos, Y=Y, Yspdf=Y.spdf, Ypos.sldf = Ypos.sldf)
Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
A      <- inla.row.kron(Ahd, Aosc)
Atildeobs <- inla.spde.make.A(mesh=mesh1d, data$Y$firing_times)
Aosc.obs  <- inla.spde.make.A(mesh=mesh, loc=as.matrix(data$Y %>% dplyr:: select(position_x, position_y)))
Ahd.obs   <- inla.spde.make.A(mesh=mesh.hd, data$Y$hd)
Aobs      <- inla.row.kron(Ahd.obs, Aosc.obs)


dGamma <- c(do.call("c", Ypos$Li))
dT  <- diff(T.data)

## spatial basis functions
Aosctmp <- as(Aosc, "dgTMatrix")
Aosc.indices <- cbind(cbind(Aosctmp@i+1, Aosctmp@j+1), Aosctmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
Aosc.indices <- Aosc.indices[order(Aosc.indices[,1]),] %>% as.data.frame # 

## spatial-directional basis functions
Atmp <- as(A, "dgTMatrix")
A.indices <- cbind(cbind(Atmp@i+1, Atmp@j+1), Atmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
A.indices <- A.indices[order(A.indices[,1]),] %>% as.data.frame #

## temporal basis functions
Attmp <- as(Atilde, "dgTMatrix")
At.indices <- cbind(cbind(Attmp@i+1, Attmp@j+1), Attmp@x) # (i,j, A[i,j]) for which Atilde[i,j] is non-zero (Time)
At.indices <- At.indices[order(At.indices[,1]),] %>% as.data.frame

names(Aosc.indices) <- c("tk", "i", "psi.o") #ot: omega 
names(A.indices) <- c("tk", "i", "psi.ot") #ot: omega x theta
names(At.indices) <- c("tk", "l", "psi.t") 


df.prism.M0    <- df.prism.M0.wrapper(Aosc.indices = Aosc.indices, dGamma=dGamma, T.data=T.data, HD.data=HD.data,
                                      coords.trap=coords.trap) %>% unnest(cols=c(val.M0))
df.prism.M1_M2 <- df.prism.M1.M2.wrapper(At.indices= At.indices, A.indices=A.indices, dGamma=dGamma, T.data=T.data, HD.data=HD.data, coords.trap=coords.trap)
df.prism.M1    <- df.prism.M1_M2 %>% dplyr::select(-val.M2) %>% unnest(cols=c(val.M1))
df.prism.M2    <- df.prism.M1_M2 %>% dplyr::select(-val.M1) %>% unnest(cols=c(val.M2))


## ------------------------------------------------
## M0 model:  Integration weights integration knots 
## ------------------------------------------------
df.W.M0 <- rbind(df.prism.M0 %>% mutate(group=tk, dGamma.lag=0) %>%
                   dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0),
                 df.prism.M0 %>% 
                   filter(tk!=1) %>%
                   mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                          group=tk-1,
                          dGamma=0) %>%
                   dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0)) %>%
  arrange(group) %>%
  mutate(dGamma.trap = dGamma + dGamma.lag) 

tol <- 0


df.dGamma.sum.k.kplus1.M0 <- df.W.M0 %>% group_by(group, i) %>%
  summarize(val = sum(max(dGamma.trap*val.M0, tol))/2,
            time = unique(time),
            direction=unique(direction),
            coords=unique(coords))  %>%
  ungroup %>% group_by(i) %>%
  summarize(val = sum(val))

W.M0 <- sparseVector(i=df.dGamma.sum.k.kplus1.M0$i,
                     x=df.dGamma.sum.k.kplus1.M0$val,
                     length=mesh$n)

W.M0.vector <- sparseVector(i=df.dGamma.sum.k.kplus1.M0$i,
                            x=df.dGamma.sum.k.kplus1.M0$val,
                            length=mesh$n)


W.ipoints.M0 <- as(W.M0, "sparseMatrix")
W.ipoints.M0 <- data.frame(coords.x1 = mesh$loc[W.ipoints.M0@i+1,1],
                           coords.x2 = mesh$loc[W.ipoints.M0@i+1,2],
                           weight=W.ipoints.M0@x) 

df.W.M1 <- rbind(df.prism.M1 %>% mutate(group=tk, dGamma.lag=0) %>%
                   dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1),
                 df.prism.M1 %>% 
                   filter(tk!=1) %>%
                   mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                          group=tk-1,
                          dGamma=0) %>%
                   dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1)) %>%
  arrange(group) %>%
  mutate(dGamma.trap = dGamma + dGamma.lag) 

tol <- 0
df.dGamma.sum.k.kplus1.M1 <- df.W.M1 %>% group_by(group, i) %>%
  summarize(val = sum(max(dGamma.trap*val.M1, tol))/2,
            time = unique(time),
            direction=unique(direction),
            coords=unique(coords))  %>%
  ungroup %>% group_by(i) %>%
  summarize(val = sum(val))

W.M1 <- sparseVector(i=df.dGamma.sum.k.kplus1.M1$i,
                     x=df.dGamma.sum.k.kplus1.M1$val,
                     length=mesh$n * mesh.hd$n)

W.M1.vector <- sparseVector(i=df.dGamma.sum.k.kplus1.M1$i,
                            x=df.dGamma.sum.k.kplus1.M1$val,
                            length=mesh$n * mesh.hd$n)

df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)), space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
## TODO: add description for following functions
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
W.ipoints.M1 <- data.frame(hd=mapindex2space.direction_basis(W.ipoints.M1@i+1)[,1],
                           coords.x1 =mapindex2space.direction_basis(W.ipoints.M1@i+1)[,2],
                           coords.x2 =mapindex2space.direction_basis(W.ipoints.M1@i+1)[,3],
                           weight=W.ipoints.M1@x) 


df.W.M2 <- rbind(df.prism.M2 %>% mutate(group=tk, dGamma.lag=0) %>%
                   dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val.M2),
                 df.prism.M2 %>% 
                   filter(tk!=1) %>%
                   mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                          group=tk-1,
                          dGamma=0) %>%
                   dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val.M2)) %>%
  arrange(group) %>%
  mutate(dGamma.trap = dGamma + dGamma.lag) ## %>%

tol <- 0

df.dGamma.sum.k.kplus1.M2 <- df.W.M2 %>% group_by(group, l, i) %>%
  summarize(val = sum(max(dGamma.trap*val.M2, tol)),
            time = unique(time),
            direction=unique(direction),
            coords=unique(coords))

W <- sparseMatrix(i=df.dGamma.sum.k.kplus1.M2$l,
                  j=df.dGamma.sum.k.kplus1.M2$i,
                  x=df.dGamma.sum.k.kplus1.M2$val/2)

W         <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(A)-ncol(W)))
W.ipoints.M2 <- as(W, "dgTMatrix")
W.ipoints.M2 <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2@i+1], hd=mapindex2space.direction_basis(W.ipoints.M2@j+1)[,1],
                           coords.x1 =mapindex2space.direction_basis(W.ipoints.M2@j+1)[,2],
                           coords.x2 =mapindex2space.direction_basis(W.ipoints.M2@j+1)[,3],
                           weight=W.ipoints.M2@x) %>% arrange(firing_times)

## B.phi0.matern = matrix(c(0,1,0), nrow=1)
## B.phi1.matern = matrix(c(0,0,1), nrow=1)
## B.phi0.oscillating = matrix(c(0,1,0,0), nrow=1)
## B.phi1.oscillating = matrix(c(0,0,1,0), nrow=1)
## B.phi2.oscillating = matrix(c(0,0,0,1), nrow=1)
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




## ------------------------------------------------------
## specification of prior distribution of hyperparameters
## ------------------------------------------------------
## spatial model
sigma.range.spatial.oscillating <- .5
mu.range.spatial.oscillating    <- 20
sigma.spatial.oscillating       <- 1/5
a.par.phi.prior.spatial.oscillating <- 2
b.par.phi.prior.spatial.oscillating <- 20
## directional model
rho.directional   <- 1
sigma.directional <- .1
## 
rho.temporal  <- 1/100
sigma.temporal <- 1/3
initial.space <- list(theta1=log(21-5), theta2=-2, theta3=-3)
initial.space.direction <- list(theta1=log(21-5),theta2=-2, theta3=-3, theta4=log(3), theta5=-3)
l = -0.98
u = 1
weights.domain <- ipoints(domain=mesh)
## plot(seq(-.99, .99, len=100), prior.phi_osc(seq(-.99, .99, len=100), 5, 30) %>% exp %>% log)
## source all custom-made built models for inla.rgeneric.define
source("rgeneric_models.R")

## ----------------
## Fitting M0 model
## ----------------
space.rgeneric     <- inla.rgeneric.define(oscillating.model,
                                           M = list(M0=M0.space, M1=M1.space, M2=M2.space),
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

cmp.space <- firing_times ~
    f(cbind(coords.x1, coords.x2), model=space.rgeneric, mapper=bru_mapper(mesh,indexed=TRUE),
      extraconstr=list(A=t(weights.domain$weight), e=0)) + Intercept

fit.space <- lgcp(cmp.space,
                  data = Y.spdf,
                  ips     = W.ipoints.M0,
                  domain  = list(firing_times = mesh1d),
                  options = list( num.threads=8,verbose = TRUE, bru_max_iter=1) )

## ----------------
## Fitting M1 model
## ----------------
space.direction.rgeneric <- inla.rgeneric.define(space.direction.model,
                                                 M=list(M0.space=M0.space, M1.space=M1.space, M2.space=M2.space,
                                                        M0.direction=M0.hd, M1.direction=M1.hd, M2.direction=M2.hd),
                                                 theta.functions = list(theta.2.phi           = theta.2.phi,
                                                                        theta.2.sigma         = theta.2.sigma,
                                                                        theta.2.rho           = theta.2.rho,
                                                                        theta.2.rho.direction = theta.2.rho.direction,
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
                                                 initial.space.direction=initial.space.direction)

A.spatial.field_constr     <- as.matrix(kronecker(Diagonal(mesh.hd$n), t(weights.domain$weight)))
A.directional.field_constr <- circulant(rep(rep(c(diff(mesh.hd$loc)[1], rep(0,(mesh$n)-1))), mesh.hd$n))[1:mesh$n, ]
A.directional.field.and.spatial.field_constr <- rbind(A.spatial.field_constr, A.directional.field_constr)


cmp.space.direction <- firing_times ~
    f(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
      mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE))),
      extraconstr=list(A=A.spatial.field_constr, e=rep(0, nrow(A.spatial.field_constr)))) + Intercept

## cmp.space.direction <- firing_times ~
##     f(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
##       mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE))),
##       extraconstr=list(A=as.matrix(A.directional.field.and.spatial.field_constr),
##                        e=rep(0,nrow(A.directional.field.and.spatial.field_constr)))) + Intercept

fit.space.direction <- lgcp(cmp.space.direction, data = Y.spdf,
                            ips     = W.ipoints.M1,
                            domain  = list(firing_times = mesh1d),
                            options = list( num.threads=8, verbose = TRUE, bru_max_iter=1) )


## ----------------
## Fitting M2 model (computationally expensive)
## ----------------
##
## TODO: update time.rgeneric to include theta mappings, priors,
## etc. as arguments.
time.rgeneric            <- inla.rgeneric.define(temporal.model,
                                                 M=list(M0.temporal=M0.temporal, M1.temporal=M1.temporal, M2.temporal=M2.temporal),
                                                 hyperpar = list(
                                                     rho.temporal   = rho.temporal,
                                                     sigma.temporal = sigma.temporal
                                                 ))

if(FALSE){
    ## current implementation below is correct
    cmp.space.direction.time <- firing_times ~
        f(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
              mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE)))) +
        time(firing_times, model=time.rgeneric, mapper=bru_mapper(mesh1d, indexed=TRUE)) + Intercept

    fit.space.direction.time <- lgcp(cmp.space.direction.time, data = as.data.frame(Y.spdf),
                                     ips=W.ipoints.M2,
                                     domain = list(firing_times = mesh1d),
                                     options=list(
                                         num.threads=8,
                                         verbose = TRUE, bru_max_iter=1))
}

## ----------------------------------------------------------------
## Predictions for the expected number of firing events in test set
## ----------------------------------------------------------------
## calculation of integration weight for connected line segments ("clumps"); 
## note that adjacent clumps in test data are joined together so that the number 
## of test clumps is less than the original training/test split
weights_test <- weights_line_segments_in_train(X.test=X.test, Y.test = Y.test, mesh=mesh, mesh.hd=mesh.hd, mesh1d=mesh1d)

# unique(Y.test$index.CV) 
# =  1  2  3  4  6  8 11 13 14 15 16 19 20 21 22 23 24 29 33 37 38 43 44 46 47 49 55 56 57 61
# length(list(1:4, c(6), c(8), c(11), 13:16, 
#       19:24, c(29), c(33), 37:38, 43:44, 46:47, c(49), 55:57, c(61))) = 14
# i.e. there are 14 test clumps after grouping adjacent selected clumps
# so for example the first clump in Y.test can be obtained by filtering for index.CV %in% 1:4

# weights_test[[1]]$index.CV 1  6  8 11 13 19 29 33 37 43 46 49 55 60 which is of length 14
# I think 1:4 from unique(Y.test$index.CV) are now just labelled as 1 here etc..

# use split function to take unique(Y.test$index.CV) and returns a list of consecutive integers
clumps <- split(unique(Y.test$index.CV), cumsum(c(1, diff(unique(Y.test$index.CV)) != 1)))
# can get number of firing events in clump 1 (after grouping) of test set 
# using nrow(Y.test %>% filter(index.CV %in% clumps[[1]])) = 513
#===============================================================================
#===============================================================================
# Function that calculates the mean and variance for the predictive distribution of the number 
# of firing events on segments/clump of test set.
# we calculate the mean via Monte Carlo; we can calculate the conditional predictive mean
# (conditioning on latent params; Intercept + GRF). Then get unconditional mean via iterated expectation (tower law).
# For the variance of predictive distribution use Var(Y) = E(Var(Y | X)) + Var(E(Y | X ))
# where in our case Var(Y | X) =  E(Y | X) as conditional distribution is Poisson. 

# Inputs; 
# mod: fitted inlabru model object
# newdata: data of covariates used for sampling (not in use)
# weights.mat: weights matrix (will be extracted from weights_test object)
# n.samples: number of samples from posterior of latent params to be used

pred.mean.var <- function(mod, newdata=NULL, weights.mat, n.samples=10^3){
  # posterior sample
  post.samp <- generate(object = mod, 
                        data = newdata,
                        n.samples = n.samples)
  Int <- lapply(1:n.samples, function(i) post.samp[[i]]$Intercept) # simulated intercepts
  GRF <- lapply(1:n.samples, function(i) post.samp[[i]]$f)     # simulated GRFs 
  cond.mean <- sapply(1:n.samples, function(i) exp(Int[[i]])*sum(exp(GRF[[i]])*weights.mat)) # posterior means conditioned on latent parameters
  post.mean <- mean(cond.mean)
  post.var  <- mean(cond.mean) + var(cond.mean)   
  list(post.mean=post.mean, post.var=post.var)
  }
#===============================================================================
#===============================================================================
# test it on segment/clump 1
test1 <- pred.mean.var(mod=fit.space, 
               weights.mat=weights_test$df.W.M0$W.ipoints.M0[[1]]$W.M0,
               n.samples=100)
test1 
#===============================================================================
#===============================================================================
# repeat test on all segments/clumps 
# when running on Eddie we will use a much larger value for n.samples to make Monte Carlo error negligible
# can also speed things up by parallelizing some lapply calls 
clumps.mean.var <- lapply(1:length(clumps), function(i) pred.mean.var(mod=fit.space, 
                                                                      weights.mat=weights_test$df.W.M0$W.ipoints.M0[[i]]$W.M0,
                                                                      n.samples=100))
pred.means <- sapply(clumps.mean.var, function(i) i[[1]]) # predictive means on each segment/clump
pred.vars  <- sapply(clumps.mean.var, function(i) i[[2]]) # .......... variance ..........
# compare means with observed firings
obs.firings <- sapply(1:length(clumps), function(i) nrow(Y.test %>% filter(index.CV %in% clumps[[i]])))
# compare observed and expected firings
obs.firings
floor(pred.means)

# squared error score per segment/clump
se.score <- (pred.means - obs.firings)^2
# Dawid-Sebastiani score
ds.score <- se.score / pred.vars + log(pred.vars)
