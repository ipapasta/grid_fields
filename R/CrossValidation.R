## most recent implementation by Graeme

set.seed(111086) 
library(tidyverse)
library(purrr)
library(INLA)  
library(inlabru)
library(sp)
library(fields)
library(pals)
                                        #inla.setOption(pardiso.license = "/Users/ipapasta/pardiso.lic")
source("load_data_CV.R")
source("Functions.R")
                                        # source("osc_precision.R")
                                        # source("hd_precision.R")
                                        # source("temp_precision.R")
bru_options_set(inla.mode = "experimental")
## option "experimental" implements Variational Bayes correction this
## is useful for Poisson point process likelihoods.
## This option is implemented in latest INLA testing version
##if(inla.pardiso.check() != "SUCCESS: PARDISO IS INSTALLED AND WORKING"){
## enable openmp for parallel computing using strategy "huge".
## this is going to engage all RAM and core resources of the
## computer needs some care when setting this up on Eddie.
##bru_options_set(control.compute = list(openmp.strategy="huge"))
##}else{
##  bru_options_set(control.compute = list(openmp.strategy="pardiso"))
                                        #}
                                        #===============================================================================
## mycoords  <- SpatialPoints(cbind(Y$position_x, Y$position_y))
k         <- 5
mesh      <- inla.mesh.2d(mycoords.CV, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)

X <- X.train
Y <- Y.train

##
## size of discretized fields 
p           <- mesh$n

## circular mesh
mesh.hd <- inla.mesh.1d(seq(0, 2*pi, len=30), boundary="cyclic", degree=1)
p.hd    <- mesh.hd$n

## line splits
Ypos.ls      <- split.segments.wrapper.function(X=X, mesh=mesh, mesh.hd =mesh.hd)
                                        #Ypos         <- Ypos.tmp.ls$Ypos     # guess this is an error; shoul dbe Ypos.ls
                                        #filter.index <- Ypos.tmp.ls$filter.index
Ypos  <- Ypos.ls$Ypos 
filter.index <- Ypos.ls$filter.index

## ------------------
## Integration points
## ------------------
## 
## CHECK always that: dim(coords.trap) == length(HD.data); length(HD.data) == length(T.data)
##
coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))

## df.path.split <- data.frame(x.coord = coords.trap[,1], y.coord = coords.trap[,2], hd = HD.data, time = T.data)

## mesh for temporal process
## print(paste("mesh1d:", head(diff(mesh1d$loc))))
mesh1d  <- inla.mesh.1d(loc=c(T.data[seq(1, length(T.data), by = 300)], T.data[length(T.data)]), order=2)

## Matrix of basis function evaluations for positional data (INTEGRAL term of Poisson likelihood)
Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
A      <- inla.row.kron(Ahd, Aosc)

## ## Matrix of basis function evaluations for observed firing events (SUM term of Poisson likelihood)
Atildeobs <- inla.spde.make.A(mesh=mesh1d, Y$firing_times)
Aosc.obs  <- inla.spde.make.A(mesh=mesh, loc=as.matrix(Y %>% dplyr:: select(position_x, position_y)))
Ahd.obs   <- inla.spde.make.A(mesh=mesh.hd, Y$hd)
Aobs      <- inla.row.kron(Ahd.obs, Aosc.obs)

## Line segment lengths and Time interval lengths
## ----------------------------------------------

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

## !!
df.prism.M1_M2 <- df.prism.M1.M2.wrapper2(At.indices= At.indices, Aosc.indices=Aosc.indices, A.indices=A.indices,
                                          dGamma=dGamma, T.data=T.data, HD.data=HD.data, coords.trap=coords.trap)
df.prism.M1                      <- df.prism.M1_M2 %>% dplyr::select(-c(val.M2.space.time, val.M2.space.direction.time)) %>% unnest(cols=c(val.M1))
df.prism.M2.space.time           <- df.prism.M1_M2 %>% dplyr::select(-c(val.M1, val.M2.space.direction.time)) %>% unnest(cols=c(val.M2.space.time))
df.prism.M2.space.direction.time <- df.prism.M1_M2 %>% dplyr::select(-c(val.M1, val.M2.space.time)) %>% unnest(cols=c(val.M2.space.direction.time))
    

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
    summarize(val = sum(pmax(dGamma.trap*val.M0, tol))/2,
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

## 
## Finally, the W.ipoints.M0 matrix is created below which is in format appropriate to be used in lgcp ipoints arg
## 
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
    summarize(val = sum(pmax(dGamma.trap*val.M1, tol))/2,
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


## space-time
df.W.M2.space.time <- rbind(df.prism.M2.space.time %>% mutate(group=tk, dGamma.lag=0) %>%
                 dplyr::select(group, time, coords, dGamma, dGamma.lag, l, i, val.M2),
                 df.prism.M2.space.time %>% 
                 filter(tk!=1) %>%
                 mutate(time=time.lag, coords=coords.lag,
                        group=tk-1,
                        dGamma=0) %>%
                 dplyr::select(group, time, coords, dGamma, dGamma.lag, l, i, val.M2)) %>%
    arrange(group) %>%
    mutate(dGamma.trap = dGamma + dGamma.lag) ## %>%

tol <- 0
df.dGamma.sum.k.kplus1.M2.space.time <- df.W.M2.space.time %>% group_by(group, l, i) %>%
    summarize(val = sum(pmax(dGamma.trap*val.M2, tol)),
              time = unique(time),
              ## direction=unique(direction),
              coords=unique(coords))

W <- sparseMatrix(i=df.dGamma.sum.k.kplus1.M2.space.time$l,
                  j=df.dGamma.sum.k.kplus1.M2.space.time$i,
                  x=df.dGamma.sum.k.kplus1.M2.space.time$val/2)

W <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(Aosc)-ncol(W)))

## 
## Finally, the W.ipoints.M2.space.direction matrix is created below which a format appropriate to be used in inlabru
## 
W.ipoints.M2.space.time <- as(W, "dgTMatrix")
W.ipoints.M2.space.time <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2.space.time@i+1],
                                      coords.x1 = mesh$loc[W.ipoints.M2.space.time@j+1,1],
                                      ## mapindex2space.direction_basis(W.ipoints.M2.space.time@j+1)[,2]
                                      coords.x2 =mesh$loc[W.ipoints.M2.space.time@j+1,2],
                                      ## mapindex2space.direction_basis(W.ipoints.M2.space.time@j+1)[,3]
                                      weight=W.ipoints.M2.space.time@x) %>% arrange(firing_times)


## current implementation for M_{Omega,T} uses segmentation of path
## based on space-direction-time prism. This results in repeated
## coords and times with different weights
## Q: is 


## space-direction-time
## space-direction-time


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
    ## select(-c(dGamma, dGamma.lead, dGamma.lag))

## this matrix is formed so that we can compute the weight in the trapezoidal rule associated with
## one line/time/arc segment. Consider for the sake of simplicity the weight associated with the first line/time/arc segment.
## This has the form: (Li denoting the length of the ith line segment which is stored in dGamma)
## (L_1/2) * sum_{k*=0}^{1} [ psi_i (s(t_k*)) * psi_j (theta(t_k*)) * psi_l (t_k*) ],
## and for the second:
## (L_2/2) * sum_{k*=1}^{2} [ psi_i (s(t_k*)) * psi_j (theta(t_k*)) * psi_l (t_k*) ],
## etc.
## this operation is done next by

tol <- 0
## df.dGamma.sum.k.kplus1.M2 <- df.W.M2 %>% group_by(group, l, i) %>%
##     summarize(val = sum(max(dGamma.trap*val.M2, tol)))

df.dGamma.sum.k.kplus1.M2 <- df.W.M2 %>% group_by(group, l, i) %>%
    summarize(val = sum(pmax(dGamma.trap*val.M2, tol)),
              time = unique(time),
              direction=unique(direction),
              coords=unique(coords))

## where an additional tol argument was used to avoid numerical instability in cases
## where dGamma.trap * val.M2 was negative. The latest implementation does not suffer
## from any instability so tol is set to 0. Note that dGamma.trap is always equal to dGamma
## L_k for a weight term that comprises a sum of basis functions at the kth and (k+1)th line/time/arc segment as
## the formula above suggests. 
## There must be a cleaner way of calculating this but have given up.
## A snapshot of the data frame above is given below
## 
## > head(as.data.frame(df.dGamma.sum.k.kplus1.M2),50)
##    group l    i          val.M2        time direction  coords.1  coords.2
## 1      1 1 7660 1.056309e-02 0.003131200  1.365427  54.58986 101.65842
## 2      1 1 7726 1.065576e-01 0.003131200  1.365427  54.58986 101.65842
## 3      1 1 8037 5.540326e-02 0.003131200  1.365427  54.58986 101.65842
## 4      1 1 8652 1.272533e-16 0.003131200  1.365427  54.58986 101.65842
## 5      1 1 8932 4.572901e-03 0.003131200  1.365427  54.58986 101.65842
## 6      1 1 8998 3.970761e-02 0.003131200  1.365427  54.58986 101.65842
## 7      1 1 9309 2.367703e-02 0.003131200  1.365427  54.58986 101.65842
## 8      1 1 9924 4.741966e-17 0.003131200  1.365427  54.58986 101.65842
## 9      1 2 7660 0.000000e+00 0.003131200  1.365427  54.58986 101.65842
## 10     1 2 7726 1.284553e-04 0.003131200  1.365427  54.58986 101.65842
## 11     1 2 8037 6.678871e-05 0.003131200  1.365427  54.58986 101.65842
## 12     1 2 8652 1.534040e-19 0.003131200  1.365427  54.58986 101.65842
## 13     1 2 8932 0.000000e+00 0.003131200  1.365427  54.58986 101.65842
## 14     1 2 8998 4.786758e-05 0.003131200  1.365427  54.58986 101.65842
## 15     1 2 9309 2.488814e-05 0.003131200  1.365427  54.58986 101.65842
## 16     1 2 9924 5.716447e-20 0.003131200  1.365427  54.58986 101.65842
## 17     2 1 7726 4.774293e-01 0.008852429  1.358788  54.36973 101.69133
## 18     2 1 8037 2.303666e-01 0.008852429  1.358788  54.36973 101.69133
## 19     2 1 8652 3.098200e-01 0.008852429  1.358788  54.36973 101.69133
## 20     2 1 8998 1.651041e-01 0.008852429  1.358788  54.36973 101.69133
## 21     2 1 9309 8.584379e-02 0.008852429  1.358788  54.36973 101.69133
## 22     2 1 9924 5.214417e-02 0.008852429  1.358788  54.36973 101.69133
## 23     2 2 7726 2.983597e-03 0.008852429  1.358788  54.36973 101.69133
## 24     2 2 8037 2.777072e-04 0.008852429  1.358788  54.36973 101.69133
## 25     2 2 8652 1.936157e-03 0.008852429  1.358788  54.36973 101.69133
## 26     2 2 8998 5.021535e-04 0.008852429  1.358788  54.36973 101.69133
## 27     2 2 9309 1.034848e-04 0.008852429  1.358788  54.36973 101.69133
## 28     2 2 9924 3.258644e-04 0.008852429  1.358788  54.36973 101.69133
## 29     3 1 7643 1.282249e-03 0.032641279  1.331181  53.45440 101.82815
## 30     3 1 7726 3.378735e-02 0.032641279  1.331181  53.45440 101.82815
## 31     3 1 8037 1.021394e-17 0.032641279  1.331181  53.45440 101.82815
## 32     3 1 8652 2.193160e-02 0.032641279  1.331181  53.45440 101.82815
## 33     3 1 8915 2.001912e-04 0.032641279  1.331181  53.45440 101.82815
## 34     3 1 8998 5.686572e-03 0.032641279  1.331181  53.45440 101.82815
## 35     3 1 9309 1.719054e-18 0.032641279  1.331181  53.45440 101.82815
## 36     3 1 9924 3.690209e-03 0.032641279  1.331181  53.45440 101.82815
## 37     3 2 7643 8.473320e-06 0.032641279  1.331181  53.45440 101.82815
## 38     3 2 7726 2.185065e-04 0.032641279  1.331181  53.45440 101.82815
## 39     3 2 8037 6.382991e-20 0.032641279  1.331181  53.45440 101.82815
## 40     3 2 8652 1.449277e-04 0.032641279  1.331181  53.45440 101.82815
## 41     3 2 8915 1.322897e-06 0.032641279  1.331181  53.45440 101.82815
## 42     3 2 8998 3.553707e-05 0.032641279  1.331181  53.45440 101.82815
## 43     3 2 9309 1.074288e-20 0.032641279  1.331181  53.45440 101.82815
## 44     3 2 9924 2.306121e-05 0.032641279  1.331181  53.45440 101.82815
## 45     4 1 7643 3.608168e-01 0.034324800  1.329228  53.38962 101.83783
## 46     4 1 7726 5.672905e-01 0.034324800  1.329228  53.38962 101.83783
## 47     4 1 8652 4.003107e-01 0.034324800  1.329228  53.38962 101.83783
## 48     4 1 8915 6.673838e-02 0.034324800  1.329228  53.38962 101.83783
## 49     4 1 8998 8.856825e-02 0.034324800  1.329228  53.38962 101.83783
## 50     4 1 9924 7.404336e-02 0.034324800  1.329228  53.38962 101.83783


## once we have all such weights, we aggregate them using
## sum_{k=1}^N (L_k/2) * sum_{k*=k}^{k+1} [ psi_i (s(t_k*)) * psi_j (theta(t_k*)) * psi_l (t_k*) ]
## this is done efficiently with sparseMatrix

W <- sparseMatrix(i=df.dGamma.sum.k.kplus1.M2$l,
                  j=df.dGamma.sum.k.kplus1.M2$i,
                  x=df.dGamma.sum.k.kplus1.M2$val/2)

## if there are columns in the matrix W above that are everywhere 0, then the sparseMatrix function drops them.
## I think this may happen due to the knots placed outside the domain 
## where the animal never enters (cf a plot of spatial mesh - outer boundary). In this case, there is no contribution as there is
## no line segment in any  these triangles so I manually add them by appending 0 columns in the matrix W

W         <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(A)-ncol(W)))

## 
## Finally, the W.ipoints.M2 matrix is created below which a format appropriate to be used in inlabru
## 
W.ipoints.M2 <- as(W, "dgTMatrix")
W.ipoints.M2 <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2@i+1],
                           hd=mapindex2space.direction_basis(W.ipoints.M2@j+1)[,1],
                           coords.x1 =mapindex2space.direction_basis(W.ipoints.M2@j+1)[,2],
                           coords.x2 =mapindex2space.direction_basis(W.ipoints.M2@j+1)[,3],
                           weight=W.ipoints.M2@x) %>% arrange(firing_times)


df.W.M2.space.direction.time <- rbind(df.prism.M2.space.direction.time %>% mutate(group=tk, dGamma.lag=0) %>%
                 dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val.M2),
                 df.prism.M2.space.direction.time %>% 
                 filter(tk!=1) %>%
                 mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                        group=tk-1,
                        dGamma=0) %>%
                 dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, l, i, val.M2)) %>%
    arrange(group) %>%
    mutate(dGamma.trap = dGamma + dGamma.lag) ## %>%

tol <- 0
df.dGamma.sum.k.kplus1.M2.space.direction.time <- df.W.M2.space.direction.time %>% group_by(group, l, i) %>%
    summarize(val = sum(pmax(dGamma.trap*val.M2, tol)),
              time = unique(time),
              direction=unique(direction),
              coords=unique(coords))

W <- sparseMatrix(i=df.dGamma.sum.k.kplus1.M2.space.direction.time$l,
                  j=df.dGamma.sum.k.kplus1.M2.space.direction.time$i,
                  x=df.dGamma.sum.k.kplus1.M2.space.direction.time$val/2)

W         <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(A)-ncol(W)))

## 
## Finally, the W.ipoints.M2.space.direction matrix is created below which a format appropriate to be used in inlabru
## 
W.ipoints.M2.space.direction.time <- as(W, "dgTMatrix")
W.ipoints.M2.space.direction.time <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2.space.direction.time@i+1],
                                                hd=mapindex2space.direction_basis(W.ipoints.M2.space.direction.time@j+1)[,1],
                                                coords.x1 =mapindex2space.direction_basis(W.ipoints.M2.space.direction.time@j+1)[,2],
                                                coords.x2 =mapindex2space.direction_basis(W.ipoints.M2.space.direction.time@j+1)[,3],
                                                weight=W.ipoints.M2.space.direction.time@x) %>% arrange(firing_times)

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
initial.direction <- list(theta4=log(pi), theta5=0)
initial.time      <- list(theta1=log(100), theta2=log(3))
l = -0.98
u = 1

Pl.Omega       <- Polygon(expand.grid(c(min(X$position_x),max(X$position_x)), c(min(X$position_y),max(X$position_y)))[c(1,2,4,3),])
ID.Omega       <- "Omega"
Pls.Omega      <- Polygons(list(Pl.Omega), ID=ID.Omega)
SPls.Omega     <- SpatialPolygons(list(Pls.Omega))
weights.domain <- ipoints(domain=mesh, samplers=SPls.Omega)
locs           <- weights.domain@coords
rownames(locs) <- NULL


## ----------------
## Fitting M0 model
## ---------------- 
source("rgeneric_models.R")

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
                                           initial.space=initial.space)

A.spatial.field_constr <- inla.spde.make.A(mesh=mesh, loc=locs,
                                           weights=weights.domain@data[,1],
                                           block=rep(1, nrow(weights.domain@coords)))
cmp.space <- firing_times ~
    spde2(cbind(coords.x1, coords.x2), model=space.rgeneric, mapper=bru_mapper(mesh,indexed=TRUE)) + Intercept

fit.space <- lgcp(cmp.space,
                  data = Y.spdf,
                  ips     = W.ipoints.M0,
                  domain  = list(firing_times = mesh1d),
                  options = list(num.threads=8,verbose = TRUE, bru_max_iter=1) )

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
    Intercept

fit.space.direction <- lgcp(cmp.space.direction,
                            data    = Y.spdf,
                            ips     = W.ipoints.M1,
                            domain  = list(firing_times = mesh1d),
                            options = list(num.threads=8, verbose = FALSE, bru_max_iter=1)) %>%
    bru_rerun()


summary(fit.space)
summary(fit.space.direction)


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

## A.spatial.field_constr_along.directions     <- as.matrix(kronecker(Diagonal(mesh.hd$n),
##                                                                    as.matrix(inla.spde.make.A(mesh=mesh, loc=locs,
##                                                                                               weights=weights.domain@data[,1],
##                                                                                               block=rep(1, nrow(weights.domain@coords))), nrow=1)))


cmp.space.time <- firing_times ~
    spde2(cbind(coords.x1, coords.x2), model=space.rgeneric, mapper=bru_mapper(mesh,indexed=TRUE)
          ## ,
         ##  extraconstr=list(A=as.matrix(A.spatial.field_constr,nrow=1), e=0)
          ) + 
    spde1(firing_times, model=time.rgeneric, mapper=bru_mapper(mesh1d, indexed=TRUE)) + Intercept

fit.space.time <- lgcp(cmp.space.time, data = as.data.frame(Y.spdf),
                       ips=W.ipoints.M2.space.time,
                       domain = list(firing_times = mesh1d),
                       options=list(
                           num.threads=8,
                           verbose = TRUE, bru_max_iter=1)) %>%
    bru_rerun()



##
## ----------------------------------
## Fitting space-direction-time model
## ----------------------------------

cmp.space.direction.time <- firing_times ~
    spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
          mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE)))) +
    time(firing_times, model=time.rgeneric, mapper=bru_mapper(mesh1d, indexed=TRUE)) + Intercept

fit.space.direction.time <- lgcp(cmp.space.direction.time, data = as.data.frame(Y.spdf),
                                 ips=W.ipoints.M2,
                                 domain = list(firing_times = mesh1d),
                                 options=list(
                                     num.threads=8,
                                     verbose = TRUE, bru_max_iter=1))



## ----------------------------------------------------------------
## Predictions for the expected number of firing events in test set
## ----------------------------------------------------------------
## calculation of integration weight for connected line segments
## ("clumps"); note that adjacent clumps in test data are joined
## together so that the number of test clumps is less than the
## original training/test split

## unique(Y.test$index.CV) 
## =  1  2  3  4  6  8 11 13 14 15 16 19 20 21 22 23 24 29 33 37 38 43 44 46 47 49 55 56 57 61
## length(list(1:4, c(6), c(8), c(11), 13:16, 
##       19:24, c(29), c(33), 37:38, 43:44, 46:47, c(49), 55:57, c(61))) = 14
## i.e. there are 14 test clumps after grouping adjacent selected clumps
## so for example the first clump in Y.test can be obtained by filtering for index.CV %in% 1:4
weights_test <- weights_line_segments_in_train(X.test=X.test, Y.test = Y.test, mesh=mesh, mesh.hd=mesh.hd, mesh1d=mesh1d)
## weights_test[[1]]$index.CV 1 6 8 11 13 19 29 33 37 43 46 49 55 60
## which is of length 14 I think 1:4 from unique(Y.test$index.CV) are
## now just labelled as 1 here etc..

## use split function to take unique(Y.test$index.CV) and returns a list of consecutive integers

clumps       <- split(unique(Y.test$index.CV), cumsum(c(1, diff(unique(Y.test$index.CV)) != 1)))
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

pred.mean.var <- function(weights.mat, post.sample){
                                        # posterior sample
    n <- length(post.sample) 
    Int <- lapply(1:n, function(i) post.sample[[i]]$Intercept) # simulated intercepts
    GRF <- lapply(1:n, function(i) post.sample[[i]]$spde2)     # simulated GRFs 
    cond.mean <- sapply(1:n, function(i) exp(Int[[i]])*sum(exp(GRF[[i]])*weights.mat)) # posterior means conditioned on latent parameters
    post.mean <- mean(cond.mean)
    post.var  <- mean(cond.mean) + var(cond.mean)   
    list(post.mean=post.mean, post.var=post.var)
}

##===============================================================================
##===============================================================================
## test it on segment/clump 1
##samp1 <- generate(object = fit.space,
##                  n.samples = 100)
##test1 <- pred.mean.var(weights.mat=weights_test$df.W.M0$W.ipoints.M0[[1]]$W.M0,
## post.sample=samp1)
##test1 
##===============================================================================
##===============================================================================
## repeat test on all segments/clumps 
## posterior samples;
samp1 <- generate(object = fit.space, n.samples = 1000)

clumps.mean.var.M0 <- lapply(1:length(clumps), function(i) pred.mean.var(weights.mat=weights_test$df.W.M0$W.ipoints.M0[[i]]$W.M0,
                                                                         post.sample=samp1))
pred.means.M0 <- sapply(clumps.mean.var.M0, function(i) i[[1]]) # predictive means on each segment/clump
pred.vars.M0  <- sapply(clumps.mean.var.M0, function(i) i[[2]]) # .......... variance ..........
                                        # compare means with observed firings
obs.firings <- sapply(1:length(clumps), function(i) nrow(Y.test %>% filter(index.CV %in% clumps[[i]])))

save(pred.means.M0, file="/exports/eddie/scratch/s0233535/pred.means.M0.RData")
save(pred.vars.M0, file="/exports/eddie/scratch/s0233535/pred.vars.M0.RData")
save(obs.firings, file="/exports/eddie/scratch/s0233535/obs.firings.RData")

## repeat for M1
samp2 <- generate(object = fit.space.direction, 
                  n.samples = 1000)

clumps.mean.var.M1 <- lapply(1:length(clumps), function(i) pred.mean.var(weights.mat=weights_test$df.W.M1$W.ipoints.M1[[i]]$W.M1,
                                                                         post.sample=samp2))
pred.means.M1 <- sapply(clumps.mean.var.M1, function(i) i[[1]]) # predictive means on each segment/clump
pred.vars.M1  <- sapply(clumps.mean.var.M1, function(i) i[[2]]) 

save(pred.means.M1, file="/exports/eddie/scratch/s0233535/pred.means.M1.RData")
save(pred.vars.M1,  file="/exports/eddie/scratch/s0233535/pred.vars.M1.RData")


## repeat for space-time model
samp.space.time <- generate(object = fit.space.time, 
                            n.samples = 1000)

clumps.mean.var.space.time <- lapply(1:length(clumps), function(i) pred.mean.var(weights.mat=weights_test$df.W.M1$W.ipoints.M1[[i]]$W.M1,
                                                                         post.sample=samp.space.time))
pred.means.M1 <- sapply(clumps.mean.var.M1, function(i) i[[1]]) # predictive means on each segment/clump
pred.vars.M1  <- sapply(clumps.mean.var.M1, function(i) i[[2]]) 




## 
plot(obs.firings, pred.means.M0, pch=16,
     ylim=c(0,max(c(pred.means.M0, pred.means.M1, obs.firings))))
points(obs.firings, pred.means.M1, pch=1)
abline(a=0,b=1)

plot(log(obs.firings), log(pred.means.M0), pch=16, cex=0.8, 
     ylim=c(0,max(c(log(pred.means.M0), log(pred.means.M1), log(obs.firings)))))
points(log(obs.firings), log(pred.means.M1), pch=1, cex=1.1)
abline(a=0,b=1)

plot(obs.firings, abs(obs.firings - pred.means.M0), pch=16)
points(obs.firings, abs(obs.firings - pred.means.M1), pch=1)
## abline(a=0,b=1)



cbind(obs.firings, pred.means.M0, pred.means.M1)
cbind(obs.firings, abs(pred.means.M0-obs.firings), abs(pred.means.M1-obs.firings))
## compare observed and expected firings
## obs.firings

## floor(pred.means)
## squared error score per segment/clump
se.score.M0 <- (pred.means.M0 - obs.firings)^2
se.score.M1 <- (pred.means.M1 - obs.firings)^2
## Dawid-Sebastiani score
ds.score.M0 <- se.score.M0 / pred.vars.M0 + log(pred.vars.M0)
ds.score.M1 <- se.score.M1 / pred.vars.M1 + log(pred.vars.M1)
## -----------------------------------------------------------------------------

##
S.se    <- se.score.M0-se.score.M1
S.ds    <- ds.score.M0-ds.score.M1
Tobs.se <- mean(S.se)
Tobs.ds <- mean(S.ds)
randT.se <- NULL
randT.ds <- NULL
K       <- 1
while(K <= 10000){
    rand  <- 2*rbinom(length(se.score.M0), 1, 0.5)-1
    randT.se[K] <- mean(rand*S.se)
    randT.ds[K] <- mean(rand*S.ds)
    K <- K+1
}

pval.se <- mean(randT.se > Tobs)
pval.ds <- mean(randT.ds > Tobs)

hist(randT,breaks=50)
abline(v=Tobs, lwd=2)

plot(se.score.M0, se.score.M1)
abline(a=0, b=1)
plot(log(ds.score.M0), log(ds.score.M1))
abline(a=0, b=1)

