set.seed(111086) 
library(tidyverse)
library(purrr)
## ensure latest INLA testing version is installed
## INLA is not on CRAN:
## https://www.r-inla.org/download-install
library(INLA)
library(inlabru)
bru_options_set(inla.mode = "experimental")
## option "experimental" seems to implement Variational Bayes
## correction this is useful for Poisson point process
## likelihoods. This option is implemented in latest INLA testing version
if(inla.pardiso.check() != "SUCCESS: PARDISO IS INSTALLED AND WORKING"){
    ## enable openmp for parallel computing using strategy "huge".
    ## this is going to engage all RAM and core resources of the
    ## computer needs some care when setting this up on Eddie.
     bru_options_set(control.compute = list(openmp.strategy="huge"))
}else{
    bru_options_set(control.compute = list(openmp.strategy="pardiso"))
}
library(sp)
library(fields)
library(nloptr)
library(pals)
inla.setOption(pardiso.license = "/Users/ipapasta/pardiso.lic")
## Unfortunately pardiso license is discontinued for academic users.
source("load_data.R")
source("Functions.R")
source("osc_precision.R")
source("hd_precision.R")
source("temp_precision.R")

k    <- 5
mesh      <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
## plot(mesh, asp=1)
##
## size of discretized fields 
p           <- mesh$n
p.theta     <- 30
theta.nodes <- seq(0, 2*pi, len=p.theta)
mesh.hd     <- inla.mesh.1d(theta.nodes, boundary="cyclic", degree=1)
## x1 <- cos(mesh.hd$loc)
## y1 <- sin(mesh.hd$loc)
## plot(x1, y1)
nodes       <- c(mesh.hd$loc, 2*pi)
intervals   <- head(cbind(nodes, lead(nodes)), -1)
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

Ypos.tmp <- data.frame(
    hd=X$hd, time=X$synced_time,
    coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(coords.lead = lead(coords)) %>%
    mutate(time.lead = lead(X$synced_time)) %>%
    mutate(hd.lead = lead(X$hd)) %>%
    head(-1)



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

line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
                             filter.zero.length=FALSE,
                             ep=Ypos.tmp$coords.lead, tol=.0)


df <- data.frame(origin=line.segments$split.origin,
                 filter.index=line.segments$filter.index,
                 sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                 ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
    group_by(origin) %>%
    summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
    mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
    mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 


## attribute named _data_ stores length of line segments, time differences and arclengths
Ypos <- inner_join(Ypos.tmp %>%
                   mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
    mutate(Li = map2(sp, ep,
                     function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
    mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
    mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 

filter.index  <- do.call("c", Ypos$filter.index)


## ------------------
## Integration points
## ------------------
## 
## CHECK always that: dim(coords.trap) == length(HD.data); length(HD.data) == length(T.data)
##
coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))





## ---------------------------------------
## SpatialPointsDataFrame and SpatialLines
## ---------------------------------------
## inlabru accepts below formats for the data
## Y.spdf is Y data frame except that coords are encoded as SpatialPoints
## Ypos.sldf Ypos data frame except for coords and coords.lead are encoded as SpatialLines
Y.spdf    <- SpatialPointsDataFrame(coords = SpatialPoints(cbind(Y$position_x, Y$position_y)),
                                    data   = as.data.frame(Y%>%dplyr::select(-c(position_x, position_y))))
Ypos.sldf <- SpatialLinesDataFrame(sl   = SpatialLines(lapply(as.list(1:nrow(Ypos)),
                                                              function(k) Lines(list(Line(cbind(c(Ypos$coords[k,1],
                                                                                                  Ypos$coords.lead[k,1]),
                                                                                                c(Ypos$coords[k,2],
                                                                                                  Ypos$coords.lead[k,2])))), ID=k))),
                                   data = Ypos %>% dplyr::select(-c(coords, coords.lead)))

data <- list(Ypos=Ypos, Y=Y, Yspdf=Y.spdf, Ypos.sldf = Ypos.sldf)


## mesh for temporal process
## print(paste("mesh1d:", head(diff(mesh1d$loc))))
mesh1d  <- inla.mesh.1d(loc=c(T.data[seq(1, length(T.data), by = 300)], T.data[length(T.data)]), order=2)


## Matrix of basis function evaluations for positional data (INTEGRAL term of Poisson likelihood)
Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
A      <- inla.row.kron(Ahd, Aosc)


## ## Matrix of basis function evaluations for observed firing events (SUM term of Poisson likelihood)
Atildeobs <- inla.spde.make.A(mesh=mesh1d, data$Y$firing_times)
Aosc.obs  <- inla.spde.make.A(mesh=mesh, loc=as.matrix(data$Y %>% dplyr:: select(position_x, position_y)))
Ahd.obs   <- inla.spde.make.A(mesh=mesh.hd, data$Y$hd)
Aobs      <- inla.row.kron(Ahd.obs, Aosc.obs)


## nrow(A)==nrow(At); 92264 each row of above matrices contains
## non-zero values at knots wrapping a distinct line segment.
## Ck <- sapply(dGamma, function(x) rep(x, 6))
## (coords.trap, HD.data, T.data)


## Line segment lengths and Time interval lengths
## ----------------------------------------------
## NOTE: head direction arclengths are not used since dGamma(t)/d(t) is defined as:
## ((d x(t)/dt)^2 + (d y(t)/dt)^2)^(1/2)
## 
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

## 
## dim(A)[1] == dim(Atilde)[1] is TRUE, both matrices are basis function evaluations at
## the starting coordinates (A) and the starting times (Atilde) of the line segments 
## (positional data), that is, rows are line/time/arc segments (properly speaking: line segments, time intervals and arcs?)
## for every starting coord/time/angle of a  line/time/arc segment, there are 6 spatio-temporal basis functions that give a non-zero contribution,
## that is, 3*2 (3 knots of a triangle * 2 knots of an arc) and 2 temporal basis functions that give non-zero contributions,
## that is, 2 time interval knots.
## A.indices and At.indices: first column is renamed to tk which stands for index of line/time/arc segment
## Each tk appears 6 times, e.g., length(A.indices[,1])/6 = N, where N is the number of line/time/arc segments
## A.indices: second column is renamed to i which stands for the index of the spatio-directional basis function
## At.indices: second row is l which stands for the index of the temporal basis function
##
names(Aosc.indices) <- c("tk", "i", "psi.o") #ot: omega 
names(A.indices) <- c("tk", "i", "psi.ot") #ot: omega x theta
names(At.indices) <- c("tk", "l", "psi.t") 



## the code below used to define df.prism.M1_M2 first groups At.indices and A.indices by line/time/angle segment
## then nests them so each row of the nested data frame contains all basis function evaluation data for line/time/arc segments.
## During nesting, information on the index of the basis functions and its associated value is stored
## in new column variables named as
## data.x: for the temporal basis functions, and as
## data.y for the spatio.directional basis functions.
## Then information on times, head directions, and coordinates is appended to the nested data frame:
## for every line/time/arc segment (indexed by variable tk), information about the coordinate, the time and the head direction is appended.
## Information on lags and leads is also included because these are subsequently used to get the weights for
## the numerical integration scheme based on the trapezoidal rule. For the computation of the weights, the
## lengths of the line segments (dGamma) together with their lags and leads are also appended.
## Finally, for the code in the last column variable named _val_. Fix a line segment, say the first one tk=1.
## Then, for example, the first elements of the column variables data.x and data.y (these are lists due to the nest operation) are:
## 
## > data.y[[1]]
## # A tibble: 6 x 2
##       i psi.ot
##   <dbl>  <dbl>
## 1  7660 0.0475
## 2  7726 0.405 
## 3  8037 0.246 
## 4  8932 0.0205
## 5  8998 0.175 
## 6  9309 0.106
## 
## > data.x[[1]]
##       l psi.t
##   <dbl> <dbl>
## 1     1     1
## 2     2     0
##

## 
## The code in this last function creates the Cartesian product {1,2} X {7660, 7726, 8037, 8932, 8998, 9309},
## that is, the set of all ordered pairs (a,b) where a \in {1,2} and b \in {7660, 7726, 8037, 8932, 8998, 9309}
## with expand.grid, and calculates, for each pair, the product of the basis functions \psi_{T} * \psi_{Omega x Theta}
## Lastly, the function returns a data.frame that has the index of the temporal basis function l, the index of the
## spatio-directional basis function i, and the product of the basis functions _val_
## this data framed is stored in the column variable named val. The final commands discard data.x and data.y
## which are no longer used and unnest the data frame to bring it back to standard form
## 

## df.prism <- full_join(At.indices %>% group_by(tk) %>% nest(),
##                 A.indices %>% group_by(tk) %>% nest(), by="tk") %>%
##     arrange(tk) %>%
##     ungroup %>% 
##     mutate(
##         time          = T.data,
##         time.lag      = c(0, time[-length(time)]), #issue with NAs, use 0 (will be discarded later)
##         direction     = HD.data,
##         direction.lag = c(0, HD.data[-length(direction)]),
##         coords = I(coords.trap),
##         coords.lag = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
##         dGamma=c(dGamma,0),
##         dGamma.lead = lead(dGamma),
##         dGamma.lag = lag(dGamma),
##         val.M1 = pmap(list(data.y, dGamma), function(x, y, z){
##             ooo <- unlist(lapply(1:nrow(y)), function(k) {
##                 y$psi.ot[oo[k,2]]})
##             oooo <- data.frame(i=y$i[oo[,2]],val=ooo)
##             oooo
##         }),
##         val = pmap(list(data.x, data.y, dGamma), function(x, y, z){
##             oo  <- expand.grid(1:nrow(x), 1:nrow(y))
##             ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
##                 x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
##             oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val=ooo)
##             oooo
##         })) %>%
##     dplyr::select(-c("data.x", "data.y")) %>%
##     unnest(val)

df.prism.M0 <- Aosc.indices %>% group_by(tk) %>% nest() %>%
    arrange(tk) %>%
    ungroup %>% 
    mutate(
        time          = T.data,
        time.lag      = c(0, time[-length(time)]),
        direction     = HD.data,
        direction.lag = c(0, HD.data[-length(direction)]),
        coords        = I(coords.trap),
        coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
        dGamma=c(dGamma,0),
        dGamma.lead = lead(dGamma),
        dGamma.lag = lag(dGamma),
        val.M0 = pmap(list(data, dGamma), function(y, z) {
            oo  <- 1:nrow(y)
            ooo <- unlist(lapply(1:nrow(y), function(k) {
                y$psi.o[oo[k]]}))
            oooo <- data.frame(i=y$i[oo],val.M0=ooo)
            oooo
        })) %>% 
    dplyr::select(-c("data")) 

df.prism.M1_M2 <- full_join(At.indices %>% group_by(tk) %>% nest(),
                A.indices %>% group_by(tk) %>% nest(), by="tk") %>%
    arrange(tk) %>%
    ungroup %>% 
    mutate(
        time          = T.data,
        time.lag      = c(0, time[-length(time)]),
        direction     = HD.data,
        direction.lag = c(0, HD.data[-length(direction)]),
        coords        = I(coords.trap),
        coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
        dGamma=c(dGamma,0),
        dGamma.lead = lead(dGamma),
        dGamma.lag = lag(dGamma),
        val.M1 = pmap(list(data.y, dGamma), function(y, z) {
            oo  <- 1:nrow(y)
            ooo <- unlist(lapply(1:nrow(y), function(k) {
                y$psi.ot[oo[k]]}))
            oooo <- data.frame(i=y$i[oo],val.M1=ooo)
            oooo
        }),
        val.M2 = pmap(list(data.x, data.y, dGamma), function(x, y, z){
            oo  <- expand.grid(1:nrow(x), 1:nrow(y))
            ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
            oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
            oooo
        })) %>%
    dplyr::select(-c("data.x", "data.y")) 

df.prism.M0 <- df.prism.M0 %>% unnest(cols=c(val.M0))
df.prism.M1 <- df.prism.M1_M2 %>% dplyr::select(-val.M2) %>% 
    unnest(cols=c(val.M1))
df.prism.M2 <- df.prism.M1_M2 %>% dplyr::select(-val.M1) %>%
    unnest(cols=c(val.M2))



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
## df.dGamma.sum.k.kplus1.M0 <- df.W.M0 %>% group_by(group, i) %>%
##     summarize(val = sum(max(dGamma.trap*val.M0, tol)))
## tmp <- df.W.M0 %>% group_by(group, i) %>% nest

df.dGamma.sum.k.kplus1.M0 <- df.W.M0 %>% group_by(group, i) %>%
    summarize(val = sum(max(dGamma.trap*val.M0, tol))/2,
              time = unique(time),
              direction=unique(direction),
              coords=unique(coords))  %>%
    ungroup %>% group_by(i) %>%
    summarize(val = sum(val))

## fit.space.direction$summary.random$spde2$mean
    

W.M0 <- sparseVector(i=df.dGamma.sum.k.kplus1.M0$i,
                     x=df.dGamma.sum.k.kplus1.M0$val,
                     length=mesh$n)

W.M0.vector <- sparseVector(i=df.dGamma.sum.k.kplus1.M0$i,
                            x=df.dGamma.sum.k.kplus1.M0$val,
                            length=mesh$n)

## 
## Finally, the W.ipoints.M2 matrix is created below which is in format appropriate to be used in lgcp ipoints arg
## 
W.ipoints.M0 <- as(W.M0, "sparseMatrix")
W.ipoints.M0 <- data.frame(coords.x1 = mesh$loc[W.ipoints.M0@i+1,1],
                           coords.x2 = mesh$loc[W.ipoints.M0@i+1,2],
                        weight=W.ipoints.M0@x) 


## ------------------------------------------------
## M1 model:  Integration weights integration knots 
## ------------------------------------------------
## for ips argument in lgcp function, we need to have the integration weights in the following format
## !! 4th Feb for M1 firing_times are not needed as there is no temporal basis function? Validate
## If firing_times are kept, does code run?
## |--------------+-----------+-----------+-----------+--------|
## | firing_times | hd        | coords.x1 | coords.x2 | weight |
## |--------------+-----------+-----------+-----------+--------|
## | t_{j(k)}     | hd_{i(k)} | s1_{i(k)} | s1_{i(k)} | w_k    |
## | .            | .         | .         | .         | .      |
## | .            | .         | .         | .         | .      |
## | .            | .         | .         | .         | .      |
## 
## 
## In what follows: two copies of the df.prism.M1_M2 are created
## the first copy has all line/time/arc segments and sets dGamma.lag = 0 everywhere
## the second copy removes the first line segment and relabels the line/time/arc segments to start from 1.
## this way, data from adjacent line/time/arc segments are given the same label.
## data are then grouped by line/time/arc segment. dGamma is set to 0 everywhere for the second copy.
## a snapshot of the data frame is shown below

##     group        time direction  coords.1  coords.2     dGamma dGamma.lag       i       val.M1  dGamma.trap
## 1       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    7660 4.745711e-02  0.22258192
## 2       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    7660 4.047041e-01  0.22258192
## 3       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    7726 2.457178e-01  0.22258192
## 4       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    7726 2.054480e-02  0.22258192
## 5       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    8037 1.752017e-01  0.22258192
## 6       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    8037 1.063744e-01  0.22258192
## 7       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    8932 4.793113e-01  0.22258192
## 8       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    8932 2.492118e-01  0.22258192
## 9       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    8998 5.724036e-16  0.22258192
## 10      1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    8998 1.786106e-01  0.22258192
## 11      1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    9309 9.286630e-02  0.22258192
## 12      1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000    9309 2.133004e-16  0.22258192
## 13      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    7726 4.793113e-01  0.22258192
## 14      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    7726 2.492118e-01  0.22258192
## 15      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    8037 5.724036e-16  0.22258192
## 16      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    8037 1.786106e-01  0.22258192
## 17      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    8652 9.286630e-02  0.22258192
## 18      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    8652 2.133004e-16  0.22258192
## 19      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    8998 5.190876e-01  0.22258192
## 20      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    8998 1.569205e-16  0.22258192
## 21      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    9309 3.368535e-01  0.22258192
## 22      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    9309 8.736490e-02  0.22258192
## 23      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    9924 2.641046e-17  0.22258192
## 24      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192    9924 5.669404e-02  0.22258192
##

## 

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


## this matrix is formed so that we can compute the weight in the trapezoidal rule associated with
## one line/time/arc segment. Consider for the sake of simplicity the weight associated with the first line/time/arc segment.
## This has the form: (Li denoting the length of the ith line segment which is stored in dGamma)
## (L_1/2) * sum_{k*=0}^{1} [ psi_i (s(t_k*)) * psi_j (theta(t_k*)) ],
## and for the second:
## (L_2/2) * sum_{k*=1}^{2} [ psi_i (s(t_k*)) * psi_j (theta(t_k*)) ],
## etc.
## this operation is done next by

tol <- 0
## df.dGamma.sum.k.kplus1.M1 <- df.W.M1 %>% group_by(group, i) %>%
##     summarize(val = sum(max(dGamma.trap*val.M1, tol)))
## tmp <- df.W.M1 %>% group_by(group, i) %>% nest

df.dGamma.sum.k.kplus1.M1 <- df.W.M1 %>% group_by(group, i) %>%
    summarize(val = sum(max(dGamma.trap*val.M1, tol))/2,
              time = unique(time),
              direction=unique(direction),
              coords=unique(coords))  %>%
    ungroup %>% group_by(i) %>%
    summarize(val = sum(val))

## fit.space.direction$summary.random$spde2$mean
    

W.M1 <- sparseVector(i=df.dGamma.sum.k.kplus1.M1$i,
                     x=df.dGamma.sum.k.kplus1.M1$val,
                     length=mesh$n * mesh.hd$n)

W.M1.vector <- sparseVector(i=df.dGamma.sum.k.kplus1.M1$i,
                            x=df.dGamma.sum.k.kplus1.M1$val,
                            length=mesh$n * mesh.hd$n)

## 
## Finally, the W.ipoints.M2 matrix is created below which is in format appropriate to be used in lgcp ipoints arg
## 
W.ipoints.M1 <- as(W.M1, "sparseMatrix")
W.ipoints.M1 <- data.frame(hd=mapindex2space.direction_basis(W.ipoints.M1@i+1)[,1],
                           coords.x1 =mapindex2space.direction_basis(W.ipoints.M1@i+1)[,2],
                        coords.x2 =mapindex2space.direction_basis(W.ipoints.M1@i+1)[,3],
                        weight=W.ipoints.M1@x) 



## df.indices labels the indices of the knots for the directional, the spatial and the spatio-directional basis functions
## this data frame is created to create correspondences
## between spatio-directional basis knots with spatial basis knots, an
## between spatio-directional basis knots with head directional basis knots, respectively.



## So for example, if the spatio-directional basis knots are labeled as 1, 2, ..., p_Omega * p_Theta
## then the function mapindex2space.direction_basis takes as argument the label of spatio-diretional basis knot
## and returns the coordinates and the head direction associated with the spatial basis function and the
## head directional basis function. This function uses mapindex2space.direction_basis which works similarly but
## instead of returning coords and angles, it returns the indices of the associated basis functions.




## ------------------------------------------------
## M2 model:  Integration weights integration knots 
## ------------------------------------------------
## NOTE: Suppose there are N line segments
## integration points are (t_i, s_i, hd_i), i=0, .., N where
## t_i =  initial time of line segment i
## s_i =  position at initial time of line segment i
## s_i =  head direction of animal at initial time of line segment i
## |--------------+-----------+-----------+-----------+--------|
## | firing_times | hd        | coords.x1 | coords.x2 | weight |
## |--------------+-----------+-----------+-----------+--------|
## | t_{j(k)}     | hd_{i(k)} | s1_{i(k)} | s1_{i(k)} | w_k    |
## | .            | .         | .         | .         | .      |
## | .            | .         | .         | .         | .      |
## | .            | .         | .         | .         | .      |
##

## 
## In what follows: two copies of the df.prism.M2 are created
## the first copy has all line/time/arc segments and sets dGamma.lag = 0 everywhere
## the second copy removes the first line segment and relabels the line/time/arc segments to start from 1, that is,
## data are grouped by line/time/arc segment. dGamma is set to 0 everywhere for the second copy
## a snapshot of the data frame is shown below

##     group        time direction  coords.1  coords.2     dGamma dGamma.lag l      i       val.M2 dGamma.trap
## 1       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 1   7660 4.745711e-02  0.22258192
## 2       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 2   7660 0.000000e+00  0.22258192
## 3       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 1   7726 4.047041e-01  0.22258192
## 4       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 2   7726 0.000000e+00  0.22258192
## 5       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 1   8037 2.457178e-01  0.22258192
## 6       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 2   8037 0.000000e+00  0.22258192
## 7       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 1   8932 2.054480e-02  0.22258192
## 8       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 2   8932 0.000000e+00  0.22258192
## 9       1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 1   8998 1.752017e-01  0.22258192
## 10      1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 2   8998 0.000000e+00  0.22258192
## 11      1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 1   9309 1.063744e-01  0.22258192
## 12      1 0.003131200  1.365427  54.58986 101.65842 0.22258192 0.00000000 2   9309 0.000000e+00  0.22258192
## 13      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 1   7726 4.787342e-01  0.22258192
## 14      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 2   7726 5.771147e-04  0.22258192
## 15      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 1   8037 2.489118e-01  0.22258192
## 16      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 2   8037 3.000635e-04  0.22258192
## 17      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 1   8652 5.717144e-16  0.22258192
## 18      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 2   8652 6.892025e-19  0.22258192
## 19      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 1   8998 1.783955e-01  0.22258192
## 20      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 2   8998 2.150560e-04  0.22258192
## 21      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 1   9309 9.275448e-02  0.22258192
## 22      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 2   9309 1.118157e-04  0.22258192
## 23      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 1   9924 2.130436e-16  0.22258192
## 24      1 0.003131200  1.365427  54.58986 101.65842 0.00000000 0.22258192 2   9924 2.568244e-19  0.22258192
##

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
    summarize(val = sum(max(dGamma.trap*val.M2, tol)),
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
## be used to obtain the M matrices  (defined in
## spde2_implementation.spdf) which are used both in
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

## some information on lgcp arguments
## domain: means integration domain. When integration domain is supplied together with samplers,
## then lgcp constructs the integration scheme. There are options for
## how it does that. By default, .. it will first take domain and
## samplers, say for example samplers specify subintervals (time),
## then it constructs the intersection of domains and samplers and
## then it will place integration points on the remaining knots, so
## essentially, original knots inside intervals and at the endpoints
## of the intervals.  Then it takes samplers which will have say one
## interval.

## ----------------
## Fitting M0 model
## ----------------
## current implementation below is correct. Integration weights computed directly from ipoints function
## ipoints is invoked since samplers and domain are given

if(FALSE){
    cmp.oscillating.rgeneric <- coordinates ~ 
        spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
        Intercept
    fit.oscillating.rgeneric <- lgcp(cmp.oscillating.rgeneric,
                                     data = Y.spdf,
                                     samplers = Ypos.sldf,
                                     domain = list(coordinates = mesh),
                                     options=list(verbose = TRUE))
}

## below, samplers is not provided and fit is made using from integration weights computed from trapezoidal approximation
## weights are provided via W.ipoints.M0 in ips argument of lgcp. TODO: check the two implementations give
## roughly the same output. 
cmp.space <- firing_times ~
    spde2(cbind(coords.x1, coords.x2), model=oscillating.rgeneric, mapper=bru_mapper(mesh,indexed=TRUE)) + Intercept
fit.space <- lgcp(cmp.space,
                  data = Y.spdf,
                  ips     = W.ipoints.M0,
                  domain  = list(firing_times = mesh1d),
                  options = list( num.threads=8,verbose = TRUE, bru_max_iter=1) )



## fit.space$summary.hyperpar
##                        mean       mean         mean
## Theta1 for spde2  2.9597178  2.9626358    2.9644598
## Theta2 for spde2  0.5297674  0.5308157    0.5306628
## Theta3 for spde2 -2.6339410  -2.7731437  -2.7649695

## ----------------
## Fitting M1 model
## ----------------
## NOTES: the integration points that are supplied are incorrect here but were used to verify that the code runs.
## A correct implementation would need to compute the W.ipoints for M1 correctly-Work in progress
cmp.space.direction <- firing_times ~
    spde2(list(spatial=cbind(coords.x1, coords.x2), direction=hd), model=space.direction.rgeneric,
          mapper=bru_mapper_multi(list(spatial=bru_mapper(mesh,indexed=TRUE), direction=bru_mapper(mesh.hd, indexed=TRUE)))) +
    Intercept

fit.space.direction <- lgcp(cmp.space.direction, data = Y.spdf,
                            ips     = W.ipoints.M1,
                            domain  = list(firing_times = mesh1d),
                            options = list( num.threads=8,verbose = TRUE, bru_max_iter=1) )


## model based estimates of expected number of events on the entire path from M0 and M1
EN.M0 <- exp(fit.space$summary.fixed$mean) * sum(W.M0.vector*exp(fit.space$summary.random$spde2$mean))
EN.M1 <- exp(fit.space.direction$summary.fixed$mean) * sum(W.M1.vector*exp(fit.space.direction$summary.random$spde2$mean))
EN.M0 - nrow(Y)
EN.M1 - nrow(Y)



## ----------------
## Fitting M2 model (computationally expensive)
## ----------------
##
## 

if(FALSE){
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
}


## Model based estimate of expected number of points on the entire path


## exp(fit.space.direction$summary.fixed$mean) * sum(W.M1.vector*exp(fit.space.direction$summary.random$spde2$mean))


