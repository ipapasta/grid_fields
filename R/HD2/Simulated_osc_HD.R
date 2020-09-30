## -------------------------------
## Ogata's method
## -------------------------------

## ------------------------------------------
## simulation from fitted intensity function
## ------------------------------------------

load("mesh.RData")
options(warn=-1)
library(INLA)
library(gridExtra)
library(RColorBrewer)
library(scales)
library(Matrix)
source("Functions.R")
load("/home/ipapasta/Software/R/grid_fields/R/data/Simulated_data/simulated_cn_same_30var_v_grid5trials7_simulated/fitted_model.osc_HD.RData")
counter <- 0
Xestimates<-pthetapc.prop.marg.post_osc_HD (par.theta=opt.theta$par, hyperpar=hyperpar, data=data, X=Xinit, beta=betainit, mesh=mesh, A=A, tA=tA, W=W, tW=tW, DW=DW,
                                Aobs=Aobs, tAobs=tAobs, order.HD=order.HD, type="biased", acc=1e-7, print.verbose=FALSE, return.X=TRUE)
p            <- mesh$n
alpha        <- 2
X     <- Xest
beta  <- betaest
Xbeta <- rBind(t(t(X)), t(t(beta)))      
n         <- nrow(data$Y)
one.n     <- matrix(rep(1, n), ncol=1)
lambdafit <- exp((as.numeric(beta)*one.n) + Aobs%*% Xest)  

##
## Spatial coordinates
##
coords.pos <- do.call("rbind", Ypos$coords)
hd.pos     <- Ypos$hd

## ---------------------------------------------------------
## lambda
## ---------------------------------------------------------
one.vec <- matrix(rep(1, nrow(coords.pos)), ncol=1)
order.HD <- 3
proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords.pos))
Aosc.pred      <- proj.s.pred$proj$A
Ahd.pred       <- circular.make.A(hd.pos, order=order.HD)
Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
lambdapred  <- exp(matrix((as.numeric(beta)*one.vec), ncol=1) + Apred%*%Xest)

## -----------------------------------------------
## simulation of homogeneous poisson point process
## -----------------------------------------------

lambda_star  <- max(lambdapred)
Ti           <- NULL
tot <- 0
T.pos <- Ypos$time

coord.sim <- NULL
T.sim <- NULL
hd.sim <- NULL

## Simulation
while(tot < tail(tail(Ypos$time),1)){
    Ti <- c(Ti, rexp(1, 1/lambda_star)) #interarrival times
    T.star <- tail(cumsum(Ti),1)
    wh.prev <- tail(which(T.star > T.pos) , 1)
    wh.next <- head(which(T.star < T.pos) , 1)
    T.prev  <- T.pos[tail(which(T.star > T.pos) , 1)]
    T.next  <- T.pos[head(which(T.star < T.pos) , 1)]
    coord.star <- unlist(Ypos$coords[wh.prev]) + (T.star-T.prev)/(T.next-T.prev)*(unlist(Ypos$coords[wh.next]) - unlist(Ypos$coords[wh.prev]))
    hd.star    <- unlist(Ypos$hd[wh.prev]) + (T.star-T.prev)/(T.next-T.prev)*(unlist(Ypos$hd[wh.next]) - unlist(Ypos$hd[wh.prev]))
    ## evaluation of intensity function
    proj.s.pred    <- inla.mesh.projector(mesh, loc=matrix(coord.star, nrow=1))
    Aosc.pred      <- proj.s.pred$proj$A
    Ahd.pred       <- circular.make.A(hd.star, order=order.HD)
    Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
    lambdapred  <- exp(matrix((as.numeric(beta)), ncol=1) + Apred%*%Xest)
    U.star      <- runif(1, 0, lambda_star)
    if(U.star < as.numeric(lambdapred)){
        T.sim <- c(T.sim, T.star)
        coord.sim <- rbind(coord.sim, coord.star)
    }
    ## next step requires spatial interpolation according to the simulated spike
    ## Ui <-    
    tot <- sum(Ti)
}


Pl   <- Polygon(cbind(X$position_x, X$position_y))
ID   <- "[0,1]x[0,1]"
Pls  <- Polygons(list(Pl), ID=ID)
SPls <- SpatialPolygons(list(Pls))
df   <- data.frame(value=1, row.names=ID)
SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         # str(df)
trajectory <- SpatialPolygonsDataFrame(SPls, df) 


p1 <- ggplot() + geom_point(data=df2, size=2.5,
    colour=gray(1-df2$scaledX), aes(x,y)) + geom_path(data=pos.coords,
    size=.1, colour="blue", aes(x=x,y=y)) +
    geom_point(aes(x=coord.sim[,1], y=coord.sim[,2]), color="red") + coord_fixed() +
    ylim(0,100) + xlim(0,100) + theme_classic()
