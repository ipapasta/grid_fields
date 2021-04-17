## 14:33 - counter = 7
## 14:43 - counter = 11
## 14:54 - counter = 14
## 60mins- counter = 18
## 16:28 - counter = 38
if(FALSE){
    source(file="load_data.R")
    source(file="Functions.R")
    ## 
    load(file="current_intensity.RData")
    load(file="mle_oscillating.RData")
    load(file="mle_oscillating_HD.RData")
    load(file="/home/ipapasta/Desktop/mesh.RData")
    load(file="/home/ipapasta/Desktop/mle_new.RData")
    load(file="/home/ipapasta/Desktop/mle_oscillating_biased.RData")
    load(file="/home/ipapasta/Desktop/mle_oscillating_origin.RData")
    load(file="/home/ipapasta/Desktop/data/biased/Xbeta29.RData")
}
## load("mesh.RData")

## output1 needs order.HD=1
## output2 needs order.HD=5

options(warn=-1)
library(INLA)
library(gridExtra)
library(RColorBrewer)
library(scales)
library(Matrix)
## library(pals)
##
## 
##


## theta1 <- theta[1]
## theta2 <- theta[2]
## theta3 <- theta[3]
## rho     <- 19.60
## kappa   <- sqrt(8)/rho
## sigma   <- 0.13
## phi     <- -0.95
## sincpth <- sqrt(1-phi^2)/acos(phi)

## tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
## source("Functions.R")
## load("/home/ipapasta/Software/R/grid_fields/R/data/Simulated_data/simulated_cn_same_30var_v_grid5trials7_simulated/fitted_model.osc_HD.RData")
## load("/home/ipapasta/Software/R/grid_fields/R/data/Simulated_data/simulated_cn_same_30var_v_grid5trials7_simulated/mesh.RData")
## load("/data/ipapasta/biased/Xbeta0.RData")

counter <- 0
Xestimates<-pthetapc.prop.marg.post_osc_HD (par.theta=opt.theta$par, hyperpar=hyperpar, data=data, X=Xinit, beta=betainit, mesh=mesh, A=A, tA=tA, W=W, tW=tW, DW=DW,
                                Aobs=Aobs, tAobs=tAobs, order.HD=order.HD, type="biased", acc=1e-7, print.verbose=FALSE, return.X=TRUE)

p            <- mesh$n
## k            <- 12
alpha        <- 2

mycoords   <- SpatialPoints(cbind(Y$position_x, Y$position_y))
## -----------------------------------------------------------------
## Get Spatial effects and Temporal effects

## X <- Xhat[[length(Xhat)]]

X     <- Xest## (obj.to.save$Xbeta.optimizer)$mx.thetay
beta  <- betaest## (obj.to.save$Xbeta.optimizer)$mbeta.thetay
beta  <- fit.space.direction.time$summary.fixed$mean
Xest <- Matrix(fit.space.direction.time$summary.random[[1]]$mean, ncol=1)
X   <- Matrix(fit.space.direction.time$summary.random[[1]]$mean, ncol=1)
Xbeta <- rBind(t(t(X)), t(t(beta)))



## load("mesh.RData")


## Aobs: Matrix projector for spatial effects (observed)
## Bobs: Matrix projector for temporal effects (observed)
## B: design matrix for integral evaluation over time domain
## A    <- (obj.to.save$A.Aobs.L)[[1]]
## Aobs <- (obj.to.save$A.Aobs.L)[[2]]
## L    <- (obj.to.save$A.Aobs.L)[[3]]

n         <- nrow(data$Y)
one.n     <- matrix(rep(1, n), ncol=1)
## lambdafit <- exp((as.numeric(beta)*one.n) + Aobs%*% Xest)  ## %>% unlist %>% as.numeric %>% exp
## plot(trajectory)
## points(mycoords, pch=16, col=2, cex=2*(rank(lambdafit)/(1+n)))

## 
## predicted lambda
##
pos.coords <- data.frame(data$Ypos$coords)
names(pos.coords) <- c("x", "y")


slack <- 5
maxx  <- max(pos.coords[,"x"]+slack)
maxy  <- max(pos.coords[,"y"]+slack)
minx  <- min(pos.coords[,"x"]-slack)
miny  <- min(pos.coords[,"y"]-slack)

N      <- 1000
coords  <- expand.grid(seq(minx, maxx, len=N), seq(miny, maxy, len=N))

## coords  <- matrix(rep(c(33, 64), 100), byrow = TRUE, ncol=2)
one.vec <- matrix(rep(1, nrow(coords)), ncol=1)

## proj.s.pred  <- inla.mesh.projector(mesh, loc=matrix(coords, nrow=1))
Aosc.pred         <- inla.mesh.projector(mesh, loc=as.matrix(coords))$proj$A
Ahd.pred     <- inla.mesh.projector(mesh.hd, loc=rep(pi, dim(coords)[1]))$proj$A
## Ahd.pred       <- circular.make.A(seq(0.2, 2*pi-0.2, length = 100), order=order.HD)
## ---------------------------------------------------------------

Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
lambdapred  <- exp(matrix((as.numeric(beta)*one.vec), ncol=1) + Apred%*%Xest)
lambdapred  <- exp(Apred%*%Xest)
## grid(100, 100, lwd = 2)
## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 





library(pals)
X.scaled <- rank((Xbeta[1:mesh$n]))/(mesh$n)
X.scaled <- rank((X[1:mesh$n]))/(mesh$n)
df       <- data.frame(x=coords[,1], y=coords[,2], intensity=as.matrix(lambdapred))
df2      <- data.frame(x=mesh$loc[,1], y=mesh$loc[,2], scaledX=X.scaled)


## plot log intensity with base R image (use image.plot from fields for colour scale)
intens.matrix <- df %>% pivot_wider(names_from=c("x"), values_from=c("intensity")) %>% select(-c("y")) %>% unname
par(mfrow=c(1,2))
image.plot(seq(minx, maxx, len=N), seq(minx, maxx, len=N), t(intens.matrix), col=ocean.balance(100), zlim=c(exp(-4),exp(4)),
           xlab="XXX", ylab="YYY")
image.plot(seq(minx, maxx, len=N), seq(minx, maxx, len=N), t(log(intens.matrix)), col=ocean.balance(100), zlim=c(-4,4),
           xlab="XXX", ylab="YYY")

## logintensitytmp <- as(as.matrix(df), "dgTMatrix")
## logintensity.triplets <- cbind(cbind(logintensity@i+1, logintensity@j+1), logintensity@x)


pos.coords <- data.frame(data$Ypos$coords)
names(pos.coords) <- c("x", "y")

p1  <- ggplot()  +
    geom_point(data=df2, size=4.5, colour=gray(1-df2$scaledX), aes(x,y)) +
    geom_path(data=pos.coords, size=.1, colour="blue", aes(x=x,y=y)) +
    ## geom_path(data=, size=.1, colour="blue", aes(x=position_x,y=position_y)) +
    ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)), aes(x=position_x, y=position_y),color="red") + coord_fixed() +
    geom_point(data=data$Y, aes(x=position_x, y=position_y),color="red") + coord_fixed() +
    ylim(0,100) + xlim(0,100) + theme_classic()
p2  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=intensity), interpolate=TRUE) +
    ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)), aes(x=position_x, y=position_y),color="red") +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(df$intensity),max(df$intensity)))+
    ## scale_fill_viridis_c(option="Viridis")+
    ## scale_fill_gradient2(low = "darkblue", mid = "green", high = "red", midpoint = 0.4)  +
    coord_fixed()+ ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
                   ##            aes(x=position_x, y=position_y),color="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))  
p3  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=log(intensity)), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(log(df$intensity)),max(log(df$intensity))))+
    ## scale_fill_viridis_c(option="Viridis")+
    coord_fixed()+ ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
               ## aes(x=position_x, y=position_y),color="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))

p3  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=log(intensity)), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(-4,4))+
    ## scale_fill_viridis_c(option="Viridis")+
    coord_fixed()+ ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
               ## aes(x=position_x, y=position_y),color="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))






a <- grid.arrange(p1, p2, p3, nrow=1)



p3  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=log(intensity)), interpolate=TRUE)  +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(-4,4))+
    xlim(0,100)+
    ylim(0,100)+
    coord_fixed()+
        theme_classic()

p4  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=log(intensity)), interpolate=FALSE)  +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(-4,4))+
    xlim(0,100)+
    ylim(0,100)+
    coord_fixed()+
    theme_classic() 

plot(p3)

grid.arrange(p3,p4,nrow=1)

Hesstmp <- as(Hessianest, "dgTMatrix")
Hesstmp.triplets <- cbind(cbind(Hesstmp@i+1, Hesstmp@j+1), Hesstmp@x)

quilt.plot(Hesstmp.triplets[,1], Hesstmp.triplets[,2], Hesstmp.triplets[,3],
           col=ocean.balance(100), asp=1)

## ggsave(filename="/home/ipapasta/Desktop/preliminary_intensity_2.pdf", a)
## ggsave(filename="/home/ipapasta/Dropbox/org/Research/TeX/grid_fields/intensity.s", a)


## poster plot
p1  <- ggplot()  +
    geom_point(data=df2, size=5.5, colour=gray(1-df2$scaledX), aes(x,y)) +
    geom_path(data=pos.coords, size=.1, colour="blue", aes(x=x,y=y)) +
    ## geom_path(data=, size=.1, colour="blue", aes(x=position_x,y=position_y)) +
    ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)), aes(x=position_x, y=position_y),color="red") + coord_fixed() +
    geom_point(data=data$Y, aes(x=position_x, y=position_y),color="red") + coord_fixed() +
    ylim(5,100) + xlim(0,100) + theme_void()

library(ggpubr)
p1 <- p1  +  ggpubr::theme_transparent()


ggsave(filename="poster_rodent_statistics.pdf", p1)
ggsave(filename="poster_rodent_statistics.png", plot = p1,
       bg = "transparent")



plot(mesh, asp=1, main= "spatial effect - model fit")
points(mesh$loc, col=gray(1-X.scaled), pch=16, cex=2,asp=1)
points(mycoords, cex=.8, col=2, pch=16)
lines(as.matrix(pos.coords), lwd=.5, col="blue" )

plot(mesh, asp=1, main= "spatial effect - model fit")
points(mesh$loc, pch=16, cex=as.numeric(10*exp(Xest)/sum(exp(Xest))),asp=1)
points(mycoords, cex=.8, col=2, pch=16)
lines(as.matrix(pos.coords), lwd=.5, col="blue" )



## ------------------------
## effect of time
## ------------------------

Zest <- fit.space.direction.time$summary.random[[2]]$mean
time.grid <- seq(0.000, max(T.data), len=4000)
proj.t.pred  <- inla.mesh.projector(mesh1d, loc=as.matrix(time.grid))
Atemp.pred      <- proj.t.pred$proj$A
## Aosc.pred      <- proj.s.pred$proj$A
## Ahd.pred     <- inla.mesh.projector(mesh.hd, loc=theta.seq)$proj$A
one.vec <- matrix(rep(1, length(time.grid)), ncol=1)
## ---------------------------------------------------------------
## Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
lambdapred.temp  <- exp(matrix(Atemp.pred%*%Zest))
loglambdapred.temp  <- Atemp.pred%*%Zest

plot(time.grid, loglambdapred.temp, type="l")
plot(time.grid, lambdapred.temp, type="l")
abline(h=1)

## grid(00, 100, lwd = 2)
## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 


## ------------------------
## effect of head-direction
## ------------------------
## credible intervals
Sigma <- solve(Hessianest)

coords  <- matrix(rep(c(90,60), 50), byrow=T, ncol=2)
coords  <- matrix(rep(c(65,20), 50), byrow=T, ncol=2)
coords  <- matrix(rep(c(10,10), 50), byrow=T, ncol=2)
theta.seq <- seq(0., 2*pi, len=50)
proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
Aosc.pred      <- proj.s.pred$proj$A
Ahd.pred     <- inla.mesh.projector(mesh.hd, loc=theta.seq)$proj$A
one.vec <- matrix(rep(1, length(theta.seq)), ncol=1)

## ---------------------------------------------------------------
Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
lambdapred  <- exp(matrix(Apred%*%Xest))
loglambdapred  <- Apred%*%Xest

## varloglambdapred <- Apred %*% (Sigma[1:length(Xest), 1:length(Xest)] %*% t(Apred))

par(mfrow=c(1,2))
plot(theta.seq, as.numeric(loglambdapred), type="l")
plot(theta.seq, as.numeric(lambdapred), type="l")


df.hd <- data.frame(theta=theta.seq, log.hd.effect=as.numeric(loglambdapred))

                    ## log.hd.effect.low = as.numeric(loglambdapred) - 1.96 * sqrt(diag(varloglambdapred)),
                    ## log.hd.effect.upp = as.numeric(loglambdapred) + 1.96 * sqrt(diag(varloglambdapred))

## g = ggplot(df.hd,aes(x=factor(theta),y=log.hd.effect, group=1)) +
##     geom_polygon(fill=NA, color="black") +
##     geom_line() +
##     coord_polar(start=pi) +
##     coord_polar() + theme_minimal()

## plot(g)


g = ggplot(df.hd, aes(x=theta, y=log.hd.effect)) +
    ## geom_errorbar(aes(ymin=log.hd.effect.low, ymax=log.hd.effect.upp), width=.1) +
    geom_polygon(fill=NA, color="black") +
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    scale_y_continuous(breaks=seq(0.,3.5, by=0.15), limits=c(0., 3.5))+
    coord_polar(start = pi, direction=1)  +
    theme_minimal()

plot(g)










## ## visualisation
## x11()


par(mfrow=c(1,1))
## diagonal
plot(theta.seq, rep(1,length(theta.seq)), ylim=c(exp(-4),exp(4)), col="white", axes=FALSE, xlab="", ylab="")
for(i in 1:100){
    coords  <- matrix(rep(c(i,i), 1000), byrow=T, ncol=2)
    theta.seq <- seq(0.1, 2*pi, len=1000)
    ## coords  <- matrix(rep(c(33, 64), 100), byrow = TRUE, ncol=2)
    order.HD <- 5
    ## proj.s.pred  <- inla.mesh.projector(mesh, loc=matrix(coords, nrow=1))
    proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
    Aosc.pred      <- proj.s.pred$proj$A
    Ahd.pred     <- inla.mesh.projector(mesh.hd, loc=theta.seq)$proj$A
    ## Ahd.pred       <- circular.make.A(rep(pi, dim(coords)[1]), order=order.HD)
    one.vec <- matrix(rep(1, length(theta.seq)), ncol=1)
    ## Ahd.pred       <- circular.make.A(seq(0.2, 2*pi-0.2, length = 100), order=order.HD)
    ## ---------------------------------------------------------------
    Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
    lambdapred  <- exp(Apred%*%Xest)
    loglambdapred  <- Apred%*%Xest
    ## grid(00, 100, lwd = 2)
    ## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 
    lines(theta.seq, exp(as.numeric((loglambdapred))), col=gray(level=i/100))
}

dev.off()


par(mfrow=c(1,1))
## diagonal
plot(theta.seq, rep(0,length(theta.seq)), add=TRUE, ylim=c(-5*7.128448e-06,  15*1.714683e-06), col="white")

plot(theta.seq, rep(0,length(theta.seq)), add=TRUE, ylim=c(0.13555, 0.1359), col="white")
for(i in 1:100){
    coords  <- matrix(rep(c(i,i), 1000), byrow=T, ncol=2)
    theta.seq <- seq(0.1, 2*pi, len=1000)
    ## coords  <- matrix(rep(c(33, 64), 100), byrow = TRUE, ncol=2)
    order.HD <- 3
    ## proj.s.pred  <- inla.mesh.projector(mesh, loc=matrix(coords, nrow=1))
    proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
    Aosc.pred      <- proj.s.pred$proj$A
    Ahd.pred       <- circular.make.A(theta.seq, order=order.HD)
    ## Ahd.pred       <- circular.make.A(rep(pi, dim(coords)[1]), order=order.HD)
    one.vec <- matrix(rep(1, length(theta.seq)), ncol=1)
    ## Ahd.pred       <- circular.make.A(seq(0.2, 2*pi-0.2, length = 100), order=order.HD)
    ## ---------------------------------------------------------------
    Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
    lambdapred  <- exp(matrix((as.numeric(beta)*one.vec), ncol=1) + Apred%*%Xest)
    loglambdapred  <- Apred%*%Xest
    ## grid(00, 100, lwd = 2)
    ## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 
    lines(theta.seq, as.numeric(lambdapred), col=gray(i/1170))
}


dev.off()


## middle horizontal
plot(theta.seq, rep(10,length(theta.seq)), add=TRUE, ylim=c(-.1,.1))
for(i in 1:100){
    coords  <- matrix(rep(c(i,50), 1000), byrow=T, ncol=2)
    theta.seq <- seq(0.1, 2*pi, len=1000)
    ## coords  <- matrix(rep(c(33, 64), 100), byrow = TRUE, ncol=2)
    order.HD <- 7
    ## proj.s.pred  <- inla.mesh.projector(mesh, loc=matrix(coords, nrow=1))
    proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
    Aosc.pred      <- proj.s.pred$proj$A
    Ahd.pred       <- circular.make.A(theta.seq, order=order.HD)
    ## Ahd.pred       <- circular.make.A(rep(pi, dim(coords)[1]), order=order.HD)
    one.vec <- matrix(rep(1, length(theta.seq)), ncol=1)
    ## Ahd.pred       <- circular.make.A(seq(0.2, 2*pi-0.2, length = 100), order=order.HD)
    ## ---------------------------------------------------------------
    Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
    lambdapred  <- exp(matrix((as.numeric(beta)*one.vec), ncol=1) + Apred%*%Xest)
    loglambdapred  <- Apred%*%Xest
    ## grid(00, 100, lwd = 2)
    ## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 
    lines(theta.seq, as.numeric(loglambdapred), col=gray(i/1170))
}




## intensity
plot(theta.seq, rep(10,length(theta.seq)), add=TRUE, ylim=c(0.9,1.1))
for(i in 1:100){
    coords  <- matrix(rep(c(i,50), 1000), byrow=T, ncol=2)
    theta.seq <- seq(0.1, 2*pi, len=1000)
    ## coords  <- matrix(rep(c(33, 64), 100), byrow = TRUE, ncol=2)
    order.HD <- 7
    ## proj.s.pred  <- inla.mesh.projector(mesh, loc=matrix(coords, nrow=1))
    proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
    Aosc.pred      <- proj.s.pred$proj$A
    Ahd.pred       <- circular.make.A(theta.seq, order=order.HD)
    ## Ahd.pred       <- circular.make.A(rep(pi, dim(coords)[1]), order=order.HD)
    one.vec <- matrix(rep(1, length(theta.seq)), ncol=1)
    ## Ahd.pred       <- circular.make.A(seq(0.2, 2*pi-0.2, length = 100), order=order.HD)
    ## ---------------------------------------------------------------
    Apred    <- inla.row.kron(Ahd.pred, Aosc.pred)
    lambdapred  <- exp(matrix((as.numeric(beta)*one.vec), ncol=1) + Apred%*%Xest)
    loglambdapred  <- Apred%*%Xest
    ## grid(00, 100, lwd = 2)
    ## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 
    lines(theta.seq, as.numeric(lambdapred), col=gray(i/1170))
}



plot(theta.seq, exp(as.numeric(lambdapred)),)




if(FALSE)
    {
        ## -----------------------------------------------------
        ## Spatial effects plots
        ## -----------------------------------------------------
        pdf(file="~/Desktop/mfit_hotspots_blue.pdf")
        ## 
        X.scaled <- rank((Xbeta[1:nrow(mesh$loc)]))/
            (mesh$loc %>% nrow)
        ## mesh with spatial effects
        par(mar=0 %>% rep(4))
        plot(mesh, asp=1,
             main= "spatial effect - model fit")
        points(mesh$loc, col=gray(1-X.scaled),
               pch=16, cex=2,asp=1)
        points(coo, cex=.8, col=2, pch=16)
        ## HotSpots
        hotspots <- which(lambdafit>0.005)
        points(data[hotspots, c("Longitude", "Latitude")],
               col=4, cex=1.2, pch=16)

        
        points((new.data%>% select(Longitude,Latitude)) %>%
               tail(n=3), col=4, pch=16, cex=4)
        ## Linde accidents

        dev.off()

        plot(mesh, asp=1,main= "Linde data - red points")    
        points(mesh$loc,col=gray(1-X.scaled),
               pch=16, cex=1.2,asp=1)
        points(coo, cex=0.5, col=2, pch=16)
        ##

        plot(mesh, asp=1,
             main= "Locations to predict - blue points")    
        points(mesh$loc, col=gray(1-X.scaled),
               pch=16, cex=1.2,asp=1)
        points(new.data$Long,new.data$Lat,
               col=4, pch=16, cex=1)



        ## plot(density(lambdafit),main="")
        ## abline(v=lambdapred)
        
        ## plot(sort(lambda),sort(rank(lambda)/length(lambda)))
        ## abline(v=lambdapred)
        ## abline(h=.65)
    }
