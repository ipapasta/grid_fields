## ---------
## Model M0
## ---------

## plot estimated intensity
pxl         <- pixels(mesh, nx=500, ny=500)    
Aosc.test   <- inla.mesh.projector(mesh, loc=pxl)$proj$A
quilt.plot(pxl@coords[,1], pxl@coords[,2], as.numeric(Aosc.test %*% fit.space$summary.random$spde2$mean), ylim=c(0,100),xlim=c(0,100))
lambda.osc         <- predict(fit.space, data=NULL, formula = ~ spde2_eval(cbind(pxl@coords[,1], pxl@coords[,2])) )
lambda.osc.spatial <- predict(fit.oscillating.rgeneric, pxl, ~ spde2)

p1  <- ggplot(lambda.osc %>% mutate(coords.x1=pxl@coords[,1],coords.x2=pxl@coords[,2]),
              aes(x=coords.x1,y=coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
    geom_point(data=data$Y, size=0.5, aes(x=position_x, y=position_y),color="red") +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar") +
    xlim(0,100) + ylim(0,100) +
    coord_fixed() +
    theme_classic()

p2 <- ggplot() + gg(lambda.osc.spatial) +
    gg(mycoords, color="red", size=0.2) + #mycoords is object containing the firing events. This object is defined in load_data.R
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar") +
    xlim(0,100) + ylim(0,100) + 
    coord_equal() + theme_classic()

grid.arrange(p1,p2,nrow=1)


p3 <- ggplot() + gg(lambda.hd.fixed) +
    gg(mycoords, color="red", size=0.2) + #mycoords is object containing the firing events. This object is defined in load_data.R
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar") +
    xlim(0,100) + ylim(0,100) + 
    coord_equal() + theme_classic()


## ---------
## Model M1
## ---------
N            <- 100
coords       <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100)) %>% unname
dir.fixed    <- rep(pi/3, length=100*100)
coords.fixed <- matrix(rep(c(40,45),N*N), ncol=2,byrow=TRUE)
dir          <- seq(0,2*pi, length=N*N)
predict.data <- data.frame(coords.x1 = coords[,1], coords.x2 = coords[,2], hd = dir.fixed)

## 
## lambda.hd.fixed <- data.frame(firing_times=times, hd=dir.fixed, coords.x1=coords[,1], coords.x2=coords[,2])
## pr.int.full <- predict(fit.space.direction, dat.2.predict, ~ Intercept + spde2 + time)
lambda.hd.fixed     <- predict(fit.space.direction, data=NULL, formula = ~ spde2_eval(cbind(coords[,1],coords[,1]), rep(pi, N*N)) ) 
lambda.hd.fixed     <- predict(fit.space.direction, data=NULL, formula = ~ spde2_eval(cbind(mesh$loc[,1],mesh$loc[,2]), rep(pi, mesh$n)) ) 
lambda.hd.fixed     <- predict(fit.space.direction, data=NULL, formula = ~ spde2_eval(list(spatial=cbind(coords[,1],coords[,2]), direction=dir.fixed)) )

lambda.hd.fixed     <- predict(fit.space.direction, data=NULL, formula = ~spde2_eval(spatial=cbind(predict.data$coords.x1, predict.data$coords.x2), direction=predict.data$hd))

## lambda.coord.fixed  <- predict(fit.space.direction, NULL, ~ spde2_eval(coords.fixed, dir))

Aosc.test   <- inla.mesh.projector(mesh, loc=SpatialPoints(as.data.frame(coords)))$proj$A
Ahd.test    <- inla.mesh.projector(mesh.hd, loc=rep(pi, N*N))$proj$A
A.test      <- inla.row.kron(Ahd.test, Aosc.test)
quilt.plot(coords[,1], coords[,2], as.numeric(A.test %*% fit.space.direction$summary.random$spde2$mean))
points(mycoords,cex=0.3,pch=16,col=2)

quilt.plot(coords[,1], coords[,2], lambda.hd.fixed$mean)

lambda.hd.fixed$coords.x1 <- coords[,1]
lambda.hd.fixed$coords.x2 <- coords[,2]
quilt.plot(lambda.hd.fixed$coords.x1, lambda.hd.fixed$coords.x2, lambda.hd.fixed$mean, asp=1)
points(mycoords, cex=0.5, col="red")

p1  <- ggplot(lambda.hd.fixed, aes(x=coords.x1,y=coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
    ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)), aes(x=position_x, y=position_y),color="red") +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.hd.fixed$mean),max(lambda.hd.fixed$mean)))+
    coord_fixed()


p1 <- ggplot(lambda.coord.fixed) + 
    geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    scale_y_continuous(breaks=seq(-3.5,3.5, by=0.5), limits=c(-3.5, 3.5))+
    geom_line(aes(x=hd, y=mean)) +
    geom_hline(yintercept=0, colour="grey")+
    coord_polar(start = pi, direction=1) + theme_minimal()

## ---------
## Model M2
## ---------
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
