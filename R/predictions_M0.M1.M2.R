## ---------
## Model M0
## ---------
N       <- 100
coords  <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100)) %>% unname
predict.data.space  <- data.frame(coords.x1=coords[,1], coords.x2=coords[,2])
lambda.space        <- predict(fit.space.direction, predict.data.dir.fixed, ~ spde2)

p.space  <- ggplot(lambda.space, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(pr.int.full$mean),max(pr.int.full$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))




## ---------
## model M1
## ---------
N            <- 100
fixed.dir    <- pi
fixed.coord  <- c(40,45)
## 
coords       <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100)) %>% unname
coords.fixed <- matrix(rep(fixed.coord,N*N), ncol=2, byrow=TRUE)
dir          <- seq(0, 2*pi, length=N*N)
dir.fixed    <- rep(fixed.dir, length=100*100)

predict.data.space.direction.dir.fixed     <- data.frame(hd=dir.fixed, coords.x1=coords[,1], coords.x2=coords[,2])
predict.data.space.direction.coords.fixed  <- data.frame(hd=dir, coords.x1=coords.fixed[,1], coords.x2=coords.fixed[,2])
lambda.space.direction.dir.fixed       <- predict(fit.space.direction, predict.data.space.direction.dir.fixed, ~ spde2)
lambda.space.direction.coords.fixed    <- predict(fit.space.direction, predict.data.space.direction.coords.fixed, ~ spde2)


p.dir.fixed  <- ggplot(lambda.space.direction.dir.fixed, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(pr.int.full$mean),max(pr.int.full$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))
p.coord.fixed <- ggplot(lambda.space.direction.coords.fixed) + 
    geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    scale_y_continuous(breaks=seq(min(lambda.space.direction.coords.fixed$q0.025) %>% floor,
                                  max(lambda.space.direction.coords.fixed$q0.975) %>% ceiling, by=0.5),
                       limits=c(min(lambda.space.direction.coords.fixed$q0.025) %>% floor,
                                max(lambda.space.direction.coords.fixed$q0.975) %>% ceiling))+
    geom_line(aes(x=hd, y=mean)) +
    geom_hline(yintercept=0, colour="grey")+
    coord_polar(start = pi, direction=1) + theme_minimal()

grid.arrange(p.dir.fixed, p.coord.fixed, nrow=1)

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
## pr.int.full <- predict(fit.space.direction.time, NULL, ~ spde2_eval(coords, dir) + time_eval(times))
## pr.int.full <- predict(fit.space.direction.time, NULL, ~ spde2_eval(coords, dir) + time_eval(times))



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

## grid.arrange(p1,p2, nrow=1)

## pr.int.full <- predict(fit.space.direction.time, NULL, ~ c(Intercept = Intercept_latent))
## pr.int <- predict(fit.space.direction.time, NULL, ~ spde2_eval(coords, dir) + time_eval(times))
## pr.int <- predict(fit.space.direction.time,NULL, ~ Intercept + spde2(coords))
## p3  <- ggplot(pr.int.full, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
##     scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
##                          limits=c(min(pr.int.full$mean),max(pr.int.full$mean)))+
##     coord_fixed()+ 
##     theme_classic() + theme(legend.text=element_text(size=11))
## ## 


