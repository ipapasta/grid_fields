library(ggplot2)
library(gridExtra)
## ******************************************
## Model M0
## ******************************************
N       <- 100
coords  <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100)) %>% unname
predict.data.space  <- data.frame(coords.x1=coords[,1], coords.x2=coords[,2])
lambda.space        <- predict(fit.space, predict.data.space, ~ Intercept + spde2)

p.space  <- ggplot(lambda.space, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space$mean),max(lambda.space$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))

## ******************************************
## Model M1
## ******************************************
## ------------------------------------------------------------
## log intensity across space for fixed direction and
## log intensity across direction for fixed coordinates
## ------------------------------------------------------------

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
lambda.space.direction.dir.fixed       <- predict(fit.space.direction, predict.data.space.direction.dir.fixed, ~ Intercept + spde2)
lambda.space.direction.coords.fixed    <- predict(fit.space.direction, predict.data.space.direction.coords.fixed, ~ Intercept + spde2)


p.space.direction.dir.fixed  <- ggplot(lambda.space.direction.dir.fixed, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.direction.dir.fixed$mean),max(lambda.space.direction.dir.fixed$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))

p.space.direction.coord.fixed <- ggplot(lambda.space.direction.coords.fixed) + 
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

grid.arrange(p.space.direction.dir.fixed, p.space.direction.coord.fixed, nrow=1)



## animation

## ------------------------------------------------------------
## Intensity across space averaged over head direction and
## Intensity across direction averaged over spatial coordinates
## ------------------------------------------------------------
coords.dir   <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100), seq(0,2*pi,len=50)) %>% unname
predict.data.space.direction.full  <- data.frame(hd=coords.dir[,3], coords.x1=coords.dir[,1], coords.x2=coords.dir[,2])
lambda.space.direction.full        <- predict(fit.space.direction, predict.data.space.direction.full, ~ Intercept + spde2)

lambda.space.direction.average.dir <- lambda.space.direction.full %>% mutate(hd=coords.dir[,3], coords.x1=coords.dir[,1], coords.x2=coords.dir[,2]) %>%
    group_by(coords.x1, coords.x2) %>%
    summarize(mean = mean(mean))

lambda.space.direction.average.coord <- lambda.space.direction.full %>% mutate(hd=coords.dir[,3], coords.x1=coords.dir[,1], coords.x2=coords.dir[,2]) %>%
    group_by(hd) %>%
    summarize(mean = mean(mean))

p.space.direction.average.dir  <- ggplot(lambda.space.direction.average.dir, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.direction.average.dir$mean),max(lambda.space.direction.average.dir$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))


p.space.direction.average.coord <- ggplot(lambda.space.direction.average.coord) + 
    ## geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    ## scale_y_continuous(breaks=seq(min(lambda.space.direction.average.coord$q0.025) %>% floor,
    ##                               max(lambda.space.direction.average.coord$q0.975) %>% ceiling, by=0.5),
    ##                    limits=c(min(lambda.space.direction.average.coord$q0.025) %>% floor,
    ##                             max(lambda.space.direction.average.coord$q0.975) %>% ceiling))+
    geom_line(aes(x=hd, y=mean)) +
    geom_hline(yintercept=0, colour="grey")+
    coord_polar(start = pi, direction=1) + theme_minimal()


grid.arrange(p.space.direction.average.dir, p.space.direction.average.coord, nrow=1)


## comparison of fitted log intensities -  M0 vs M1 averaged across head direction 
grid.arrange(p.space, p.space.direction.average.dir, nrow=1)



## 
## code for animation
## 
df.animation <- lambda.space.direction.full %>%
    mutate(hd=as.factor(hd))
p.space.direction  <- ggplot(df.animation, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.direction.full$mean),max(lambda.space.direction.full$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))


p.space.direction + transition_states(hd, transition_length = 2, state_length = 1)




## ******************************************
## Model M2
## ******************************************
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


## ******************************************
## Point process intensity plotted over time (path)
## ******************************************
## extract times, coordinates and directions along the path

predict.intensity.on.path  <- data.frame(coords.x1=Ypos$coords[,1], coords.x2=Ypos$coords[,2], hd=Ypos$hd, firing_times=Ypos$time)
lambda.M0.path        <- predict(fit.space, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Ypos$time)
lambda.M1.path        <- predict(fit.space.direction, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Ypos$time)
lambda.M2.path        <- predict(fit.space.direction.time, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Ypos$time)

T1 <- 450
T2 <- 500

p.M0.path.time  <- ggplot(lambda.M0.path) +
    geom_ribbon(aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") +
    geom_line(aes(x=time,y=mean)) +
    xlim(T1, T2) +
    geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025, length(firing_times)))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))

p.M1.path.time  <- ggplot(lambda.M1.path, ) +
    geom_ribbon(aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") + 
    geom_line(aes(x=time,y=mean)) +
    xlim(T1, T2) +
        geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025, length(firing_times)))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))

p.hd  <- ggplot(lambda.M1.path, aes(x=time,y=hd)) + geom_line() +
    xlim(T1, T2) +
    theme_classic() + theme(legend.text=element_text(size=11))

grid.arrange(p.M0.path.time, p.M1.path.time, p.hd)

plot(lambda.M0.path$mean, lambda.M1.path$mean, pch=16, cex=0.1)

abline(a=0, b=1)


