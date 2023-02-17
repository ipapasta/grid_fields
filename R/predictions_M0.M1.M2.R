library(ggplot2)
library(gridExtra)
library(pals)
## ******************************************
## Model M0
## ******************************************

N       <- 100
coords  <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100)) %>% unname
predict.data.space  <- data.frame(coords.x1=coords[,1], coords.x2=coords[,2])
lambda.space        <- predict(fit.space, predict.data.space, ~ Intercept + spde2)
## lambda.space.exp    <- predict(fit.space, predict.data.space, ~ exp(Intercept + spde2))
lambda.space.M2.space.time        <- predict(fit.space.time, predict.data.space, ~ Intercept + spde2) 

p.space  <- ggplot(lambda.space, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space$mean),max(lambda.space$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))

p.space.exp  <- ggplot(lambda.space.exp, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.exp$mean),max(lambda.space.exp$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))

library(gridExtra)
grid.arrange(p.space, p.space.exp, nrow=1)

p.space.M2.space.time  <- ggplot(lambda.space.M2.space.time, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space$mean),max(lambda.space$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))

grid.arrange(p.space, p.space.M2.space.time, nrow=1)

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



## ------------------------------------------------------------
## Intensity across space averaged over head direction and
## Intensity across direction averaged over spatial coordinates
## ------------------------------------------------------------
N <- 100
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
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    geom_line(aes(x=hd, y=mean)) +
    ylim(-4,0)+
    ## geom_hline(yintercept=0, colour="grey")+
    coord_polar(start = pi, direction=1) + theme_minimal()


grid.arrange(p.space.direction.average.dir, p.space.direction.average.coord, nrow=1)




## comparison of fitted log intensities -  M0 vs M1 averaged across head direction 
grid.arrange(p.space, p.space.direction.average.dir, nrow=1)






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
lambda.M0.path             <- predict(fit.space, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Ypos$time)
lambda.M1.path             <- predict(fit.space.direction, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Ypos$time)
## lambda.M1.2.path           <- predict(fit.space.direction2, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Ypos$time)
lambda.M2.space.time.path  <- predict(fit.space.time, predict.intensity.on.path, ~ Intercept + spde2 + spde1) %>% mutate(time = Ypos$time)

## TEST SET
## for vertical ribbons use geom_rect
df.vert.ribbons <- dat$X %>% dplyr::filter(!(index.CV %in% train.index)) %>% group_by(index.CV) %>%
    summarize(minT = min(synced_time), maxT = max(synced_time)) %>%
    mutate(midT=(minT + maxT)/2,
           clump = setdiff(1:K, train.index))


## 
## cross validation predictions
## 
predict.intensity.on.path  <- data.frame(coords.x1=Yposraw$coords[,1], coords.x2=Yposraw$coords[,2], hd=Yposraw$hd, firing_times=Yposraw$time)
lambda.M0.path             <- predict(fit.space, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Yposraw$time)
lambda.M1.path             <- predict(fit.space.direction, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Yposraw$time)
## lambda.M1.2.path           <- predict(fit.space.direction2, predict.intensity.on.path, ~ Intercept + spde2) %>% mutate(time = Yposraw$time)
lambda.M2.space.time.path  <- predict(fit.space.time, predict.intensity.on.path, ~ Intercept + spde2 + spde1) %>% mutate(time = Yposraw$time)
plot(lambda.M1.path$mean, lambda.M1.2.path$mean, pch=16, cex=.2)
abline(a=0,b=1)
plot(lambda.M0.path$mean, lambda.M1.path$mean, pch=16, cex=0.2)
abline(a=0, b=1)
plot(lambda.M0.path$mean, lambda.M2.space.time.path$mean, pch=16, cex=.2)
abline(a=0,b=1)





i  <- 1
M  <- 50
T1 <- (i-1)*M + 1
T2 <- i*M
p.M2.space.time.path.time  <- ggplot() +
    geom_ribbon(data=lambda.M2.space.time.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="blue") +
    geom_ribbon(data=lambda.M0.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.1, colour="black") +
    ## geom_line(data=lambda.M2.space.time.path, aes(x=time,y=mean), colour="gray") +
    ## geom_line(data=lambda.M0.path, aes(x=time,y=mean), colour="black") +
    geom_point(data=dat$Y, aes(x=firing_times, rep(-7.8, length(firing_times))), colour="red") +
    ## geom_point(data=Y, aes(x=firing_times, rep(-7.8, length(firing_times))), colour="red") +
    geom_rect(data=df.vert.ribbons, aes(ymin=-10, ymax=10, xmin=minT, xmax=maxT), alpha=0.5)+
    annotate("text", x=df.vert.ribbons$midT, y=9, label= as.character(df.vert.ribbons$clump))+
    ## annotate("text", x=df.vert.ribbons$midT, y=7, label= as.character(round(pred.means.M2.space.time,1)))+
    ## annotate("text", x=df.vert.ribbons$midT, y=8, label= as.character(round(pred.means.M0,1)))+
    xlim(T1, T2) +
    ## ylim(-8,1.5)+
    theme_classic() + theme(legend.text=element_text(size=11))
p.M2.space.time.path.time


## including predictions from space.direction model
i  <- 2
M <- 100
T1 <- (i-1)*M + 1
T2 <- i*M

p.M2.space.time.path.time  <- ggplot() +
    geom_rect(data=df.vert.ribbons, aes(ymin=-10, ymax=10, xmin=minT, xmax=maxT), alpha=0.1, size=1, linetype=2)+
    geom_ribbon(data=lambda.M2.space.time.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="brown") +
    geom_ribbon(data=lambda.M0.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.2, colour="black") +
    geom_ribbon(data=lambda.M1.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.1, colour="khaki3") +
    ## geom_line(data=lambda.M2.space.time.path, aes(x=time,y=mean), colour="gray") +
    ## geom_line(data=lambda.M0.path, aes(x=time,y=mean), colour="black") +
    geom_point(data=dat$Y, aes(x=firing_times, rep(-7.8, length(firing_times))), colour="red") +
    ## geom_point(data=Y, aes(x=firing_times, rep(-7.8, length(firing_times))), colour="red") +
    ## annotate("text", x=df.vert.ribbons$midT, y=9, label= paste0("test set id: ", as.character(df.vert.ribbons$clump)))+
    annotate("text", x=df.vert.ribbons$midT, y=5, label= paste0("O: ",as.character(round(out.list$obs.firings))))+
    annotate("text", x=df.vert.ribbons$midT, y=8, label= paste0("M2: ", as.character(round(out.list$pred.means.M2.space.time))))+
    annotate("text", x=df.vert.ribbons$midT, y=7, label= paste0("M1: ", as.character(round(out.list$pred.means.M1.space.direction))))+
    annotate("text", x=df.vert.ribbons$midT, y=6, label= paste0("M0: ",as.character(round(out.list$pred.means.M0))))+
    xlim(T1, T2) +
    ylab("log intensity")+
    ## ylim(-8,1.5)+
    theme_classic() + theme(legend.text=element_text(size=11))

## hide test
p.M2.space.time.path.time.hide  <- ggplot() +
    geom_ribbon(data=lambda.M2.space.time.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="brown") +
    geom_ribbon(data=lambda.M0.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.2, colour="black") +
    geom_ribbon(data=lambda.M1.path, aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.1, colour="khaki3") +
    geom_point(data=dat$Y, aes(x=firing_times, rep(-7.8, length(firing_times))), colour="red") +
    annotate("text", x=df.vert.ribbons$midT, y=5, label= paste0("O: ",as.character(round(out.list$obs.firings))))+
    annotate("text", x=df.vert.ribbons$midT, y=8, label= paste0("M2: ", as.character(round(out.list$pred.means.M2.space.time))))+
    annotate("text", x=df.vert.ribbons$midT, y=7, label= paste0("M1: ", as.character(round(out.list$pred.means.M1.space.direction))))+
    annotate("text", x=df.vert.ribbons$midT, y=6, label= paste0("M0: ",as.character(round(out.list$pred.means.M0))))+
    geom_rect(data=df.vert.ribbons, aes(ymin=-10, ymax=10, xmin=minT, xmax=maxT), alpha=0.1, size=1, linetype=2)+
    xlim(T1, T2) +
    ylab("log intensity")+
    ## ylim(-8,1.5)+
    theme_classic() + theme(legend.text=element_text(size=11))

p.covar.hd <- ggplot() + geom_point(data=dat$X, aes(x=synced_time, y=hd)) +
    xlim(T1, T2) +
    xlab("time") +
    ylab("Direction of navigation")+
    theme_classic() + theme(legend.text=element_text(size=11))

a <- grid.arrange(p.M2.space.time.path.time, p.covar.hd, nrow=2)
ggsave(plot=a, filename="CV_preds_2_minutes_swap0.pdf", height=10)
ggsave(plot=p.covar.hd, filename="direction_of_navigation.pdf", height=10)
ggsave(plot=p.M2.space.time.path.time,
       filename="/Users/ipapasta/Documents/Research/Talks/Oscillating_fields_talk/CV_preds_10_minutes_swap0.pdf", width=13, height=7)
ggsave(plot=p.covar.hd,
       filename="/Users/ipapasta/Documents/Research/Talks/Oscillating_fields_talk/direction_of_navigation.pdf", height=10)


p.M0.path.time  <- ggplot(lambda.M0.path) +
    geom_ribbon(aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") +
    geom_line(aes(x=time,y=mean)) +
    xlim(T1, T2) +
    ylim(-8,1.5)+
    geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025, length(firing_times)))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))
##
##
## grid.arrange(p.M0.path.time, p.M2.space.time.path.time, ncol=1)
## 
p.M1.path.time  <- ggplot(lambda.M1.path) +
    geom_ribbon(aes(x= time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") + 
    geom_line(aes(x=time,y=mean)) +
    xlim(T1, T2) +
    ylim(-8,1.5)+
    geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025), length(firing_times))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))
##
p.M1.2.path.time  <- ggplot(lambda.M1.2.path) +
    geom_ribbon(aes(x=time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") + 
    geom_line(aes(x=time,y=mean)) +
    xlim(T1, T2) +
    ylim(-8,1.5)+
    geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025), length(firing_times))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))
## 
p.hd  <- ggplot(lambda.M1.path, aes(x=time,y=hd)) + geom_line() +
    xlim(T1, T2) +
    theme_classic() + theme(legend.text=element_text(size=11))



grid.arrange(p.M0.path.time, p.M1.path.time, p.M1.2.path.time, p.hd, ncol=1)

save(predict.data.space.direction.full, predict.intensity.on.path, lambda.M0.path, lambda.M1.path, lambda.space.direction.full, file="lambda_M0_M1.RData")

load("lambda_M0_M1.RData")


## ---------------------------------------
## Predictive distribution and CRPS scores
## ---------------------------------------
samp.M0.space_2           <- generate(object = fit.space_2, n.samples = 5000,num.threads=ncores)
samp.M2.space.time_2      <- generate(object = fit.space.time_2,  n.samples = 5000, num.threads=ncores)
samp.M1.space.direction_2 <- generate(object = fit.space.direction_2,  n.samples = 5000, num.threads=ncores)

samp.M0.space_5           <- generate(object = fit.space_5, n.samples = 5000,num.threads=ncores)
samp.M2.space.time_5      <- generate(object = fit.space.time_5,  n.samples = 5000, num.threads=ncores)
samp.M1.space.direction_5 <- generate(object = fit.space.direction_5,  n.samples = 5000, num.threads=ncores)

samp.M0.space_10           <- generate(object = fit.space_10, n.samples = 5000,num.threads=ncores)
samp.M2.space.time_10      <- generate(object = fit.space.time_10,  n.samples = 5000, num.threads=ncores)
samp.M1.space.direction_10 <- generate(object = fit.space.direction_10,  n.samples = 5000, num.threads=ncores)

samp.M0.space_20           <- generate(object = fit.space_20, n.samples = 5000,num.threads=ncores)
samp.M2.space.time_20      <- generate(object = fit.space.time_20,  n.samples = 5000, num.threads=ncores)
samp.M1.space.direction_20 <- generate(object = fit.space.direction_20,  n.samples = 5000, num.threads=ncores)
