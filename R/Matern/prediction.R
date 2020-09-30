## 14:33 - counter = 7
## 14:43 - counter = 11
## 14:54 - counter = 14
## 60mins- counter = 18
## 16:28 - counter = 38
source(file="load_data.R")
source(file="Functions.R")
## load(file="lgcp_spde_fit_shape_2_mesh_edge_size12_initial_values_0.28-0.45.RData")
## load(file="test.RData")
load(file="current_intensity.RData")
load(file="test.RData")
load("mesh.RData")
options(warn=-1)
library(gridExtra)
library(RColorBrewer)
library(scales)
library(Matrix)
##
## 
##
p            <- mesh$n
## k            <- 12
alpha        <- 2


## -----------------------------------------------------------------
## Get Spatial effects and Temporal effects
X     <- Xest## (obj.to.save$Xbeta.optimizer)$mx.thetay
beta  <- betaest## (obj.to.save$Xbeta.optimizer)$mbeta.thetay
Xbeta <- rBind(t(t(X)),t(t(beta)))      #transposes needed?


## Aobs: Matrix projector for spatial effects (observed)
## Bobs: Matrix projector for temporal effects (observed)
## B: design matrix for integral evaluation over time domain
## A    <- (obj.to.save$A.Aobs.L)[[1]]
## Aobs <- (obj.to.save$A.Aobs.L)[[2]]
## L    <- (obj.to.save$A.Aobs.L)[[3]]

n         <- nrow(data$Y)
one.n     <- matrix(rep(1, n), ncol=1)
lambdafit <- exp((as.numeric(beta)*one.n) + Aobs%*% Xest)  ## %>% unlist %>% as.numeric %>% exp
## plot(trajectory)
## points(mycoords, pch=16, col=2, cex=2*(rank(lambdafit)/(1+n)))

## 
## predicted lambda
##
pos.coords <- data.frame(do.call("rbind", data$Ypos$coords))
names(pos.coords) <- c("x", "y")

maxx <- max(pos.coords[,"x"])
maxy <- max(pos.coords[,"y"])
minx <- min(pos.coords[,"x"])
miny <- min(pos.coords[,"y"])

N      <- 5000
coords <- expand.grid(seq(minx, maxx, seq=N), seq(miny, maxy, seq=N))
one.vec <- matrix(rep(1, nrow(coords)), ncol=1)
Apred       <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coords))
lambdapred  <- exp(matrix((as.numeric(beta)*one.vec), ncol=1) + Apred%*%Xest)
## grid(100, 100, lwd = 2)
## %>% unlist %>% as.numeric %>% exp## matrix(Apred %*% X + Bpred %*% beta) %>% unlist %>% as.numeric %>% exp 



X.scaled <- rank((Xbeta[1:mesh$n]))/(mesh$n)
df       <- data.frame(x=coords[,1], y=coords[,2], intensity=as.matrix(lambdapred))
df2      <- data.frame(x=mesh$loc[,1], y=mesh$loc[,2], scaledX=X.scaled)


p1  <- ggplot() + geom_path(data=pos.coords, size=.1, colour="blue", aes(x=x,y=y)) +
    gg(mesh) + geom_point(data=df2, size=3, colour=gray(1-df2$scaledX), aes(x,y)) +
    geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
               aes(x=position_x, y=position_y),color="red") + coord_fixed() +theme_classic()
p2  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=intensity), interpolate=TRUE) +
    scale_fill_gradient2(low = "darkblue", mid = "green", high = "red", midpoint = 0.4)  +
    coord_fixed()+ ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
                   ##            aes(x=position_x, y=position_y),color="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))  
p3  <- ggplot(df, aes(x,y)) + geom_raster(aes(fill=log(intensity)), interpolate=TRUE) +
    scale_fill_gradient2(low = "darkblue", mid = "green", high = "red", midpoint = -1.5) +
    coord_fixed()+ ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
               ## aes(x=position_x, y=position_y),color="red") + 
    theme_classic() + theme(legend.text=element_text(size=11))  


a <- grid.arrange(p1, p2, p3, nrow=1)

ggsave(filename="/home/ipapasta/Dropbox/org/Research/TeX/grid_fields/intensity.ps", a)

plot(mesh, asp=1, main= "spatial effect - model fit")
points(mesh$loc, col=gray(1-X.scaled), pch=16, cex=2,asp=1)
points(mycoords, cex=.8, col=2, pch=16)
lines(as.matrix(pos.coords), lwd=.5, col="blue" )

plot(mesh, asp=1, main= "spatial effect - model fit")
points(mesh$loc, pch=16, cex=as.numeric(10*exp(Xest)/sum(exp(Xest))),asp=1)
points(mycoords, cex=.8, col=2, pch=16)
lines(as.matrix(pos.coords), lwd=.5, col="blue" )

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
