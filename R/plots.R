source("load_data.R")
## 
## X matrix (head directions?)
## 
X %>% dplyr::select(position_x_pixels, position_y_pixels) %>% plot(type="l", lwd=3)
## for(i in 1:nrow(X))
## {
##     points((X %>% dplyr::select(position_x_pixels, position_y_pixels))[i,], cex=1.5, col=2, pch=19)
## }
X %>% dplyr::select(position_x_pixels, position_y_pixels) %>% plot(type="l", lwd=3)
X %>% dplyr::select(position_x, position_y) %>% points(pch=16) 
X %>% dplyr::select(position_x, position_y) %>% plot


## 
## Y matrix (neuron firings?)
## 

pdf(height=7, width=11,file="/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/diff_synced_time.pdf")
plot(diff(X$synced_time)[1:100], type="l", ylab="time to next measurement (whole session) (diff of synced time)")
dev.off()

names(Y)

X %>% dplyr::select(position_x_pixels, position_y_pixels) %>% plot(type="l", lwd=3, col="white")

df1 <- X %>% dplyr::select(position_x_pixels, position_y_pixels, synced_time) %>% mutate(fire=0)
df2 <- Y %>% dplyr::select(position_x_pixels, position_y_pixels, firing_times)  %>% mutate(fire=1)

names(df1) <- c("x", "y", "time", "fire")
names(df2) <- c("x", "y", "time", "fire")

Z <- bind_rows(df1,df2) %>% arrange(time)



plot((Z %>% dplyr::select(x, y)), col = "white")
for(i in 1:nrow(Z))
{
    Sys.sleep(diff(Z$time)[i]/(200*mean(diff(Z$time))))
    points((Z %>% dplyr::select(x, y))[i,], col=ifelse(Z$fire[i]==0, "gray", "red"), pch=19, cex=Z$fire[i]+1)
}

plot(diff(Y$firing_times), type="l")



dft <- diff(Y$firing_times)
lambdahat <- mean(dft)

plot(Y$firing_times[1:200], rep(0, 200), pch=4)


dftnorm <- dft/lambdahat

plot(qexp((1:length(dft))/(1+length(dft))), sort(dftnorm))
abline(a=0, b=1)

plot(qexp((1:length(dft))/(1+length(dft)), lambdahat), sort(dft))
abline(a=0,b=1)





##
## X and Y
##
X %>% dplyr::select(position_x_pixels, position_y_pixels) %>% plot(type="l", lwd=3)
Y %>% dplyr::select(position_x_pixels, position_y_pixels) %>% points(pch=16, col=gray(rank(Y$hd)/(1+length(Y$hd))), cex=2)






## 
##
##


data("gorillas", package = "inlabru")
gonests <- gorillas$nests
goboundary <- gorillas$boundary
mycoords   <- SpatialPoints(cbind(Y$position_x_pixels, Y$position_y_pixels))
mycoordsDF <- SpatialPointsDataFrame(coords=mycoords, data=data.frame(hd=Y$hd))

Pl <- Polygon(cbind(X$position_x_pixels, X$position_y_pixels))
## tmp <- SpatialPolygons(Srl=list(tmp))
## SpatialPolygonsDataFrame(list(tmp))

ID   <- "[0,1]x[0,1]"
Pls  <- Polygons(list(Pl), ID=ID)
# str(Pls)
SPls <- SpatialPolygons(list(Pls))
# str(SPls)
df   <- data.frame(value=1, row.names=ID)
# str(df)
SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         # str(df)
trajectory <- SpatialPolygonsDataFrame(SPls, df) 



pdf(height=8, width=11,
    file="/home/ipapasta/Dropbox/org/Research/TeX/grid_fields/mouse_trajectory_triangulation.pdf")
par(oma=c(0,0,0,0))
par(mar=c(1,1,1,1))
par(mfrow=c(1,2))
## plot(trajectory)
## axis(1)
## axis(2)
## box()
plot(trajectory, asp=1)
points(mycoords, pch=16, col=2, cex=0.7)
plot(mesh, asp=1, main="")
points(mycoords, col=2, pch=16, cex=0.7)
dev.off()


postscript(height=600, width=800,
    file="/home/ipapasta/Dropbox/org/Research/TeX/grid_fields/mouse_trajectory_triangulation.ps")
par(oma=c(0,0,0,0))
par(mar=c(1,1,1,1))
par(mfrow=c(1,2))
## plot(trajectory)
## axis(1)
## axis(2)
## box()
plot(trajectory, asp=1)
points(mycoords, pch=16, col=2, cex=0.7)
plot(mesh, asp=1, main="")
points(mycoords, col=2, pch=16, cex=0.7)
dev.off()



mesh <- inla.mesh.2d(mycoords, max.edge=20, cutoff=5)

plot(mesh, asp=1, lwd=4)
## points(mycoords)
plot(trajectory, add=TRUE)

plot(mesh, asp=1)
points(mycoords, col=2, pch=16)
## plot(trajectory, add=TRUE)

ggplot() + gg(mesh) + gg(mycoords) + coord_fixed() + ggtitle("Firing events")


## specify correlation structure
boundary.samplers <- SpatialPointsDataFrame(coords=cbind(c(0, 0, 1, 1), c(0, 1, 0, 1)), data=data.frame(cbind(c(0, 0, 1, 1), c(0, 1, 0, 1))))
## SpatialPolygons(list(Polygons(list(Polygon(cbind(c(0, 0, 1, 1), c(0, 1, 0, 1)))), ID="[0,1]x[0,1]")))
mycoordsdf        <- SpatialPointsDataFrame(coords=cbind(Y$position_x_pixels, Y$position_y_pixels), data=Y)




matern <- inla.spde2.pcmatern(mesh, 
                              prior.sigma = c(.5, 0.01), 
                              prior.range = c(0.5, 0.01))

cmp    <- coordinates ~ mySmooth(map = coordinates,
                                 model = matern) +
    Intercept

fit    <- lgcp(cmp, mycoordsdf, samplers=boundary.samplers)

lambda <- predict(fit, pixels(mesh), ~ exp(mySmooth + Intercept))

loglambda <- predict(fit, pixels(mesh), ~ mySmooth + Intercept)

colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}

pl1 <- ggplot() + 
    gg(lambda) + 
    gg(boundary.samplers) +
    gg(mycoords) +
    ggtitle("LGCP fit to Points", subtitle = "(Response Scale)") + 
    coord_fixed() +
    colsc(lambda$median)

pl2 <- ggplot() + 
    gg(loglambda) + 
    gg(boundary.samplers) +
    ## gg(mycoords) +
    ggtitle("LGCP fit to Points", subtitle = "(Linear Predictor Scale)") + 
    coord_fixed() +
    colsc(loglambda$median)

multiplot(pl1, pl2, cols = 2)



## ----------------------------------
