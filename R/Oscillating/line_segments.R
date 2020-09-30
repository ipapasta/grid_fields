library(Rcpp)
library(tidyverse)
library(purrr)
library(INLA)   #install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(inlabru)
library(sp)
library(pryr)
library(fields)
library(nloptr)
## 
source("load_data.R")
source("Functions.R")
source("osc_precision.R")
source("priorXbeta_osc.R")
source("priortheta_osc.R")              #think about choice of prior distribution.
source("gradient_osc.R")
source("hessian_osc.R")
source("llik.R")
source("pthetamargipost_osc_stable.R")
sourceCpp("hessian_functions.cpp")

##
## create/load mesh object:
## 
load("mesh.RData")                     # ! rgdal missing from maths
                                        # computing server so
                                        # currently saved locally and
                                        # uploaded for practical
                                        # reasons. Request to be
                                        # installed.
## p <- mesh$n
options(warn=-1)                        #suppress warnings



## convert a line specified by co-ordinates x and y into a SpatialPolygon object
SPL.xy <- function(x, y){
    Pl  <- Polygon(rbind(x, y, x))
    Pls <-  Polygons(list(Pl), ID = "a")
    SPls <- SpatialPolygons(list(Pls))
    ## plot(SPls)
    return(SPls)
}

## create a SpatialPolygons object from specified by co-ordinates x and y into a SpatialPolygon object
SPL.xyz <- function(x, y, z, ID){
    Pl  <- Polygon(rbind(x, y, z))
    Pls <-  Polygons(list(Pl), ID = ID)
    ## SPls <- SpatialPolygons(list(Pls))
    ## plot(SPls)
    return(Pls)
}

if(FALSE){
    k    <- 3
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 20*k), offset=c(0.03,50), cutoff=k/2)
    par(mfrow=c(1,2))
    plot(mesh, asp=1)
    plot(mycoords, add=TRUE, col=2, pch=16, cex=0.5)
    plot(trajectory)
    points(mycoords, pch=16, cex=0.7, col=2)
}

mesh.polygons <- vector('list', nrow(mesh$graph$tv))
for(i in 1:nrow(mesh$graph$tv)){
    mesh.polygons[[i]] <- SPL.xyz((mesh$loc[mesh$graph$tv[i,],-3])[1,],
    (mesh$loc[mesh$graph$tv[i,],-3])[2,], (mesh$loc[mesh$graph$tv[i,],-3])[3,], ID=i)
}
mesh.polygons <- SpatialPolygons(mesh.polygons)

Ypos <- data.frame(speed=X$speed, synced_time=X$synced_time, coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(lead = lead(coords)) %>% head(-1) %>%
    ## compute the unit vector in the direction of the next point
    mutate(dir      = map2(coords, lead, function(x, y) (y-x))) %>%
    mutate(SPl      = map2(coords, lead, SPL.xy)) %>%
    mutate(SPcoords = lapply(coords, function(x) SpatialPoints(matrix(x, ncol=2)))) %>%
    mutate(SPlead   = lapply(lead, function(x) SpatialPoints(matrix(x, ncol=2)))) %>%
    mutate(Li       = map2(coords, lead, function(x, y) as.numeric(sqrt(sum((y-x)^2))))) %>%
    mutate(endpoints = pmap(list(SPcoords, SPlead, coords, lead), function(x, y, w, z) {
        triangle.coords <- over(x, mesh.polygons)
        triangle.lead   <- over(y, mesh.polygons)
        if(triangle.coords == triangle.lead){
            print(paste("nocrossing"))
            return(rbind(w, z))
        }
        if(triangle.coords != triangle.lead){
            intersection.point <- NA
            coords.triangle.index <- triangle.coords          
            pts1 <- line.line.intersection(P1 = w,
                                           P2 = z,
                                           P3 = ((mesh.polygons@polygons[coords.triangle.index])[[1]])@Polygons[[1]]@coords[1,],
                                           P4 = ((mesh.polygons@polygons[coords.triangle.index])[[1]])@Polygons[[1]]@coords[2,])
            ## print(paste("pts1: ", pts1))
            if(!is.logical(na.omit(pts1))) intersection.point <- pts1
            pts2 <- line.line.intersection(P1 = w,
                                           P2 = z,
                                           P3 = ((mesh.polygons@polygons[coords.triangle.index])[[1]])@Polygons[[1]]@coords[1,],
                                           P4 = ((mesh.polygons@polygons[coords.triangle.index])[[1]])@Polygons[[1]]@coords[3,])
            ## print(paste("pts2: ", pts2))
            if(!is.logical(na.omit(pts2))) intersection.point <- pts2
            pts3 <- line.line.intersection(P1 = w,
                                           P2 = z,
                                           P3 = ((mesh.polygons@polygons[coords.triangle.index])[[1]])@Polygons[[1]]@coords[2,],
                                           P4 = ((mesh.polygons@polygons[coords.triangle.index])[[1]])@Polygons[[1]]@coords[3,])
            ## print(paste("pts3: ", pts3))
            if(!is.logical(na.omit(pts3))) intersection.point <- pts3
            print(paste("Intersection points is: ", round(intersection.point[1], 2), ", ", round(intersection.point[2], 2)))
            return(rbind(w, intersection.point, z))            
        }
    }))

save(Ypos, file="Ypos.RData")

triangles.coords <- over(SpatialPoints(do.call("rbind", Ypos$coords)), mesh.polygons)
triangles.lead   <- over(SpatialPoints(do.call("rbind", Ypos$lead)), mesh.polygons)


diff.triangles   <- which(as.logical(apply(cbind(triangles.coords, triangles.lead),1,function(x) x[1]==x[2])==FALSE))

for(i in 1:length(diff.triangles)){
    plot(mesh.polygons[as.numeric(triangles.coords[diff.triangles[i]]),], asp=1)
    plot(mesh.polygons[as.numeric(triangles.lead[diff.triangles[i]]),], add=TRUE)
    Sys.sleep(1)
}


Ypos <- Ypos %>% mutate(index=1:nrow(Ypos)) %>% 
    mutate(slines = pmap(list(coords, lead, index), function(x, y, z) Lines(Line(rbind(x, y)), ID=as.character(z))))

## o <- SpatialLines(Ypos$slines)

o <- raster::intersect(SpatialLines(Ypos$slines), mesh.polygons)  

o <- over(SpatialLines(Ypos$slines), mesh.polygons)

raster::intersect(mesh.polygons,
                  SpatialLines(list(Lines(Line(rbind(Ypos$coords[[1]], Ypos$lead[[1]])), ID="1"))))  

SpatialLines(list(Lines(Line(cbind(Ypos$coords[[1]], Ypos$lead[[1]])), ID="1"),
                  Lines(Line(cbind(Ypos$coords[[4]], Ypos$lead[[4]])), ID="4"))) %over% y=mesh.polygons

plot(mesh, asp=1)

plot(raster::intersect(mesh.polygons,
                       SpatialLines(list(Lines(Line(rbind(Ypos$coords[[1]], Ypos$lead[[1]])), ID="1")))), col=2, ad=FALSE)

plot(SpatialLines(list(Lines(Line(rbind(Ypos$coords[[1]], Ypos$lead[[1]])), ID="1"))), add=TRUE, lwd=2)



if(FALSE){
    Ypos <- Ypos %>% as_tibble %>%
        mutate(midpoints = lapply(endpoints, function(endpoints){
            out <- matrix(apply(endpoints, 2,
                                function(col.vec) {
                                    out <- NULL
                                    for(i in 1:(length(col.vec)-1)) out[i] <- (col.vec[i]+col.vec[i+1])/2
                                    return(as.numeric(out))
                                }), nrow=(nrow(endpoints)-1))
            return(out)
        })) %>%
        mutate(midpoints.weights = map2(midpoints, sLi, function(midpoints, sLi){
            out <- matrix(rep(sLi, nrow(midpoints)), ncol=1)
            return(out)
        }))
}




if(FALSE){
    ((mesh.polygons@polygons[597])[[1]])@Polygons[[1]]@coords[2,]

    plot(mesh.polygons[597,])
    plot(mesh.polygons[5553,], add=TRUE)
    plot(Ypos$SPcoords[[3]], add=TRUE)
    plot(Ypos$SPlead[[3]], add=TRUE)
    pts <- SpatialPoints(matrix(line.line.intersection(P1 = Ypos$coords[[3]],
                                                       P2 = Ypos$lead[[3]],
                                                       P3 = ((mesh.polygons@polygons[597])[[1]])@Polygons[[1]]@coords[2,],
                                                       P4 = ((mesh.polygons@polygons[597])[[1]])@Polygons[[1]]@coords[3,]), nrow=1))

    function(x, y, z, w){
        c(line.line.intersection(P1 = x, P2 = y, P3 = w, P4 = z)
        }

        plot(pts, add=TRUE, pch=16)


        tmp1 <- lapply(Ypos$SPcoords, function(y) lapply(mesh.polygons, function(x) gIntersects(y, x)))
        tmp2 <- lapply(Ypos$SPlead, function(y) lapply(mesh.polygons, function(x) gIntersects(y, x)))


        plot(Ypos$SPl[[1]], ylim=c(0,100), xlim=c(0,100), asp=1)
        for(i in 2:nrow(Ypos)){
            plot(Ypos$SPl[[i]], ylim=c(0,100), xlim=c(0,100), add=TRUE)
        }

        do.call("rbind", as_tibble(Ypos)$coords)[1,]

        do.call("rbind", as_tibble(Ypos)$lead)[1,]




        P1 <- Polygon(rbind(c(1,0),c(0,1),c(1,0)))
        Ps1 = Polygons(list(P1), ID = "a")
        plot(SpatialPolygons(list(Ps1)))
        ##
        ## create a regular mesh of points between all line segments defined
        ## by the trajectory of the mouse (positional data)
        ##
        Ypos <- data.frame(speed=X$speed, synced_time=X$synced_time, coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
            mutate(lead = lead(coords)) %>% head(-1) %>%
            ## compute the unit vector in the direction of the next point
            mutate(dir    = map2(coords, lead, function(x, y) (y-x))) %>% ##/sqrt(sum((y-x)^2))
            ## compute the length of the lines joining the points
            mutate(Li     = map2(coords, lead, function(x, y) as.numeric(sqrt(sum((y-x)^2)))))

        ##
        ## set a threshold value to be used as a rule for the segmentation (k = ceiling of |Li|/h)
        ##
        h <- quantile(unlist(Ypos$Li), 0.00001) 
        0.001

        Ypos <- Ypos %>%
            ## the final number of line segments that compose each line in s (if 1 then line is not segmented)
            mutate(k   = ceiling(as.numeric(Li)/h)) %>%
            ## and the length of the line segments is given by sLi
            mutate(sLi = map2(Li, k, function(Li, k) Li/k)) %>%
            mutate(endpoints= pmap(list(coords, dir, k), function(coords, dir, k){
                ## -------------------------------------------------------------
                ## compute the coordinates of the endpoints of all line segments
                ## if k==1 (1 line segment),
                ## Test: code must return the actual coordinates of the
                ## initial line segment. DONE
                ## -------------------------------------------------------------
                endpoints <- coords # matrix(nrow=k+1, ncol=2)
                for(i in 1:k) endpoints <- rbind(endpoints, coords + i*(dir/k))
                return(unname(endpoints))
            })) %>% as_tibble %>%
            mutate(midpoints = lapply(endpoints, function(endpoints){
                out <- matrix(apply(endpoints, 2,
                                    function(col.vec) {
                                        out <- NULL
                                        for(i in 1:(length(col.vec)-1)) out[i] <- (col.vec[i]+col.vec[i+1])/2
                                        return(as.numeric(out))
                                    }), nrow=(nrow(endpoints)-1))
                return(out)
            })) %>%
            mutate(midpoints.weights = map2(midpoints, sLi, function(midpoints, sLi){
                out <- matrix(rep(sLi, nrow(midpoints)), ncol=1)
                return(out)
            }))

        data <- list(Ypos=Ypos, Y=Y)


        ## animation?
        if(FALSE){
            
            o <- SPL.xyz((mesh$loc[mesh$graph$tv[1,],-3])[1,], (mesh$loc[mesh$graph$tv[1,],-3])[2,], (mesh$loc[mesh$graph$tv[1,],-3])[3,])
            plot(o, ylim=c(-20,120), xlim=c(-20,120)) 
            for(i in 2:nrow(mesh$graph$tv)){
                o <- SPL.xyz((mesh$loc[mesh$graph$tv[i,],-3])[1,], (mesh$loc[mesh$graph$tv[i,],-3])[2,], (mesh$loc[mesh$graph$tv[i,],-3])[3,])
                plot(o, add=TRUE)
            }
            
        }
}







    
