## Construct meshes for all data sets

## offset : it specifies how much the domain will be extended in the
## outer and inner part. If negative it is interpreted as a factor
## relative to the approximate data diameter. If positive it is the
## extension distance on same scale unit to the coordinates provided.

## cutoff : it specifies the minimum distance allowed between
## points. It means that if the distance between two points is less
## than the supplied value then they are replaced by a single
## vertex. Its very useful in case of clustered data points because it
## avoids building many small triangles arround clustered points.

## min.angle: it specifies the minimum internal angle of the
## triangles. This could be a two-dimensional vector with the same
## meaning of the others. Take in mind that we would like to have a
## mesh with triangles as regular as possible.

## n, interior : The argument n is the initial number of points on the
## extended boundary. The interior is a list of segments to specify
## interior constraints, each one of inla.mesh.segment class.


## ----------------------------------------------------------------------------
## From spatial_firing, filter only those elements for which session_id appears
## in trajectory object.
## ----------------------------------------------------------------------------

library(sp)
library(tidyverse)
library(INLA)
## 
load("trajectory.RData")
load("spatial_firing_meshes.RData")
## str(trajectory)
## str(spatial_firing)
if(getwd()=="/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating"){
    spatial_firing  <- as_tibble(spatial_firing) %>% dplyr::filter(grid_score>0.8) %>% arrange(desc(session_id)) %>%
        mutate(Y = map2(position_x, position_y, function(x, y) na.omit(data.frame( position_x = x, position_y = y)))) %>%
        mutate(nrowY = unlist(lapply(Y, function(x) dim(as.matrix(x))[1]))) %>%    
        group_by(session_id) %>% mutate(cond = as.numeric(nrowY == min(nrowY))) %>% ungroup() %>% 
        dplyr::filter(cond==1) %>% arrange(session_id) %>%
        mutate(mycoords = lapply(Y, function(x) SpatialPoints(x))) %>% filter(any(nrowY==min(nrowY)))
}
## ids             <- unique(spatial_firing$session_id) %>% as.character

df   <- data.frame(value=1, row.names="[0,1]x[0,1]")
trajectory      <- as_tibble(trajectory) %>% dplyr::filter(session_id %in% spatial_firing$session_id) %>%    
    arrange(desc(session_id)) %>%
    mutate(X = map2(position_x, position_y, function(x,y) na.omit(data.frame( position_x=x, position_y=y)))) %>%
    mutate(mycoords = lapply(X, function(x) SpatialPoints(x))) %>%
    mutate(Pl = lapply(X, function(x) Polygon(as.matrix(x)))) %>%
    mutate(Pls  = lapply(Pl, function(x) Polygons(list(x), ID="[0,1]x[0,1]"))) %>% 
    mutate(SPls = lapply(Pls, function(x) SpatialPolygons(list(x))))  %>% 
    mutate(trajectory = lapply(SPls, function(x) SpatialPolygonsDataFrame(x, df)))

## cbind(trajectory$session_id, spatial_firing$session_id)
## create meshes for all data
if(getwd()=="/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating"){
    k <- 3
    spatial_firing <- spatial_firing %>% 
        mutate(mesh = lapply(mycoords, function(x) inla.mesh.2d(x, max.edge=c(k, 10*k), offset=c(0.03,50), cutoff=k/2)))
    save(spatial_firing, file="spatial_firing_meshes.RData")
}



##
## create a regular mesh of points between all line segments defined
## by the trajectory of the mouse (positional data)
##

trajectory  <- trajectory %>%
    mutate(Ypos = lapply(X, function(x){
        ## missing speed?
        ## uniform(h)
        ## min position_x[[1]] difference 0.007736023.
        h <- 0.0001
        o <- data.frame( coords=I(lapply(as.list(apply(cbind(x$position_x, x$position_y),1, as.list)), unlist))) %>%
            mutate(lead = lead(coords)) %>% head(-1) %>%
            ## compute the unit vector in the direction of the next point
            mutate(dir    = map2(coords, lead, function(x, y) (y-x))) %>% ##/sqrt(sum((y-x)^2))
            ## compute the length of the lines joining the points
            mutate(Li     = map2(coords, lead, function(x, y) as.numeric(sqrt(sum((y-x)^2))))) %>%
            ##
            ## set a threshold value to be used as a rule for the segmentation (k = ceiling of |Li|/h)
            ##
            ## Ypos <- Ypos %>%
            ## the final number of line segments that compose each line in s (if 1 then line is not segmented)
            ## o <- o %>%
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
        ## return(Ypos)
    }))

## par(mfrow=c(1,2))
if(FALSE){
    for(i in 1:nrow(spatial_firing)){
        Sys.sleep(3)
        plot(spatial_firing$mesh[[i]], asp=1)
        Sys.sleep(3)
        plot(spatial_firing$mycoords[[i]], add=TRUE, col=2, pch=16)
    }
}

if(FALSE){
    save(trajectory, file="trajectory_data.RData")
    save(spatial_firing, file="spatial_firing_data.RData")
}







