library(arrow)
library(tidyverse)
library(dplyr)
library(purrr)
library(readxl)
library(sp)
load("/data/spatial_firing_all_mice_hist.Rda")
load("data/trajectory_all_mice_hist(1).Rda")
## 
grid_cells_index.firing <- which(spatial_firing$grid_score>0.8)
grid_cells_session.id   <- spatial_firing$session_id[grid_cells_index.firing]
grid_cells_index.position <- which(trajectory_all_mice_hist$session_id %in% grid_cells_session.id)
## check they match
cbind(grid_cells_session.id, trajectory_all_mice_hist$session_id[grid_cells_index.position])
## 
df <- data.frame(session=(grid_cells_session.id),
                 index.firing=(grid_cells_index.firing),
                 index.position=(grid_cells_index.position))
## 
data_extract <- function(index.position, index.firing, trajectory, firing){
    ## X: trajectory
    X <- data.frame(synced_time=as.numeric(trajectory$synced_time[[index.position]]),
                    position_x = trajectory$position_x[[index.position]],
                    position_y = trajectory$position_y[[index.position]],
                    hd = (trajectory$hd[[index.position]] + 180)*(pi/180))
    ## Y: firing events
    Y <- data.frame(firing_times=as.numeric(firing$firing_times[[index.firing]])/(30*1000),
                    position_x = firing$position_x[[index.firing]],
                    position_y = firing$position_y[[index.firing]],
                    hd = (firing$hd[[index.firing]] + 180)*(pi/180))
    list.data <- list(X=X, Y=Y)
}
## 
dat <- data_extract(df$index.position[5], df$index.firing[5],
                    trajectory=trajectory_all_mice_hist,
                    firing=spatial_firing)
## 
X <- dat$X
Y <- dat$Y

##
## Firing events and trajectory
##
ID   <- "[0,1]x[0,1]"
spikes       <- SpatialPoints(cbind(Y$position_x, Y$position_y))
Pl           <- Polygon(cbind(X$position_x, X$position_y))
Pl.Omega     <- Polygon(expand.grid(c(min(X$position_x),max(X$position_x)), c(min(X$position_y),max(X$position_y)))[c(1,2,4,3),])
ID.Omega     <- "Omega"
Pls.Omega    <- Polygons(list(Pl.Omega), ID=ID.Omega)
SPls.Omega   <- SpatialPolygons(list(Pls.Omega))
Pls  <- Polygons(list(Pl), ID=ID)
SPls <- SpatialPolygons(list(Pls))
df   <- data.frame(value=1, row.names=ID)
SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         
trajectory <- SpatialPolygonsDataFrame(SPls, df) 


theta_precession <- FALSE
if(theta_precession){
    ## install.packages('arrow')
    theta_position <- read_feather(file.path(getwd(), "data", "position_theta_cluster_7.feather")) %>%
        dplyr::mutate(hd = theta_angle_3 + pi) %>% dplyr::select(-theta_angle_3)
    theta_firing <- read_feather(file.path(getwd(), "data/", "spatial_firing_theta_cluster_7.feather")) %>%
            dplyr::mutate(firing_times = firing_times/30000) %>%
        dplyr::mutate(hd = theta_angle + pi) %>% dplyr::select(-theta_angle)
    X <- theta_position
    Y <- theta_firing
    spikes       <- SpatialPoints(cbind(Y$position_x, Y$position_y))
    Pl           <- Polygon(cbind(X$position_x, X$position_y))
    Pl.Omega     <- Polygon(expand.grid(c(min(X$position_x),max(X$position_x)), c(min(X$position_y),max(X$position_y)))[c(1,2,4,3),])
    ID.Omega     <- "Omega"
    Pls.Omega    <- Polygons(list(Pl.Omega), ID=ID)
    SPls.Omega   <- SpatialPolygons(list(Pls.Omega))
    ID   <- "[0,1]x[0,1]"
    Pls  <- Polygons(list(Pl), ID=ID)
    SPls <- SpatialPolygons(list(Pls))
    df   <- data.frame(value=1, row.names=ID)
    SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         
    trajectory <- SpatialPolygonsDataFrame(SPls, df)
    ## domain <- SpatialPolygons(list(, "s1"))
}
