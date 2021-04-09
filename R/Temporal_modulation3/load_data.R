library(tidyverse)
library(dplyr)
library(readxl)
library(sp)
##
## load data
##
## --------------
## Simulated data
## --------------
if(FALSE){
    X <- read.csv("data/Simulated_data/simulated_cn_same_30var_v_grid5trials1_simulated/session.csv",
                  stringsAsFactors=FALSE)
    Y <- read.csv("data/Simulated_data/simulated_cn_same_30var_v_grid5trials1_simulated/1_firing_events.csv",
                  stringsAsFactors=FALSE)
    X <- X %>% mutate(hd=  (hd + 180)*(pi/180))
    Y <- Y %>% mutate(hd = (hd + 180)*(pi/180))
    Y$firing_times <- Y$firing_times/1000
    mycoords       <- SpatialPoints(cbind(Y$position_x, Y$position_y))

    ## 
    ## reduce size of dataset
    ##

    if(FALSE){
        max_time <- 300
        Y <- Y %>% filter(firing_times < max_time)
        X <- X %>% filter(synced_time < max_time)
    }

    mycoords       <- SpatialPoints(cbind(Y$position_x, Y$position_y))
}

experimental <- TRUE
if(experimental){
    load("data/spatial_firing_all_mice_hist.Rda")
    load("data/trajectory_all_mice_hist(1).Rda")

    grid_cells_index.firing <- which(spatial_firing$grid_score>0.8)
    grid_cells_session.id   <- spatial_firing$session_id[grid_cells_index.firing]
    grid_cells_index.position <- which(trajectory_all_mice_hist$session_id %in% grid_cells_session.id)

    ## check they match
    cbind(grid_cells_session.id, trajectory_all_mice_hist$session_id[grid_cells_index.position])
    
    df <- data.frame(session=(grid_cells_session.id),
                    index.firing=(grid_cells_index.firing),
                    index.position=(grid_cells_index.position))

    data_extract <- function(index.position, index.firing, trajectory, firing){
        X <- data.frame(synced_time=as.numeric(trajectory$synced_time[[index.position]]),
                   position_x = trajectory$position_x[[index.position]],
                   position_y = trajectory$position_y[[index.position]],
                   hd = (trajectory$hd[[index.position]] + 180)*(pi/180))
        Y <- data.frame(firing_times=as.numeric(firing$firing_times[[index.firing]])/(30*1000),
                   position_x = firing$position_x[[index.firing]],
                   position_y = firing$position_y[[index.firing]],
                   hd = (firing$hd[[index.firing]] + 180)*(pi/180))
        list.data <- list(X=X, Y=Y)
    }

    dat <- data_extract(df$index.position[5], df$index.firing[5],
                 trajectory=trajectory_all_mice_hist,
                 firing=spatial_firing)

    X <- dat$X
    Y <- dat$Y

    mycoords       <- SpatialPoints(cbind(Y$position_x, Y$position_y))
    ## plot(mycoords)
}
