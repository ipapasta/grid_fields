## most recent implementation by Graeme
args = commandArgs(trailingOnly=TRUE)
## args[1] = 20
time_splits = as.numeric(args[1])
library(tidyverse)
library(dplyr)
## library(readxl)
library(sp)
simulation.hpp <- FALSE
##
## load data
##
## --------------
## Simulated data
## --------------
if(FALSE) {
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
    ## load("spatial_firing_all_mice_hist.Rda")
    ## load("trajectory_all_mice_hist(1).Rda")
    grid_cells_index.firing <- which(spatial_firing$grid_score>0.8)
    grid_cells_session.id   <- spatial_firing$session_id[grid_cells_index.firing]
    grid_cells_index.position <- which(trajectory_all_mice_hist$session_id %in% grid_cells_session.id)
    ## check they match
    cbind(grid_cells_session.id, trajectory_all_mice_hist$session_id[grid_cells_index.position])
    df <- data.frame(session=(grid_cells_session.id),
                     index.firing=(grid_cells_index.firing),
                     index.position=(grid_cells_index.position))
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
    dat <- data_extract(df$index.position[5], df$index.firing[5],
                        trajectory=trajectory_all_mice_hist,
                        firing=spatial_firing)
    Xraw <- dat$X
    ## set.seed(111086)
    quant <- 1
    dat$X <- dat$X %>% mutate(index.CV = rep(-1, nrow(dat$X))) %>%
        dplyr::filter(synced_time <= quantile(dat$X$synced_time, quant))
    ## ------------------------------------
    ## SIMULATE HOMOGENEOUS Poisson Process
    ## ------------------------------------
    if(!simulation.hpp){
        dat$Y <- dat$Y %>% mutate(index.CV = rep(-1, nrow(dat$Y))) %>%
            dplyr::filter(firing_times <= quantile(dat$X$synced_time, quant))
    }else{
        dat$Y <- path.hpp(lambda=0.1, x=dat$X$position_x, y=dat$X$position_y, time=dat$X$synced_time, hd=dat$X$hd)
        dat$Y <- dat$Y %>% mutate(index.CV = rep(-1, nrow(dat$Y))) %>%
            dplyr::filter(firing_times <= quantile(dat$X$synced_time, quant))
    }
    if(FALSE){
        lambda.CV <- exp(.8)*mean(diff(dat$X$synced_time))    ## (exp(.5) for 30 segments)
        ## lambda.CV <- .5
        Z <- 0
        counter <- 0
        while(Z[length(Z)] < max(dat$X$synced_time)){
            print(counter)
            Z <- c(Z, Z[length(Z)]+rexp(1,lambda.CV))
            dat$X$index.CV[which(dat$X$synced_time < Z[counter+2] & dat$X$synced_time > Z[counter+1])] <- counter + 1
            dat$Y$index.CV[which(dat$Y$firing_times < Z[counter+2] & dat$Y$firing_times > Z[counter+1])] <- counter + 1
            counter <- counter + 1
        }
        K <- max(dat$X$index.CV)
        train.index <- sort(sample(1:K, size = floor(K/2), replace=FALSE))
    }
    if(TRUE){
        Z <- 0
        counter <- 0
        while(Z[length(Z)] < max(dat$X$synced_time)){
            print(counter)
            Z <- c(Z, Z[length(Z)]+time_splits)
            dat$X$index.CV[which(dat$X$synced_time < Z[counter+2] & dat$X$synced_time > Z[counter+1])] <- counter + 1
            dat$Y$index.CV[which(dat$Y$firing_times < Z[counter+2] & dat$Y$firing_times > Z[counter+1])] <- counter + 1
            counter <- counter + 1
        }
        K <- max(dat$X$index.CV)
        train.index <- seq(1, K , by = 2)
        char.to.save <- as.character(round(max(dat$X$synced_time)/K))
    }
    X.train <- dat$X %>% dplyr::filter(index.CV %in% train.index)
    Y.train <- dat$Y %>% dplyr::filter(index.CV %in% train.index)
    X.test <- dat$X %>% dplyr::filter(!(index.CV %in% train.index))
    Y.test <- dat$Y %>% dplyr::filter(!(index.CV %in% train.index))
    if(FALSE){
        clumps <- split(unique(Y.test$index.CV), cumsum(c(1, diff(unique(Y.test$index.CV)) != 1)))
        obs.firings <- sapply(1:length(clumps), function(i) nrow(Y.test %>% filter(index.CV %in% clumps[[i]])))
        plot((obs.firings))
    }
    ##
    ## Firing events
    ## 
    mycoords.CV       <- SpatialPoints(cbind(Y.train$position_x, Y.train$position_y))        
    ##
    ## trajectory
    ##
    Pl   <- Polygon(cbind(X.train$position_x, X.train$position_y))
    Pl   <- Polygon(cbind(dat$X$position_x, dat$X$position_y))
    ID   <- "[0,1]x[0,1]"
    Pls  <- Polygons(list(Pl), ID=ID)
    SPls <- SpatialPolygons(list(Pls))
    df   <- data.frame(value=1, row.names=ID)
    SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         # str(df)
    trajectory <- SpatialPolygonsDataFrame(SPls, df)
                                        #plot(trajectory)
}


if(FALSE){



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

        dat <- data_extract(df$index.position[5], df$index.firing[5],
                            trajectory=trajectory_all_mice_hist,
                            firing=spatial_firing)

        set.seed(111086)
        dat$X <- dat$X %>% mutate(index.CV = rep(-1, nrow(dat$X)))
        dat$Y <- dat$Y %>% mutate(index.CV = rep(-1, nrow(dat$Y)))
        lambda.CV <- mean(diff(dat$X$synced_time))
        Z <- 0
        counter <- 0
        while(Z[length(Z)] < max(dat$X$synced_time)){
            print(counter)
            Z <- c(Z, Z[length(Z)]+rexp(1,lambda.CV))
            dat$X$index.CV[which(dat$X$synced_time < Z[counter+2] & dat$X$synced_time > Z[counter+1])] <- counter + 1
            dat$Y$index.CV[which(dat$Y$firing_times < Z[counter+2] & dat$Y$firing_times > Z[counter+1])] <- counter + 1
            counter <- counter + 1
        }
        K <- max(dat$X$index.CV)
        train.index <- sort(sample(1:K, size = floor(K/2), replace=FALSE))
        X.train <- dat$X %>% dplyr::filter(index.CV %in% train.index)
        Y.train <- dat$Y %>% dplyr::filter(index.CV %in% train.index)

        X.test <- dat$X %>% dplyr::filter(!(index.CV %in% train.index))
        Y.test <- dat$Y %>% dplyr::filter(!(index.CV %in% train.index))
        
        ##
        ## Firing events
        ## 
        mycoords.CV       <- SpatialPoints(cbind(Y.train$position_x, Y.train$position_y))

        ##
        ## trajectory
        ##
        Pl   <- Polygon(cbind(X.train$position_x, X.train$position_y))
        ID   <- "[0,1]x[0,1]"
        Pls  <- Polygons(list(Pl), ID=ID)
        SPls <- SpatialPolygons(list(Pls))
        df   <- data.frame(value=1, row.names=ID)
        SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         # str(df)
        trajectory <- SpatialPolygonsDataFrame(SPls, df)

        plot(trajectory)
    }
}
