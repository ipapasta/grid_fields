library(dplyr)
library(readxl)
library(sp)
##
## load data
##
if(getwd()=="/Users/ipapasta/Documents/Research/Software/R/grid_fields/R/Temporal_modulation"){
X <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields_old/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
Y <- read_excel("/Users/ipapasta/Documents/Research/Software/R/grid_fields_old/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
## X <- read.csv("/home/ipapasta/Software/R/grid_fields/mouse_14_0516_cell_13/position.pkl.csv", stringsAsFactors=FALSE)
## Y <- read.csv("/home/ipapasta/Software/R/grid_fields/mouse_14_0516_cell_13/spatial_firing.pkl.csv")
}

## mouse 5 cell 6
if(getwd()=="/Users/ipapasta/Documents/Research/Software/R/grid_fields/R/Oscillating"){
X <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
Y <- read_excel("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
}


if(FALSE){
    ## mouse 12_0410 cell 1
    if(getwd()=="/Users/ipapasta/Documents/Research/Software/R/grid_fields/R/Oscillating"){
        X <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse12_0410_cell1/position.pkl.csv", stringsAsFactors=FALSE)
        Y <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse12_0410_cell1/spatial_firing.pkl.csv")
    }


    ## mouse_14_0516_cell_13
    if(getwd()=="/Users/ipapasta/Documents/Research/Software/R/grid_fields/R/Oscillating"){
        X <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse_14_0516_cell_13/position.pkl.csv", stringsAsFactors=FALSE)
        Y <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse_14_0516_cell_13/spatial_firing.pkl.csv")
    }

    ## mouse12_0410_cell1
    if(getwd()=="/Users/ipapasta/Documents/Research/Software/R/grid_fields/R/Oscillating"){
        X <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse12_0410_cell1/position.pkl.csv", stringsAsFactors=FALSE)
        Y <- read.csv("/Users/ipapasta/Documents/Research/Software/R/grid_fields/mouse12_0410_cell1/spatial_firing.pkl.csv")
    }
}
## y <- Y$firing_times
## x <- X$synced_time
mycoords   <- SpatialPoints(cbind(Y$position_x, Y$position_y))
## 
## scale positions to [0,1] \times [0,1]
## 

## X$position_x_pixels <- X$position_x_pixels/mx
## X$position_y_pixels <- X$position_y_pixels/my
## Y$position_x_pixels <- Y$position_x_pixels/mx
## Y$position_y_pixels <- Y$position_y_pixels/my

if(TRUE){
    ##
    ## firing events
    ##
    mycoordsDF <- SpatialPointsDataFrame(coords=mycoords, data=data.frame(hd=Y$hd))
    ##
    ## trajectory
    ##
    Pl   <- Polygon(cbind(X$position_x, X$position_y))
    ID   <- "[0,1]x[0,1]"
    Pls  <- Polygons(list(Pl), ID=ID)
    SPls <- SpatialPolygons(list(Pls))
    df   <- data.frame(value=1, row.names=ID)
    SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         # str(df)
    trajectory <- SpatialPolygonsDataFrame(SPls, df) 
    if(FALSE){
        plot(trajectory)
        points(mycoords, pch=16)
        plot(trajectory)
        points(mycoords, pch=16, cex=0.7, col=2)
    }
}


if(FALSE){
    load("spatial_firing.Rdata")
    load("trajectory.Rdata")

    posxf = spatial_firing$position_x   
    posyf = spatial_firing$position_y   

    posx <- trajectory$position_x
    posy <- trajectory$position_y

    ind <- 140
    ## plot(posx[[ind]], posy[[ind]], type="l", asp=1)
    Pl   <- Polygon(na.omit(cbind(posx[[ind]], posy[[ind]])))
    ID   <- "[0,1]x[0,1]"
    Pls  <- Polygons(list(Pl), ID=ID)
    SPls <- SpatialPolygons(list(Pls))
    df   <- data.frame(value=1, row.names=ID)
    SPDF       <- SpatialPolygonsDataFrame(SPls, df)                                         # str(df)
    traj       <- SpatialPolygonsDataFrame(SPls, df)
    mycoords   <- SpatialPoints(na.omit(cbind(posxf[[ind]], posyf[[ind]])))

    plot(traj)
    points(mycoords, col=2, pch=16)
    
}


if(FALSE){
    library(sp)   
    ## 
    load("trajectory.Rdata")
    load("spatial_firing.Rdata")
    ## 
    X = data.frame( position_x = trajectory$position_x[110][[1]],    # List element 110
                   position_y = trajectory$position_y[110][[1]] )
    ## rm(trajectory)

    Y = data.frame( position_x = spatial_firing$position_x[348][[1]],    # List element 348
                   position_y = spatial_firing$position_y[348][[1]],
                   hd = spatial_firing$hd[348][[1]] )
    ## rm(spatial_firing)                                        # Firing events
    mycoords   = SpatialPoints(cbind(Y$position_x, Y$position_y))
    mycoordsDF = SpatialPointsDataFrame(coords=mycoords, data=data.frame(hd=Y$hd))

                                        # Trajectory
    session_id.trajectory <- trajectory$session_id
    session_id.spatial_firing <- spatial_firing$session_id

    
    Pl   = Polygon(cbind(X$position_x, X$position_y))
    ID   = "[0,1]x[0,1]"
    Pls  = Polygons(list(Pl), ID=ID)
    SPls = SpatialPolygons(list(Pls))
    df   = data.frame(value=1, row.names=ID)
    SPDF  = SpatialPolygonsDataFrame(SPls, df) 
    trajectory = SpatialPolygonsDataFrame(SPls, df) 

    plot(trajectory)
    points(mycoords, pch=16, col=2)

    ## for(i in 1:)
}    
