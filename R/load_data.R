library(dplyr)
library(readxl)
library(sp)

##
## load data
## 

X <- read.csv("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
Y <- read_excel("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")

y <- Y$firing_times
x <- X$synced_time

## 
## scale positions to [0,1] \times [0,1]
## 


## X$position_x_pixels <- X$position_x_pixels/mx
## X$position_y_pixels <- X$position_y_pixels/my
## Y$position_x_pixels <- Y$position_x_pixels/mx
## Y$position_y_pixels <- Y$position_y_pixels/my



(Y %>% dplyr::select(position_x, position_y))[1,]



X$position_x_pixels <- X$position_x_pixels
X$position_y_pixels <- X$position_y_pixels
Y$position_x_pixels <- Y$position_x_pixels
Y$position_y_pixels <- Y$position_y_pixels

##
## firing events
##
mycoords   <- SpatialPoints(cbind(Y$position_x, Y$position_y))
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
    points(mycoords, pch=16, col=2)
}
