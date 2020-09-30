library(dplyr)
library(readxl)
library(sp)
##
## load data
##


## X <- read.csv("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
## Y <- read_excel("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
## X <- read.csv("/home/ipapasta/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
X <- read.csv("/home/ipapasta/Software/R/grid_fields/M12_2018-04-10_14-22-14_of/session.csv", stringsAsFactors=FALSE)
## Y <- read_excel("/home/ipapasta/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
Y <- read.csv("/home/ipapasta/Software/R/grid_fields/M12_2018-04-10_14-22-14_of/2_firing_events.csv", stringsAsFactors=FALSE)
X <- X %>% mutate(hd=  (hd + 180)*(pi/180))
Y <- Y %>% mutate(hd = (hd + 180)*(pi/180))



## X <- read.csv("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/M12_2018-04-10_14-22-14_of/session.csv", stringsAsFactors=FALSE)
## ## Y <- read_excel("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
## Y <- read.csv("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/M12_2018-04-10_14-22-14_of/2_firing_events.csv", stringsAsFactors=FALSE)

## X <- X %>% mutate(hd = (hd + 180)*(pi/180)) 
## Y <- Y %>% mutate(hd = (hd + 180)*(pi/180))


##
## Spatial coordinates
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
    points(mycoords, pch=16, cex=0.7, col=2)
}


## maybe remove data preprocessing in main script spde.osc.HD.R and
## include here


## request data format

## X <- read.csv("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
## Y <- read_excel("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
## X <- read.csv("/home/ipapasta/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_whole_session.csv", stringsAsFactors=FALSE)
## Y <- read_excel("/home/ipapasta/Software/R/grid_fields/mouse_5_cell_6/mouse5_spatial_data_firing_events.xlsx")
## X <- X %>% mutate(hd=  (hd + 180)*(pi/180))
## Y <- Y %>% mutate(hd = (hd + 180)*(pi/180))

