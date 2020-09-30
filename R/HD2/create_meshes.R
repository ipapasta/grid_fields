create.mesh.given.data.path <- function(path){
    firing_events.data <- paste0(path, "/1_firing_events.csv")
    k <- 2
    Y <- read.csv(firing_events.data, stringsAsFactors=FALSE)
    coords   <- SpatialPoints(cbind(Y$position_x, Y$position_y))
    mesh <- inla.mesh.2d(coords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
    save(mesh, file=paste0(path, "/mesh3.RData"))
}
df <- tibble(dirs = list.dirs("/home/ipapasta/grid-fields/R/data/Simulated_data/")[-1])
o <- df %>% mutate(mesh = lapply(dirs, create.mesh.given.data.path))



