## ---------------
## Cyclic B-splines
## ---------------
hd.bsp.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    sigma   <- exp(theta2)
    tausq   <- 1/(4*(sigma^2)*(kappa^3))
    fem       <- inla.mesh.fem(mesh=mesh, order = 2)
    tausq*((kappa^4)*(fem$c0) + (2*(kappa^2)*(fem$g1)) + fem$g2)
}

## theta <- seq(0, 2*pi, len=10)
## mesh.cyclic <- inla.mesh.1d(theta, boundary="cyclic", degree=1)
## fem <- inla.mesh.1d.fem(mesh.cyclic)

