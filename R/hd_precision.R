## ---------------
## Cyclic B-splines
## ---------------
hd.bsp.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    ## phi     <- (1-exp(-theta3))/(1+exp(-theta3))
    ## sincpth <- sqrt(1-phi^2)/acos(phi)
    sigma   <- exp(theta2)
    ## variance needs adjustment
    tausq   <- (pi+(cosh(pi*kappa)*sinh(pi*kappa))/kappa)/
        ((sigma^2)*(2*kappa)^2*sinh(pi*kappa)^2)
    fem       <- inla.mesh.fem(mesh=mesh, order = 2)
    tausq*((kappa^4)*(fem$c0) + (2*(kappa^2)*(fem$g1)) + fem$g2)
}

## theta <- seq(0, 2*pi, len=10)
## mesh.cyclic <- inla.mesh.1d(theta, boundary="cyclic", degree=1)
## fem <- inla.mesh.1d.fem(mesh.cyclic)

