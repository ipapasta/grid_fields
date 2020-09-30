##
## Input: 
##      - theta1: correlation length parameter. Note that
##                this parameter does not have the same interpretation
##                as in the usual Matern covariance family due to the oscillation
##                of the process. However, it can still be used to describe
##                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
##                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
##      - theta2: scale parameter that controls the variance of the field.
##                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
##                standard deviation of the process but is related to.
##      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
##                -1 < phi < 0: oscillating 
##                 0 < phi < 1: overdampened oscillating 
## Output:
##      - precision matrix for the weights at the mesh vertices 
##        
## 

temp.precision <- function(theta, mesh, o=2){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    sigma   <- exp(theta2)
    tausq   <- 1/(4*(sigma^2)*(kappa^3))
    ## ---------------------------------------
    o       <- inla.mesh.fem(mesh=mesh, order = o)
    ## 
    tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*(o$g1)) + o$g2)
    ## 
}

if(FALSE){
    tmp  <- NULL
    sim  <- matrix(nrow=1000, ncol=1000)
    tmp2 <- matrix(nrow=1000, ncol=1000)
    y <- seq(0.01, 1000, len=100)
    theta.Z <- theta[4:5]    
    QZ     <- temp.precision(theta=theta.Z, mesh=mesh1d)
    logpar    <- theta.Z
    mesh1d    <- inla.mesh.1d (y)
    for(i in 1:1000){
        ## logpar    <- c(log(200), log(1))
        ## Q         <- temp.precision(theta=logpar, mesh=mesh1d, o=2)
        x         <- inla.qsample(Q=QZ[1:100,1:100])
        tmp2[i,]  <- x
        A         <- inla.spde.make.A(mesh=mesh1d, loc=y)
        sim[i,]   <- as.numeric(A %*% x)
        tmp[i]    <- var(sim[i,])
    }

    o <- inla.mesh.basis(mesh1d)


    matplot(t(sim), type="l", col=1)

    plot(seq(0.01, 1000, len=1000), as.numeric(sim[i,]), type="l")



    mesh1d  <- inla.mesh.1d(seq(0.01, 99, len=100), offset = c(.1, 1), n = c(1, 2), max.edge = c(0.05, 0.5))
    mesh1d  <- inla.mesh.1d(seq(0.01, 99, len=100), max.edge=c(1, 3), cutoff=0.5)
    ggplot()  + gg(mesh1d, shape = "|", size = 5) + theme_classic() + xlim(0,10)


    ggplot() + gg(mesh) + coord_fixed() + xlim(0,100) + ylim(0,100) + theme_classic()


    n = 100
    loc = matrix(runif(n*2), n, 2)
    mesh = inla.mesh.2d(loc, max.edge=0.05)
    basis = inla.mesh.basis(mesh, n=c(4,5))

    proj = inla.mesh.projector(mesh)
    image(proj$x, proj$y, inla.mesh.project(proj, basis[,7]))
}
