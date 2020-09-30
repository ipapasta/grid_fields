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

hd.precision <- function(theta, order){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    sigma   <- exp(theta2) #!!
    tausq   <- 1/(4*(sigma^2)*(kappa^3))
    ## ---------------------------------------
    ## next line of code needs amendment (circular harmonics)
    ## o         <- inla.mesh.fem(mesh=mesh, order = o)
    ##
    C         <- C.circular.harmonics(order, inverse=FALSE)
    C.inverse <- C.circular.harmonics(order, inverse=TRUE)
    G         <- G.circular.harmonics(order) 
    tausq*((kappa^4)*C + (2*(kappa^2)*G) + (G%*%(C.inverse%*%G)))
    ## 
}



if(FALSE){
    ## theta2  <- seq(-1,-10,len=100)
    ord <- 3
    theta2 <- 0.1
    theta1  <- seq(-10,10,len=100)
    ## theta <- c(.7, 1e-100)
    for(i in seq_along(theta1)){
        Q     <- hd.precision(theta=c(theta1[i], theta2), order=ord)
        plot(diag(Q), type="l")
        Sys.sleep(.1)
    }    
    ord <- 2
    tmp  <- NULL
    sim  <- matrix(nrow=1000, ncol=100)
    tmp2 <- matrix(nrow=1000, ncol=1+2*ord)
    y    <- seq(0.01, 2*pi, len=100)
    theta <- c(pi/2, -3)
    Q     <- hd.precision(theta=theta, order=ord)
    logpar    <- theta
    for(i in 1:100){
        ## logpar    <- c(log(200), log(1))
        ## Q         <- temp.precision(theta=logpar, mesh=mesh1d, o=2)
        x         <- inla.qsample(Q=Q)
        tmp2[i,]  <- x
        A         <- circular.make.A(hd.vec=y, order=ord)
        sim[i,]   <- as.numeric(A %*% x)
        tmp[i]    <- var(sim[i,])
    }
    matplot(y, (t(sim[1:10,])), type="l", col=1)
}

if(FALSE){
    tmp  <- NULL
    sim  <- matrix(nrow=1000, ncol=100)
    tmp2 <- matrix(nrow=1000, ncol=7)
    y    <- seq(0.01, 2*pi, len=100)
    theta <- c(.00001, 1.0)
    Q     <- hd.precision(theta=theta, order=3)
    logpar    <- theta
    for(i in 1:100){
        ## logpar    <- c(log(200), log(1))
        ## Q         <- temp.precision(theta=logpar, mesh=mesh1d, o=2)
        x         <- inla.qsample(Q=Q)
        tmp2[i,]  <- x
        A         <- circular.make.A(hd.vec=y, order=3)
        sim[i,]   <- as.numeric(A %*% x)
        tmp[i]    <- var(sim[i,])
    }
    matplot(y, t(sim[1:10,]), type="l", col=1)

    plot(sim[2,])
    abline(h=sim[2,1])
}

if(FALSE){
    tmp  <- NULL
    sim  <- matrix(nrow=1000, ncol=1000)
    tmp2 <- matrix(nrow=1000, ncol=1000)
    y <- seq(0.01, 1000, len=100)
    theta.Z <- theta[4:5]    
    Q     <- theta.precision(theta=theta.Z, mesh=mesh1d)
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
