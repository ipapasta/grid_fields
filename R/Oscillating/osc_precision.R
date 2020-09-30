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

osc.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    theta3 <- theta[3]
    rho     <- exp(theta1)
    kappa   <- sqrt(8)/rho
    sigma   <- exp(theta2)
    phi     <- (1-exp(-theta3))/(1+exp(-theta3))
    sincpth <- sqrt(1-phi^2)/acos(phi)
    tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
    ## ---------------------------------------
    o       <- inla.mesh.fem(mesh, order = 2)
    ## 
    tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*phi*(o$g1)) + o$g2)
    ## 
}

## ---------------------------------------------------
## experimentation of oscillating covariance structure
## ---------------------------------------------------


if(FALSE){
    ## for(i in 1:100){
    
    theta1  <- 2.51                           # rho/kappa
    theta2  <- 2.                           # sigma
    theta3  <- -8                         # phi
    theta1  <- 2.12                           # rho/kappa
    theta2  <- 4                            # sigma
    theta3  <- -8.61                         # phi
    rho     <- exp(theta1)
    kappa   <- sqrt(8)*exp(-theta1)
    sigma   <- exp(theta2)
    phi     <- (1-exp(-theta3))/(1+exp(-theta3))
    print(paste("rho is: ", rho, " kappa is: ", kappa, " sigma is: ", sigma, " phi is: ", phi))
    tausq   <- 1/(4*pi*(sigma^2)*kappa^2)
    QO      <- osc.precision(c(theta1, theta2, theta3), mesh)    ## for some theta1, theta2, theta3 values and my mesh
    o       <- inla.mesh.fem(mesh, order = 2)
    QM      <- tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*(o$g1)) + o$g2)
    x       <- inla.qsample(Q=QO)
    xM      <- inla.qsample(Q=QM) 
    locs    <- as.matrix(expand.grid(0:100, 0:100))
    Apred   <- inla.spde.make.A(mesh=mesh, loc=locs)
    sim     <- Apred%*%x
    simM    <- Apred%*%xM
    
    df      <- data.frame(x=locs[,1], y=locs[,2], data=as.numeric(sim))
    dfM     <- data.frame(x=locs[,1], y=locs[,2], data=as.numeric(simM))
    ## 
    ##
    library(colorRamps)
    p <- ggplot(data=df) + geom_raster(aes(x=x,y=y, fill=data), interpolate=FALSE) + coord_fixed()+
        scale_fill_viridis_c(option="Viridis")+
        ## scale_color_gradientn(colours=matlab.like(10)) +
        theme_classic() + theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20),
                                text = element_text(size=24)) +
        xlab("x") + ylab("y") + ggtitle("Oscillating")
    ## 
    pM <- ggplot(data=dfM) + geom_raster(aes(x=x,y=y, fill=data), interpolate=FALSE) + coord_fixed()+
        scale_fill_viridis_c(option="Viridis")+
        ## scale_fill_brewer()+        
        theme_classic() + theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20),
                                text = element_text(size=24)) +
        xlab("x") + ylab("y") + ggtitle("Matern")
    plt <- grid.arrange(p, pM, nrow=1)
    plt
    
    a <- ggsave(plt, filename="/home/ipapasta/Desktop/sims.pdf")
    
}


