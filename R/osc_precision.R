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



osoM.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    theta3 <- theta[3]
    rho     <- exp(theta1)
    kappa   <- sqrt(8)/rho
    sigma   <- exp(theta2)
    phi     <- (1+exp(-theta3))/(exp(-theta3))
    ## sincpth <- sqrt(1-phi^2)/acos(phi)
    tausq   <- log(-1+2*phi*(phi+sqrt((phi^2)-1)))/(8*pi*(sigma^2)*(kappa^2)*sqrt((phi^2) -1))
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
    theta1 <- c(3, -1/2, -5 )
    theta2 <- c(3, -1/2, -1 )
    theta3 <- c(3, -1/2, 1)
    theta4 <- c(3, -1/2, 10)
    Q1       <- osoM.precision(theta1, mesh)
    Q1.osc   <- osc.precision(theta1, mesh)
    Q2   <- osoM.precision(theta2, mesh)
    Q3   <- osoM.precision(theta3, mesh)
    Q4   <- osoM.precision(theta4, mesh)
    x1  <- inla.qsample(n=1, Q=Q1)
    x2  <- inla.qsample(n=1, Q=Q2)
    x3  <- inla.qsample(n=1, Q=Q3)
    x4  <- inla.qsample(n=1, Q=Q4)
    ## 
    coords  <- expand.grid(seq(-20, 120, seq=100), seq(-20, 120, seq=100))
    proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
    Aosc.pred      <- proj.s.pred$proj$A
    field <- cbind(as.numeric(Aosc.pred %*% x1), as.numeric(Aosc.pred %*% x2), 
                   as.numeric(Aosc.pred %*% x3), as.numeric(Aosc.pred %*% x4))
    df.plot       <- data.frame(x=coords[,1], y=coords[,2], field1=field[,1], field2=field[,2], field3=field[,3], field4=field[,4])
    lims <- c(min(c(df.plot$field1, df.plot$field2, df.plot$field3, df.plot$field4)),max(c(df.plot$field1, df.plot$field2, df.plot$field3, df.plot$field4)))
    p1  <- ggplot(data = df.plot) + geom_raster(aes(x,y,fill=field1), interpolate=TRUE) + 
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                             limits=lims) + coord_fixed() + theme_classic()
    p2  <- ggplot(data = df.plot) + geom_raster(aes(x,y,fill=field2), interpolate=TRUE) + 
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                             limits=lims) + coord_fixed() + theme_classic()
    p3  <- ggplot(data = df.plot) + geom_raster(aes(x,y,fill=field3), interpolate=TRUE) + 
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                             limits=lims) + coord_fixed() + theme_classic()
    p4  <- ggplot(data = df.plot) + geom_raster(aes(x,y,fill=field4), interpolate=TRUE) + 
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                             limits=lims) + coord_fixed() + theme_classic()
    a <- grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)

}


if(FALSE){
    ## for(i in 1:100){
    
    theta1  <- 2.51                           # rho/kappa
    theta2  <- 2.                             # sigma
    theta3  <- -8                             # phi
    theta1  <- 2.12                           # rho/kappa
    theta2  <- 4                              # sigma
    theta3  <- -8.61                          # phi
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


