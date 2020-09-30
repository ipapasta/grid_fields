prior.beta_osc_hd <- function(beta, X, theta, mesh, order.hd)
{
    ## ---------------------------------------
    ## input:
    ##       X: vector of parameters X
    ##       beta: vector of parameters X
    ##       spde: model object for oscillating Gaussian field created with inla.spde2
    ##       theta: hyperprior parameters
    ##       mesh: Delaunay triangulation of spatial domain constructed using inla.mesh.2d from INLA package
    ##             
    ## output:
    ##        prior log density of X, beta given theta
    ## ---------------------------------------
    theta.osc <- theta[1:3]
    theta.hd  <- theta[4:5]    
    Q.osc     <- osc.precision(theta=theta.osc, mesh=mesh)
    Q.hd      <- hd.precision(theta=theta.hd, order=order.hd)
    Q         <- Q.hd %x% Q.osc
    lpbeta <- dnorm(as.numeric(beta), mean=0, sd=1/sqrt(1e-6), log=TRUE)
    lpX    <- ldmvnorm(X, matrix(rep(0, length(X)), ncol=1), Q)
    out    <- (lpbeta + lpX )
    return(out)
}



