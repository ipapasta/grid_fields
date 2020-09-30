prior.Xbeta_osc <- function(X, beta, theta, mesh)
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
    QX     <- osc.precision(theta=theta, mesh=mesh) 
    lpbeta <- dnorm(as.numeric(beta), mean=0, sd=1/sqrt(1e-6), log=TRUE)
    lpX    <- ldmvnorm(X, matrix(rep(0, length(X)), ncol=1), QX)
    out    <- (lpX + lpbeta)
    return(out)
}



