objective_osc_temp <- function(par, theta, data, mesh, mesh1d, A, Atilde, Aobs, Atildeobs, W) {
    ## input:
    ##        par: list containing vector X and scalar beta
    ##        theta: fixed hyperprior parameters!! (revisit selection)
    ##        spde: a model object for a Matern Gaussian model created with inla.spde2
    ##              using a PC prior for the parameters.        
    ##        data: list of Y and Ypos (list(Y, Ypos))
    ##              - Y: data matrix of firing events
    ##              - Ypos: data matrix of positional data
    ##        mesh: Delaunay triangulation of spatial domain
    ##              constructed using inla.mesh.2d from INLA package
    ##        L: vector of areas of line segments used for 
    ##           the approximation of the integral of the
    ##           void probability. 
    ##        A: matrix of basis functions (psi_k(s_i))_{ik} evaluated        
    ##          at midpoints of line segments overline{s}_i
    ##        Aobs: matrix of basis functions (psi_k(s_i))_{ik}
    ##              evaluated at observed points s_i
    ## output: posterior log density of X,beta | mid theta, Y. This
    ## function is optimized in a nested optimization procedure to give the
    ## mode of the posterior of X and beta at a
    ## configuration of the parameter theta (range and scale of
    ## GMRF Matern).
    beta <- par$beta
    X    <- par$X
    Z    <- par$Z
    out  <-  - prior.betaXZ_osc_temp(beta=beta, X=X, Z=Z, theta=theta, mesh = mesh, mesh1d=mesh1d) -
        pp.llik(data=data, X=X, Z=Z, mesh=mesh, mesh1d=mesh1d, beta=beta,
                A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs, W=W)
    return(out)
}
