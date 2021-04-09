pp.llik <- function(data, X, Z, mesh, mesh1d, beta, A, Atilde, Aobs, Atildeobs, W) 
{
    ## ----------------------------------------------------------------
    ## input: data: list of Y and Ypos (list(Y, Ypos))
    ##        Y: data matrix of firing events
    ##        Ypos: data matrix of positional data
    ##        mesh: Delaunay triangulation of spatial domain
    ##              constructed using inla.mesh.2d from INLA package
    ##        L: vector of areas of line segments used for 
    ##           the approximation of the integral of the
    ##           void probability. 
    ##        beta: intercept of linear predictor
    ##        A: matrix of basis functions (psi_k(s_i))_{ik} evaluated        
    ##          at midpoints of line segments overline{s}_i
    ##        Aobs: matrix of basis functions (psi_k(s_i))_{ik}
    ##              evaluated at observed points s_i
    ## output: log-likelihood from point process conditioned on X.
    ## ----------------------------------------------------------------
    Y            <- data$Y
    n            <- nrow(Y)
    ## p            <- nrow(X)
    ## ptilde       <- nrow(Z)
    tr.one.n     <- matrix(rep(1, n), nrow=1)
    beta         <- matrix(beta, ncol=1)
    X            <- matrix(X, ncol=1)
    expX         <- exp(X)
    expZ         <- exp(Z)
    tW           <- t(W)
    sumlogR      <- tr.one.n %*% ((Aobs %*% X) + (Atildeobs %*% Z)) + n * beta 
    ##
    integr       <- as.numeric(exp(beta)) * (t(expZ)%*%(W%*%expX))
    out          <- sumlogR - integr
    out
}
