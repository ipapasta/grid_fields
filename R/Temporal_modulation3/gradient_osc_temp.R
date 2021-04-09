grad.objective_osc_temp <- function(par, theta,  data,  mesh, mesh.theta, mesh1d, A, Atilde, Aobs, Atildeobs, W){
    ## input:
    ##        par: named list containing vector of parameters X and scalar parameter beta
    ##        theta: fixed hyperprior parameters!! (revisit selection)
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
    ## output:
    ##         gradient vector of posterior log density of (X, beta) | (theta, Y)
    n    <- nrow(data$Y)
    X    <- Matrix(par$X, ncol=1)
    Z    <- Matrix(par$Z, ncol=1)
    beta <- Matrix(par$beta, ncol=1)
    theta.X     <- theta[1:3]
    theta.Theta <- theta[4:5]
    theta.Z <- theta[6:7]
    ## -----------------------------
    one.n   <- Matrix(rep(1, n), ncol=1)
    tW      <- t(W)
    expbeta <- exp(beta)
    expX    <- exp(X)
    expZ    <- exp(Z)
    tAobs   <- t(Aobs)
    tAtildeobs  <- t(Atildeobs)
    Atildeobs <- t(Atildeobs)
    ## -----------------------------
    Q.Omega     <- osc.precision(theta=theta.X, mesh=mesh)
    Q.Theta <- hd.bsp.precision(theta=theta.Theta, mesh=mesh.theta)
    Q.X     <- kronecker(Q.Theta, Q.Omega)
    Q.Z     <- temp.precision(theta=theta.Z, mesh=mesh1d) 
    Q.beta  <- Diagonal(length(beta), 1e-6)
    ## -----------------------------
    grad.beta  <- - as.numeric(expbeta) * (t(expZ)%*%(W%*%expX)) + n
    grad.X     <- - as.numeric(expbeta) * ((tW%*%expZ) * expX) + tAobs %*% one.n
    grad.Z     <- - as.numeric(expbeta) * ((W %*%expX) * expZ) + tAtildeobs %*% one.n
    grad.Xbeta <- Matrix(-rBind(grad.beta, grad.X, grad.Z) + rBind(Q.beta%*%beta, Q.X%*%X, Q.Z%*%Z ),ncol=1) # grad beta, X, Z
    return(grad.Xbeta)    
}        



