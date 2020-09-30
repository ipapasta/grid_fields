grad.objective <- function(par, theta, spde,  data,  mesh, L, A, Aobs)
    ## input:
    ##        par: named list containing vector of parameters X and scalar parameter beta
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
    ## output:
    ##         gradient vector of posterior log density of X,beta | mid theta, Y. 
{
    n    <- nrow(data$Y)
    X    <- Matrix(par$X, ncol=1)
    p    <- nrow(X)
    beta <- Matrix(par$beta, ncol=1)
    tL   <- t(L)
    Q.X    <- inla.spde.precision(spde, theta=log(theta)) #!
    Q.beta <- Diagonal(length(beta), 1e-5) 
    expAx  <- exp(A%*%X)
    one.n     <- Matrix(rep(1, n), ncol=1)
    one.p     <- Matrix(rep(1, p), ncol=1)
    ## browser()
    ## ----------------------------------------------------
    grad.X     <- -as.numeric(exp(beta))*(t( (L%*%t(one.p))*A) %*% expAx) + (t(Aobs) %*% one.n)
    ## ----------------------------------------------------
    grad.beta  <- -as.numeric(exp(beta))*(tL%*%(expAx)) + n
    ## ----------------------------------------------------  
    grad.Xbeta <- -rBind(grad.X, grad.beta) + rBind(Q.X%*%X, Q.beta%*%beta)
    ##
    return(grad.Xbeta)
}        






## gradient of X alone
gradX.objective <- function(par, theta, spde,  data,  mesh, L, A, Aobs)
    ## input:
    ##        par: named list containing vector of parameters X and scalar parameter beta
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
    ## output:
    ##         gradient vector of posterior log density of X,beta | mid theta, Y. 
{
    n    <- nrow(data$Y)
    X    <- Matrix(par$X, ncol=1)
    p    <- nrow(X)
    beta <- Matrix(par$beta, ncol=1)
    tL   <- t(L)
    Q.X    <- inla.spde.precision(spde, theta=log(theta)) #!
    ## Q.beta <- Diagonal(length(beta), 1e-5) 
    expAx  <- exp(A%*%X)
    one.n     <- Matrix(rep(1, n), ncol=1)
    one.p     <- Matrix(rep(1, p), ncol=1)
    ## browser()
    ## ----------------------------------------------------
    grad.X     <- -as.numeric(exp(beta))*(t( (L%*%t(one.p))*A) %*% expAx) + (t(Aobs) %*% one.n)
    ## ----------------------------------------------------
    ## grad.beta  <- -as.numeric(exp(beta))*(tL%*%(expAx)) + n
    ## ## ----------------------------------------------------  
    ## grad.Xbeta <- -rBind(grad.X, grad.beta) + rBind(Q.X%*%X, Q.beta%*%beta)
    ##
    return(grad.X)
}        
