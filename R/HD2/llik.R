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
## A    <- inla.spde.make.A(mesh=mesh, loc=as.matrix(do.call("rbind",Ypos$midpoints)))
## W    <- Matrix(t(t(L)%*%A))

pp.llik <- function(data, X, mesh, beta, A, Aobs, W, tW, method="biased") # B, Bobs
{
    Y            <- data$Y
    n            <- nrow(Y)
    p            <- nrow(X)
    one.n        <- matrix(rep(1, n), ncol=1)
    beta         <- matrix(beta, ncol=1)
    X            <- matrix(X, ncol=1)
    expAx        <- exp(A%*%X)
    sumlogR      <- sum(one.n * (Aobs %*% X)) + n * beta ## tr.one.n %*% Bobs %*% beta        
    integr       <- sum(W*expAx)*as.numeric(exp(beta))
    out          <- -integr + sumlogR 
    print(paste("sumlogR is: ", sumlogR, "integr is: ", integr, "llik is -I+S: ", out))
    return(out)
    ## }
} 



if(FALSE){
    pp.llik <- function(data, X, mesh, beta, A, Aobs, W, tW, method="biased") # B, Bobs
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
        ## A    <- inla.spde.make.A(mesh=mesh, loc=as.matrix(do.call("rbind",Ypos$midpoints)))
        ## W    <- Matrix(t(t(L)%*%A))
        Y            <- data$Y
        n            <- nrow(Y)
        p            <- nrow(X)
        one.n        <- matrix(rep(1, n), ncol=1)
        beta         <- matrix(beta, ncol=1)
        X            <- matrix(X, ncol=1)
        sumlogR      <- sum(one.n * (Aobs %*% X)) + n * beta ## tr.one.n %*% Bobs %*% beta        
        ## integr       <- ((tW%*%(exp(X))))*as.numeric(exp(beta))
        integr       <- sum(W*exp(X))*as.numeric(exp(beta))
        out          <- -integr + sumlogR 
        print(paste("sumlogR is: ", sumlogR, "integr is: ", integr, "llik is -I+S: ", out))
        return(out)
        ## }
    } 
}


