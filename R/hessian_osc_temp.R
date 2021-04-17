hessian.objective_osc_temp<- function(par, theta,  data,  mesh, mesh.theta, mesh1d, A, Atilde, Aobs, Atildeobs, W){
    ## , W.input.to.sparseMatrix
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
    ##         hessian matrix of posterior log density of (X, beta) | (theta, Y)
    X         <- Matrix(par$X, ncol=1)
    Z         <- Matrix(par$Z, ncol=1)
    beta      <- Matrix(par$beta, ncol=1)
    theta.X     <- theta[1:3]
    theta.Theta <- theta[4:5]
    theta.Z <- theta[6:7]
    ## --------------
    n         <- nrow(data$Y)
    p         <- nrow(X)
    ptilde    <- nrow(Z)
    ## --------------
    tW      <- t(W)
    expbeta <- exp(beta)
    expX      <- Matrix(exp(X))
    expZ      <- Matrix(exp(Z))
    Q.beta  <- Diagonal(length(beta), 1e-6)
    Q.Omega     <- osc.precision(theta=theta.X, mesh=mesh)
    Q.Theta <- hd.bsp.precision(theta=theta.Theta, mesh=mesh.theta)
    Q.X     <- kronecker(Q.Theta, Q.Omega)
    Q.Z     <- temp.precision(theta=theta.Z, mesh=mesh1d)
    ## --------------------------------------------------
    gradbeta.gradbeta_T  <- -as.numeric(exp(beta))*(t(expZ) %*% (W%*%expX)) 
    ##------
    bX.T <- (t(tW %*% expZ) * (t(expX)))
    bZ.T <- (t(W %*% expX) * (t(expZ)))
    gradbeta.gradX_T  <- -as.numeric(expbeta) * bX.T
    gradbeta.gradZ_T  <- -as.numeric(expbeta) * bZ.T
    gradX.gradX_T  <-  Matrix(-as.numeric(expbeta) * Diagonal(p, as.numeric(bX.T)),sparse=TRUE)
    gradZ.gradZ_T  <-  Matrix(-as.numeric(expbeta) * Diagonal(ptilde, as.numeric(bZ.T)),sparse=TRUE)
    Wtmp <- as(W, "dgTMatrix")
    Wtmp <- cbind(Wtmp@i+1, Wtmp@j+1, Wtmp@x)
    ## 
    gradZ.gradX_T  <-  -as.numeric(expbeta) * sparseMatrix(i=Wtmp[,1], j=Wtmp[,2], x=(Wtmp[,3] * exp(Z[Wtmp[,1]]) * exp(X[Wtmp[,2]])))
    gradZ.gradX_T  <- cbind(gradZ.gradX_T, Matrix(0, nrow=nrow(gradZ.gradX_T), ncol=ncol(W)-ncol(gradZ.gradX_T)))
    ##
    out <- Matrix(rBind(cBind(-gradbeta.gradbeta_T + Q.beta, -gradbeta.gradX_T, -gradbeta.gradZ_T),
                        cBind(-t(gradbeta.gradX_T), -gradX.gradX_T + Q.X, -t(gradZ.gradX_T)),
                        cBind(-t(gradbeta.gradZ_T), -gradZ.gradX_T, -gradZ.gradZ_T + Q.Z)),sparse=TRUE)
    if(FALSE){
        ## 
        prior <- bdiag(Q.beta, Q.X, Q.Z)
        ## 
        likelihood <- Matrix(rBind(cBind(-gradbeta.gradbeta_T , -gradbeta.gradX_T, -gradbeta.gradZ_T),
                        cBind(-t(gradbeta.gradX_T), -gradX.gradX_T , -t(gradZ.gradX_T)),
                        cBind(-t(gradbeta.gradZ_T), -gradZ.gradX_T, -gradZ.gradZ_T )),sparse=TRUE)
    }
    return(out)
}


if(FALSE){
    hessianX.objective <- function(par, theta, spde, data, mesh, L, A, Aobs)
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
        ##         hessian matrix of posterior log density of X,beta | mid theta, Y. 
    {
        X         <- Matrix(par$X, ncol=1)
        n         <- nrow(data$Y)
        p         <- nrow(X)
        one.n     <- Matrix(rep(1, n), ncol=1)
        one.p     <- Matrix(rep(1, p), ncol=1)        
        expAx     <- exp(A%*%X)
        beta      <- Matrix(par$beta, ncol=1)
        Q.X       <- inla.spde.precision(spde, theta=log(theta))
        Q.beta    <- Diagonal(length(beta), 1e-6) 
        hess.X    <- Matrix(0, nrow=p, ncol=p)
        ## QU: is it possible to compute the hessian faster?
        ## options:
        ## 1) matrix algebra?
        ## 2) [X] C++ code: roughly 9 times faster than R code
        ## for(i in 1:(p-1))
        ##     for(j in i:p){
        ##         hess.X[i,j] <- -as.numeric(exp(beta))*(t(A[,i]*A[,j]*L)) %*% expAx
        ##     }   
        ## hess.X   <- Matrix::forceSymmetric(hess.X, uplo="U")
        hess.X <- hessian_diag_and_upper_triangle(beta=as.numeric(beta), p,
                                                  L=as.matrix(L), A=as.matrix(A),
                                                  expAx=as.matrix(expAx))
        hess.X <- Matrix::forceSymmetric(hess.X, uplo="U")
        ## in Rcpp dense matrices of type dgeMatrix from Matrix package
        ## are not supported yet hess.X <- hess.X + t(hess.X) -
        ## hess.beta  <- Matrix(0,nrow=length(beta),
        ##                      ncol=length(beta))
        ## hess.Xbeta <- Matrix(0,nrow=length(X),
        ##                      ncol=length(beta))
        ## --------------------------------------------------
        ## hess.beta  <- Matrix(-as.numeric(exp(beta))*(t(L)%*%expAx))
        ## ## 
        ## hess.Xbeta <- Matrix(-as.numeric(exp(beta))*(t((L%*%t(one.p))*A) %*% expAx))
        ## 
        ## 
        ## out <- rBind(cBind(-hess.X + Q.X, -hess.Xbeta),
        ##              cBind(-t(hess.Xbeta), -hess.beta + Q.beta))
        out <- -hess.X + Q.X
        return(out)
    }
}
