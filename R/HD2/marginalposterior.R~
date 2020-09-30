## par.theta <- 

## par.theta <- c(log(19.9577192834812), log(1.04979866692695), -log((1-(-0.977523485236457))/(1+(-0.977523485236457))),
##          log(1.20177352113122), log(0.n0861727221300692))


## par.theta <- 

## par.theta <- c(log(19.9577192834812), log(1.04979866692695), -log((1-(-0.977523485236457))/(1+(-0.977523485236457))),
##          log(1.20177352113122), log(0.n0861727221300692))

pthetapc.prop.marg.post_osc_HD <- function(par.theta, hyperpar, data, X, beta, mesh, A, tA, W, tW=tW, DW=DW, Aobs, tAobs, order.HD, type="biased", acc=1e-4, print.verbose=FALSE, return.X=FALSE)
{
    ## input:
    ##        par.theta: add description
    ##        spde: a model object for a Matern Gaussian model created with inla.spde2
    ##              using a PC prior for the parameters.        
    ##        data: list of Y and Ypos (list(Y, Ypos))
    ##              - Y: data matrix of firing events
    ##              - Ypos: data matrix of positional data
    ##        X: vector of parameters X
    ##        beta: scalar parameter beta
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
    ##        Approximate posterior log density of theta given Y
    ##        up to proportionality. This function is 
    ##        optimized with the polytope method (Nelder-Mead)
    method      <- type
    par.theta   <- par.theta
    par         <- list(X=X, beta=beta)
    Xbeta       <- rBind(par$X, par$beta) ## starting value
    n.X         <- length(X)
    ## n.beta      <- beta %>% length
    ## 
    ## ---------------------------------------
    ## Newton optimization
    ## ---------------------------------------
    ## starting values
    print(paste("rho is:", exp(par.theta[1]), " sigma is:", exp(par.theta[2]),  "phi is: ", (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])), "rhoHD is:", exp(par.theta[4]), "sigmaHD is:", exp(par.theta[5])))

    ## Code to replace original calls to grad.objective_osc_HD and hessian.objective_osc_HD  (Mike Sept 2019)
    n    <- nrow(data$Y)
    X    <- Matrix(par$X, ncol=1)
    p    <- nrow(X)
    beta <- Matrix(par$beta, ncol=1)
    one.n    <- Matrix(rep(1, n), ncol=1)
    one.p    <- Matrix(rep(1, p), ncol=1)
    expAx     <- exp(A%*%X)
    expbeta  <- exp(beta)
    Q.Omega  <- osc.precision(theta=par.theta[1:3], mesh=mesh)
    Q.HD     <- hd.precision(theta=par.theta[4:5], order=order.HD)
    Q.X      <- kronecker(Q.HD, Q.Omega)
    Q.beta   <- Diagonal(length(beta), 1e-6)
    grad.X     <- -as.numeric(expbeta)*tA%*%(DW%*%expAx) + (tAobs%*%one.n)
    grad.beta  <- -as.numeric(expbeta)*sum(W*expAx) + n
    grad <- matrix(-rBind(grad.X, grad.beta) + rBind(Q.X%*%X, Q.beta%*%beta), ncol=1)
    hess.beta  <- Matrix(0, nrow=length(beta), ncol=length(beta))
    hess.Xbeta <- Matrix(0, nrow=length(X), ncol=length(beta))
    hess.X     <- -as.numeric(expbeta)* (tA %*% (Diagonal(length(W), as.numeric(W*expAx)) %*% A))
    hess.beta  <- -as.numeric(expbeta)*sum(W*expAx)
    hess.Xbeta <- Matrix( -as.numeric(expbeta)*tA%*%(Diagonal(length(W), W)%*%expAx))
    hess <- Matrix(Matrix::forceSymmetric(rBind( cBind(-hess.X + Q.X, -hess.Xbeta), cBind(-t(hess.Xbeta), -hess.beta + Q.beta) )), sparse=TRUE)

    fval.p   <- fval.c <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh,
                                                      A=A, Aobs=Aobs, W=W, tW=tW, order.HD = order.HD, method=method))  
    scaled.acc  <- (max( abs(fval.c), 1)) * acc  ## (fval.c %>% abs %>% max(1)) * acc     
    tol.test    <- Inf 
    ##
    if(print.verbose){    
        print(paste("tol.test is :", tol.test, "scaled.acc is:", scaled.acc, "length grad is: ", sqrt(sum(grad^2))))
    }
    ## fdiff <- TRUE
    while(TRUE){    #!!!{
        Xbeta.p <- Xbeta
        alpha   <- 1
        fval.p  <- fval.c            
        par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1): (length(Xbeta.p))])
        ## 
        ## Newton move
        ## 
        ## ------------------------------------------------------------------------------------------------------
        ## Code to replace original calls to grad.objective_osc_HD and hessian.objective_osc_HD  (Mike Sept 2019)
        X    <- Matrix(par$X, ncol=1)
        p <- nrow(X)
        one.p      <- Matrix(rep(1, p), ncol=1)
        expAx      <- exp(A%*%X)
        beta       <- Matrix(par$beta, ncol=1)
        expbeta    <- exp(beta)
        Q.beta     <- Diagonal(length(beta), 1e-6)
        grad.X     <- -as.numeric(expbeta)*tA%*%(DW%*%expAx) + (tAobs%*%one.n)
        grad.beta  <- -as.numeric(expbeta)*sum(W*expAx) + n
        grad       <- (-rBind(grad.X, grad.beta) + rBind(Q.X%*%X, Q.beta%*%beta)) %>% Matrix(ncol=1)
        hess.beta  <- Matrix(0, nrow=length(beta), ncol=length(beta))
        hess.Xbeta <- Matrix(0, nrow=length(X), ncol=length(beta))
        hess.X     <- -as.numeric(expbeta)* (tA %*% (Diagonal(length(W), as.numeric(W*expAx)) %*% A))
        hess.beta  <- -as.numeric(expbeta)*sum(W*expAx)
        hess.Xbeta <- Matrix( -as.numeric(expbeta)*tA%*%(Diagonal(length(W), W)%*%expAx))
        hess <- Matrix(Matrix::forceSymmetric(rBind( cBind(-hess.X + Q.X, -hess.Xbeta), cBind(-t(hess.Xbeta), -hess.beta + Q.beta) )), sparse=TRUE)
        ##
        d       <- -Matrix::solve(hess, grad)
        ## d2      <- -grad  
        Xbeta   <- Xbeta.p + alpha * d            
        par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1): (length(Xbeta))])
        ## print(par$X)
        ## 
        fval.c  <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh, A=A, Aobs=Aobs, W=W, tW=tW,
                                               order.HD = order.HD, method=method))
        ## uncomment to track progress of Newton algo
        if(print.verbose){    
            print(paste("fval.c is larger than fval.p", (fval.c > fval.p)))
        }
        ##!! isSymmetric(hess)          
        while(fval.c > fval.p){                
            alpha <- alpha/2
            ## d <- grad
            Xbeta   <- Xbeta.p + alpha * d
            par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1):(length(Xbeta))])
            ## 
            fval.c  <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh, A=A,
                                                   Aobs=Aobs, W=W, tW=tW, order.HD = order.HD, method=method))
            ## 
            tmp <- (fval.c - fval.p)
            ## uncomment to track progress of Newton algo
            if(print.verbose){    
                print(paste("alpha is:", alpha,", ", "Fc-Fp is:", round(tmp,2),  "fval.c is:", fval.c))
            }
            if(fval.c < fval.p) print(paste("REDUCED STEP SIZE REACHED Fc < Fp"))
            if(alpha <= 2^(-9))
            {
                Xbeta <- Xbeta.p
                par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1):(length(Xbeta.p))])
                ## 
                fval.c  <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh, A=A, Aobs=Aobs, W=W, tW=tW,
                                                       order.HD = order.HD, method=method))
                ## uncomment to track progress of Newton algo
                if(print.verbose){    
                    print(paste("STOP MOVING"))
                }
                break
            }
        }
        scaled.acc <- (max(abs(fval.c),1)) * acc
        scaled.acc <- acc
        tol.test   <- sqrt(sum(grad^2))
        if(print.verbose){    
            print(paste("Fc is:", fval.c, "Fp is:", round(fval.p,3), "Tol test is:", round(tol.test,3), "Scaled accuracy is:", round(scaled.acc,3), "length grad is: ", sqrt(sum(grad^2))))
        }
        nmv <- (sqrt(sum((Xbeta-Xbeta.p)^2)))
        if(tol.test < scaled.acc | nmv < acc) {
            print(paste("Tol test < scaled accuracy", tol.test < scaled.acc))
            print(paste("No move < accuracy", nmv < acc))
            break
        }                    
    }
    if(print.verbose){    
        print("EXIT while loop")
    }
    lp.Xbeta <- prior.betaX_osc_HD(X=par$X, beta=par$beta, theta=c(par.theta),  mesh=mesh, order.HD = order.HD)
    Lxbeta.given.thetay    <- Matrix::expand(Matrix::Cholesky(hess))$L
    logdetQxbeta.given.thetay <- 2*sum(log(diag(Lxbeta.given.thetay)))
    lp.G     <- .5*logdetQxbeta.given.thetay
    llik     <- pp.llik(data=data, X=par$X, mesh=mesh, beta=par$beta, A=A,  Aobs=Aobs, W=W, tW=tW, method=method)
    ## output
    Xest    <<- par$X
    betaest <<- par$beta
    gradest     <<- grad
    Hessianest  <<- hess
    lp.theta    <- prior.theta_osc_HD(rho = exp(par.theta[1]),
                                      sigma = exp(par.theta[2]),
                                      ## phi = (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
                                      phi   =(1+0.99)*((exp(par.theta[3]))/(1+exp(par.theta[3]))) - 0.99, # -0.99 lower bound on phi to avoid computational issues.
                                      rhoL = hyperpar$rhoL,
                                      sigmaU=hyperpar$sigmaU,
                                      rho.HD = exp(par.theta[4]),
                                      sigma.HD = exp(par.theta[5]),
                                      alpha1 = hyperpar$alpha1,
                                      alpha2 = hyperpar$alpha2,
                                      alphaHD1 = hyperpar$alpha1,
                                      alphaHD2 =hyperpar$alpha2,
                                      rho.HDL=hyperpar$rhoL,
                                      sigma.HDU=hyperpar$sigmaU,
                                      lg=TRUE)
    ## lp.theta <- dexp(exp(par.theta[2]), rate=1/0.33, log=TRUE)
    obj <- as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G))
    print(paste("+prior.theta: ", round(lp.theta, 3), "+prior.Xbeta: ", round(lp.Xbeta, 3),
                "+llik: ", round(llik, 3), "-Xbeta.y: ", round(-lp.G, 3), "-val: ",
                round(as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G)), 3)))
    obj <- as.numeric(-(lp.Xbeta + llik - lp.G))    
    ## save(betaest, Xest, par.theta, gradest, Hessianest, obj,
    ## file = paste0("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating/data/Xbeta", counter, ".RData"))
    if(TRUE){
        save(betaest, Xest, par.theta, gradest, Hessianest, obj, file = paste0("/data/ipapasta/biased/Xbeta", counter, ".RData"))
    }
    cat("\n \n \n")
    counter <<- counter+1
    cat("\n \n \n ")               
    print(paste("counter is equal to:", counter))
    if(!return.X) return(as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G)))
    if(return.X){
        o<-list(Xest=Xest, betaest=betaest, Hessianest=Hessianest)
        return(o)
    }
}



if(FALSE){
    pthetapc.prop.marg.post_osc_HD.2 <- function(par.theta, hyperpar, data, X, beta, mesh, L, A, tA, W, tW=tW, DW=DW, Aobs, tAobs, order.HD, type="biased", acc=1e-7, print.verbose=FALSE)
    {
        ## input:
        ##        par.theta: add description
        ##        spde: a model object for a Matern Gaussian model created with inla.spde2
        ##              using a PC prior for the parameters.        
        ##        data: list of Y and Ypos (list(Y, Ypos))
        ##              - Y: data matrix of firing events
        ##              - Ypos: data matrix of positional data
        ##        X: vector of parameters X
        ##        beta: scalar parameter beta
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
        ##        Approximate posterior log density of theta given Y
        ##        up to proportionality. This function is 
        ##        optimized with the polytope method (Nelder-Mead)
        method      <- type
        par.theta   <- par.theta
        par         <- list(X=X, beta=beta)
        Xbeta       <- rBind(par$X, par$beta) ## starting value
        n.X         <- length(X)
        ## n.beta      <- beta %>% length
        ## 
        ## ---------------------------------------
        ## Newton optimization
        ## ---------------------------------------
        ## starting values
        print(paste("rho is:", exp(par.theta[1]), " sigma is:", exp(par.theta[2]),
                    "phi is: ", (1-exp(-par.theta[3]))/(1+exp(-par.theta[3]))))
        print(paste("rhoHD is:", exp(par.theta[4]), "sigmaHD is:", exp(par.theta[5])))

        grad     <- matrix(grad.objective_osc_HD(par, theta=par.theta, data=data, mesh=mesh, A=A, tA=tA, W=W, tW=tW, DW=DW, tAobs=tAobs, order.HD = order.HD, method=method), ncol=1)
        hess     <- Matrix(hessian.objective_osc_HD(par, theta=par.theta, data=data, mesh=mesh, A=A, tA=tA, W=W, tW=tW, Aobs=Aobs, order.HD = order.HD, method=method), sparse=TRUE)
        fval.p   <- fval.c <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh,
                                                          A=A, Aobs=Aobs, W=W, tW=tW, order.HD = order.HD, method=method))  
        scaled.acc  <- (max( abs(fval.c), 1)) * acc  ## (fval.c %>% abs %>% max(1)) * acc     
        tol.test    <- Inf 
        ##
        if(print.verbose){    
            print(paste("tol.test is :", tol.test,
                        "scaled.acc is:", scaled.acc))
        }
        ## fdiff <- TRUE
        while(TRUE){    #!!!{
            Xbeta.p <- Xbeta
            alpha   <- 1
            fval.p  <- fval.c            
            par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1): (length(Xbeta.p))])
            ## 
            ## Newton move
            ## 
            grad  <- grad.objective_osc_HD(par, theta=par.theta, data=data, mesh=mesh, A=A, tA=tA, W=W, tW=tW, DW=DW, tAobs=tAobs, order.HD = order.HD, method=method) %>% Matrix(ncol=1)
            hess  <- Matrix(hessian.objective_osc_HD(par, theta=par.theta, data=data, mesh=mesh, A=A, tA=tA, W=W, tW=tW, Aobs=Aobs, order.HD = order.HD, method=method), sparse=TRUE)
            ## 
            d       <- -Matrix::solve(hess, grad)
            ## d2      <- -grad  
            Xbeta   <- Xbeta.p + alpha * d            
            par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1): (length(Xbeta))])
            ## print(par$X)
            ## 
            fval.c  <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh, A=A, Aobs=Aobs, W=W, tW=tW,
                                                   order.HD = order.HD, method=method))
            ## uncomment to track progress of Newton algo
            if(print.verbose){    
                print(paste("fval.c is:", fval.c))
                ## Xbeta   <- c(par$X, par$beta) - alpha * d
                print(paste("fval.c is larger than fval.p", (fval.c > fval.p)))
            }
            ##!! isSymmetric(hess)          
            while(fval.c > fval.p){                
                alpha <- alpha/2
                ## d <- grad
                Xbeta   <- Xbeta.p + alpha * d
                par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1):(length(Xbeta))])
                ## 
                fval.c  <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh, A=A,
                                                       Aobs=Aobs, W=W, tW=tW, order.HD = order.HD, method=method))
                ## 
                tmp <- (fval.c - fval.p)
                ## uncomment to track progress of Newton algo
                if(print.verbose){    
                    print(paste("alpha is:", alpha,", ",
                                "Fc-Fp is:", round(tmp,2)))
                    print(paste("fval.c is:", fval.c))
                }
                if(fval.c < fval.p) print(paste("REDUCED STEP SIZE REACHED Fc < Fp"))
                if(alpha <= 2^(-9))
                {
                    Xbeta <- Xbeta.p
                    par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1):(length(Xbeta.p))])
                    ## 
                    fval.c  <- as.numeric(objective_osc_HD(par=par, theta=par.theta, data=data, mesh=mesh, A=A, Aobs=Aobs, W=W, tW=tW,
                                                           order.HD = order.HD, method=method))
                    ## uncomment to track progress of Newton algo
                    if(print.verbose){    
                        print(paste("STOP MOVING"))
                    }
                    break
                }
            }
            scaled.acc <- (fval.c  %>% abs %>% max(1)) * acc
            tol.test   <- (grad^2) %>% sum %>% sqrt
            if(print.verbose){    
                print(paste("f eval is:", fval.c))
                print(paste("TT is:", round(tol.test,3),
                            "SA is:", round(scaled.acc,3)))
                print(paste("Fp is:", round(fval.p,3),
                            "Fc is:", round(fval.c),3))
            }
            nmv <- (((Xbeta^2) %>% sum %>% sqrt) -
                    ((Xbeta.p^2) %>% sum %>% sqrt)) %>% abs
            if(tol.test < scaled.acc | nmv < acc) break                    
        }
        if(print.verbose){    
            print("EXIT while loop")
        }
        ## ---------------------------------------------------------
        ## ---------------------------------------------------------
        ## pp.llik(data=data, X=par$X, mesh=mesh,
        ##     L=L, beta=par$beta, Aobs=Aobs)
        ## spde  <- inla.spde2.pcmatern(mesh, alpha=3/2,
        ##                              prior.range=c(1e-2,0.05),
        ##                              prior.sigma=c(0.34,0.05))
        lp.Xbeta <- prior.betaX_osc_HD(X=par$X, beta=par$beta, theta=c(par.theta),  mesh=mesh, order.HD = order.HD)
        ##!! check Cholesky in prior 
        Lxbeta.given.thetay    <- Matrix::expand(Matrix::Cholesky(hess, Imult=1000))$L
        logdetQxbeta.given.thetay <- 2*sum(log(diag(Lxbeta.given.thetay)))
        lp.G     <- .5*logdetQxbeta.given.thetay
        llik     <- pp.llik(data=data, X=par$X, mesh=mesh, beta=par$beta, A=A,  Aobs=Aobs, W=W, tW=tW, method=method)
        ## output
        Xest    <<- par$X
        ## print(Xbeta)
        betaest <<- par$beta
        gradest     <<- grad
        Hessianest  <<- hess
        lp.theta    <- prior.theta_osc_HD(rho = exp(par.theta[1]),
                                          sigma = exp(par.theta[2]),
                                          phi = (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
                                          rhoL = hyperpar$rhoL, sigmaU=hyperpar$sigmaU,
                                          rho.HD = exp(par.theta[4]),
                                          sigma.HD = exp(par.theta[5]),
                                          alpha1 = hyperpar$alpha1, alpha2 = hyperpar$alpha2,
                                          alphaHD1 = hyperpar$alpha1, alphaHD2 =hyperpar$alpha2,
                                          rho.HDL=hyperpar$rhoL, sigma.HDU=hyperpar$sigmaU,
                                          lg=TRUE)
        ## lp.theta <- dexp(exp(par.theta[2]), rate=1/0.33, log=TRUE)
        obj <- as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G))
        print(paste("prior.theta: ", round(lp.theta, 3), "prior.Xbeta: ", round(lp.Xbeta, 3),
                    "llik: ", round(llik, 3), "Xbeta.y: ", round(lp.G, 3), "val: ",
                    round(as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G)), 3)))
        ## obj <- as.numeric(-(lp.Xbeta + llik - lp.G))    
        ## save(betaest, Xest, par.theta, gradest, Hessianest, obj,
        ## file = paste0("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating/data/Xbeta", counter, ".RData"))
        if(FALSE){
            save(betaest, Xest, par.theta, gradest, Hessianest, obj, file = paste0("/data/ipapasta/biased/Xbeta", counter, ".RData"))
        }
        cat("\n \n \n")
        counter <<- counter+1
        cat("\n \n \n ")               
        print(paste("counter is equal to:", counter))
        
        return(as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G)))
    }
}
