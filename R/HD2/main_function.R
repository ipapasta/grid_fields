pthetapc.prop.marg.post_osc_HD <- function(par.theta, hyperpar, data,
                                           X, beta, mesh, L, A, W, Aobs,
                                           type="biased", acc=1e-7, print.verbose=FALSE)
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
    grad     <- matrix(grad.objective_osc(par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, method=method), ncol=1)
    hess     <- Matrix(hessian.objective_osc(par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, method=method), sparse=TRUE)
    fval.p   <- fval.c <- as.numeric(objective_osc(par=par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs, W=W, method=method))  # !! method
    scaled.acc  <- (max( abs(fval.c), 1)) * acc  ## (fval.c %>% abs %>% max(1)) * acc     
    tol.test    <- Inf
    ## uncomment to track progress of Newton algorithm
    ## d <- -Matrix::solve(hess, grad)
    ## d2 <- -grad
    ## !!!        
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
        grad  <- grad.objective_osc(par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, method=method) %>% Matrix(ncol=1)
        hess  <- Matrix(hessian.objective_osc(par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, W=W, Aobs=Aobs, method=method), sparse=TRUE)
        ## 
        d       <- -Matrix::solve(hess, grad)
        ## d2      <- -grad
        Xbeta   <- Xbeta.p + alpha * d            
        par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1): (length(Xbeta))])
        ## print(par$X)
        ## 
        fval.c  <- as.numeric(objective_osc(par=par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs, W=W, method=method))
        ## uncomment to track progress of Newton algo
        if(print.verbose){    
            print(paste("fval.c is:", fval.c))
            ## Xbeta   <- c(par$X, par$beta) - alpha * d
            print(paste("fval.c is larger than fval.p", (fval.c > fval.p)))
        }
        if(TRUE){
            while(fval.c > fval.p){                
                alpha <- alpha/2
                ## d <- grad
                Xbeta   <- Xbeta.p + alpha * d
                par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1):(length(Xbeta))])
                ## 
                fval.c  <- as.numeric(objective_osc(par=par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs, W=W, method=method))
                ## 
                tmp <- (fval.c - fval.p)
                ## uncomment to track progress of Newton algo
                if(print.verbose){    
                    print(paste("alpha is:", alpha,", ",
                                "Fc-Fp is:", round(tmp,2)))
                    print(paste("fval.c is:", fval.c))
                }
                if(alpha <= 2^(-8))
                {
                    Xbeta <- Xbeta.p
                    par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1):(length(Xbeta.p))])
                    ## 
                    fval.c  <- as.numeric(objective_osc(par=par, theta=par.theta, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs, W=W, method=method))
                    ## uncomment to track progress of Newton algo
                    if(print.verbose){    
                        print(paste("STOP MOVING"))
                    }
                    break
                }
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
    Lxbeta.given.thetay    <- (Matrix::Cholesky(hess) %>% (Matrix::expand))$L    
    logdetQxbeta.given.thetay <- 2*sum(log(diag(Lxbeta.given.thetay)))
    lp.G     <- .5*logdetQxbeta.given.thetay
    lp.Xbeta <- prior.Xbeta_osc(X=par$X, beta=par$beta, theta=c(par.theta),  mesh=mesh)
    llik     <- pp.llik(data=data, X=par$X, mesh=mesh, L=L, beta=par$beta, A=A,  Aobs=Aobs, W=W, method=method)
    ## output
    Xest    <<- par$X
    ## print(Xbeta)
    betaest <<- par$beta
    gradest     <<- grad
    Hessianest  <<- hess
    print(paste("rho is:", exp(par.theta[1]), " sigma is:", exp(par.theta[2]),
                "phi is: ", (1-exp(-par.theta[3]))/(1+exp(-par.theta[3]))))
    lp.theta    <- prior.theta_osc_HD(rho=(2*pi)*(exp(par.theta[1])/(1+exp(par.theta[1]))),
                                   sigma=exp(par.theta[2]), phi = (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
                                   rhoL=hyperpar$rhoL, sigmaU=hyperpar$sigmaU,
                                   alpha1=hyperpar$alpha1, alpha2=hyperpar$alpha2,
                                   lg=TRUE)
    obj <- as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G))
    ## save(betaest, Xest, par.theta, gradest, Hessianest, obj,
    ## file = paste0("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating/data/Xbeta", counter, ".RData"))
    if(getwd()=="/home/ipapasta/Software/R/grid_fields/R/Oscillating"){

        if(type=="biased"){
            save(betaest, Xest, par.theta, gradest, Hessianest, obj, file = paste0("/data/ipapasta/biased/Xbeta", counter, ".RData"))
        }else{
            save(betaest, Xest, par.theta, gradest, Hessianest, obj, file = paste0("/data/ipapasta/origin/Xbeta", counter, ".RData"))
        }
    }
    ## save(betaest, Xest, par.theta, gradest, Hessianest, obj, file = paste0("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating/data/mouse14_05_16/Xbeta", counter, ".RData"))
    print(paste("prior.theta: ", round(lp.theta, 3), "prior.Xbeta: ", round(lp.Xbeta, 3),
                "llik: ", round(llik, 3), "Xbeta.y: ", round(lp.G, 3), "val: ",
                round(as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G)), 3)))
    cat("\n \n \n")
    counter <<- counter+1
    cat("\n \n \n ")               
    print(paste("counter is equal to:",counter))
    return(as.numeric(-(lp.theta + lp.Xbeta + llik - lp.G)))
}

