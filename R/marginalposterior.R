## par.theta=log.par.theta; acc=1e-7
if(FALSE){
    list.functions.in.file("marginalposterior.R")

    "objective_osc_temp"
    "prior.betaXZ_osc_temp"
    "grad.objective_osc_temp"
    "hessian.objective_osc_temp"
    "pp.llik"
    "prior.theta_osc"
}



## W.input.to.sparseMatrix, 
pthetapc.prop.marg.post_osc_temp<- function(par.theta, data, X, Z, beta, mesh, mesh.theta, mesh1d, A, Atilde, Aobs, Atildeobs, W,
                                            acc=1e-4, print.verbose=TRUE)
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
    print(paste("rho is:", exp(par.theta[1]), " sigma is:", exp(par.theta[2]),
                "phi is: ", (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
                "rho_hd is:", exp(par.theta[4]), " sigma_hd is:", exp(par.theta[5]),
                "rho_temporal is: ", exp(par.theta[6]), "sigma_temporal is: ", exp(par.theta[7])))   
    par.theta   <- par.theta
    par         <- list(X=X, Z=Z, beta=beta)
    betaXZ      <- rBind(par$beta, par$X, par$Z) ## starting value
    n.X         <- length(X)
    n.Z         <- length(Z)
    ## n.beta      <- beta %>% length
    ## 
    ## ---------------------------------------
    ## Newton optimization
    ## ---------------------------------------
    ## starting values    
    grad     <- grad.objective_osc_temp(par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs,W=W)
    hess     <- Matrix(hessian.objective_osc_temp(par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,
                                                  Atildeobs=Atildeobs,W=W),sparse=TRUE)
    ## ,W.input.to.sparseMatrix=W.input.to.sparseMatrix
    fval.p   <- fval.c <- as.numeric(objective_osc_temp(par=par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,mesh1d=mesh1d,
                                                        A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs,W=W))
    scaled.acc  <- (max( abs(fval.c), 1)) * acc  
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
    while(TRUE) {    #!!!{
        betaXZ.p <- betaXZ
        alpha   <- 1
        fval.p  <- fval.c            
        par     <- list(beta=betaXZ.p[1], X=betaXZ.p[2:(1+n.X)], Z=betaXZ.p[((1+n.X)+1):length(betaXZ.p)] )
        ## 
        ## Newton move
        ## 
        grad  <- grad.objective_osc_temp(par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,
                                                mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs,W=W)
        hess  <- hessian.objective_osc_temp(par, theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,
                                            mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs, W=W)
        d       <- -Matrix::solve(hess, grad)
        betaXZ  <- betaXZ.p  + alpha * d
        par     <- list(beta=betaXZ[1], X=betaXZ[2:(1+n.X)], Z=betaXZ[((1+n.X)+1):length(betaXZ)] )
        ## print(par$X)
        ## 
        fval.c  <- as.numeric(objective_osc_temp(par=par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,
                                                 mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs,W=W))
        ## uncomment to track progress of Newton algo
        if(print.verbose){    
            print(paste("fval.c is:", fval.c))
            ## Xbeta   <- c(par$X, par$beta) - alpha * d
            print(paste("fval.c is larger than fval.p", (fval.c > fval.p)))
        }
        if(TRUE){
            while(fval.c > fval.p) {                
                alpha <- alpha/2
                ## d <- grad
                betaXZ  <- betaXZ.p  + alpha * d
                par     <- list(beta=betaXZ[1], X=betaXZ[2:(1+n.X)], Z=betaXZ[((1+n.X)+1):length(betaXZ)] )
                ## 
                fval.c  <- as.numeric(objective_osc_temp(par=par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,
                                                         mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs,W=W))
                ## 
                tmp <- (fval.c - fval.p)
                ## uncomment to track progress of Newton algo
                if(print.verbose){    
                    print(paste("alpha is:", alpha,", ",
                                "Fc-Fp is:", round(tmp,4)))
                    print(paste("fval.c is:", fval.c))
                }
                if(alpha <= 2^(-12))
                {
                    betaXZ  <- betaXZ.p 
                    par     <- list(beta=betaXZ.p[1],
                                    X=betaXZ.p[2:(1+n.X)],
                                    Z=betaXZ.p[((1+n.X)+1):length(betaXZ.p)] )
                    ## 
                    fval.c  <- as.numeric(objective_osc_temp(par=par,theta=par.theta,data=data,mesh=mesh,mesh.theta=mesh.theta,
                                                             mesh1d=mesh1d,A=A,Atilde=Atilde,Aobs=Aobs,Atildeobs=Atildeobs,W=W))
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
        nmv <- (((betaXZ^2) %>% sum %>% sqrt) -
                ((betaXZ.p^2) %>% sum %>% sqrt)) %>% abs
        if(tol.test < scaled.acc | nmv < acc) break                    
    }
    if(print.verbose){    
        print("EXIT while loop")
    }
    ## ---------------------------------------------------------
    LbetaXZ.given.thetay    <- Matrix::expand(Matrix::Cholesky(hess, LDL=TRUE))$L
    logdetQbetaXZ.given.thetay <- 2*sum(log(diag(LbetaXZ.given.thetay)))
    ## print(paste("sum of lof diagonal elements of Cholesky is:", logdetQbetaXZ.given.thetay))
    lp.G     <- .5*logdetQbetaXZ.given.thetay
    lp.betaXZ <- prior.betaXZ_osc_temp(beta=par$beta, X=par$X, Z=par$Z, theta=c(par.theta),  mesh=mesh, mesh.theta=mesh.theta, mesh1d=mesh1d)
    llik      <- pp.llik(data=data, X=par$X, Z=par$Z, mesh=mesh, mesh1d=mesh1d, beta=par$beta, A=A, Atilde, Aobs=Aobs, Atildeobs=Atildeobs, W=W)
    betaest <<- par$beta
    Xest    <<- par$X
    Zest    <<- par$Z
    gradest     <<- grad
    Hessianest  <<- hess
    lp.theta    <- prior.theta_osc_temp(rho=exp(par.theta[1]), sigma=exp(par.theta[2]), phi = (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
                                        rho.hd=exp(par.theta[4]), sigma.hd=exp(par.theta[5]),
                                        rho.temp=exp(par.theta[6]), sigma.temp=exp(par.theta[7]))
    ## (rho=exp(par.theta[1]), sigma=exp(par.theta[2]), phi=(1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
    ##     rhoL=hyperpar$rhoL, sigmaU=hyperpar$sigmaU, alpha1=hyperpar$alpha1, alpha2=hyperpar$alpha2, # hyperparameters of spatial component
    ##     rho.temp=exp(par.theta[4]), sigma.temp=exp(par.theta[5]), rho.tempL=hyperpar$rho.tempL, sigma.tempU=hyperpar$sigma.tempU, alphaT1, alpha2T,
    ##     lg=TRUE)
    obj <- as.numeric(-(lp.theta + lp.betaXZ + llik - lp.G))
    ## save(betaest, Xest, par.theta, gradest, Hessianest, obj,
    ## file = paste0("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating/data/Xbeta", counter, ".RData"))
    ## if(FALSE){
        ## if(getwd()=="/home/ipapasta/Software/R/grid_fields/R/Temporal_modulation/"){
    ## save(betaest, Xest, Zest, par.theta, gradest, Hessianest, obj, file = paste0("/data/ipapasta/Xbeta_coarse", counter, ".RData"))
    ##     }
    ## }
    ## save(betaest, Xest, par.theta, gradest, Hessianest, obj, file = paste0("/home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/Oscillating/data/mouse14_05_16/Xbeta", counter, ".RData"))
    print(paste("prior.theta: ", round(lp.theta, 3), "prior.betaXZ: ", round(lp.betaXZ, 3),
                "llik: ", round(llik, 3), "betaXZ.y: ", round(lp.G, 3), "val: ",
                round(as.numeric(-(lp.theta + lp.betaXZ + llik - lp.G)), 3)))
    cat("\n \n \n")
    ## double assignment
    counter <<- counter+1
    ## double assignment
    cat("\n \n \n ")               
    print(paste("counter is equal to:",counter))
    return(as.numeric(-(lp.theta + lp.betaXZ + llik - lp.G)))
}



