## -----------------------------------------------------
## grid-fields 
## Laplace approximation for log-Gaussian Cox process
## using Lindgren's SPDE
## mark distribution not done yet
## -----------------------------------------------------

if(TRUE)                                
{
    ldmvnorm <- function(x, m, Q)
    {
        ## ----------------------------------------------------------------
        ## input: x (evaluation point), 
        ##        m (mean vector), 
        ##        Q (precision matrix/inverse covariance Sigma^{-1})
        ## output: log density of multivariate normal distribution
        ## ----------------------------------------------------------------
        k      <- length(x)
        Q      <- Matrix(Q)
        xminmu <- Matrix(x-m,ncol=1)
        L      <- (Matrix::Cholesky(Q) %>% (Matrix::expand))$L
        detQ   <- 2*sum(diag(L))    
        out    <- .5 * log(detQ) -
            (k/2)*log(2*pi) -
            .5*t(xminmu) %*%
            (Q %*% xminmu)    
        return(out)
    }

    
    objective <- function(par, theta, spde, data, mesh, L, A, Aobs) {
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
        X    <- par$X
        beta <- par$beta
        out <-  -prior.Xbeta(X=X, beta=beta, spde=spde, theta=theta, mesh = mesh) -
            pp.llik(data=data, X=X, mesh=mesh, L=L, beta=beta, A=A, Aobs=Aobs)
        return(out)
    }


    objective_osc <- function(par, theta, data, mesh, L, A, Aobs) {
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
        X    <- par$X
        beta <- par$beta
        out <-  -prior.Xbeta_osc(X=X, beta=beta, theta=theta, mesh = mesh) -
            pp.llik(data=data, X=X, mesh=mesh, L=L, beta=beta, A=A, Aobs=Aobs)
        return(out)
    }

    if(FALSE)                           #new version
    {
        ptheta.prop.marg.post2 <- function(par.theta, spde, data, X, beta, mesh, L, A, Aobs,  acc=1e-7)
        {
            ## !! only grad and hessian of X used. estimate of beta |
            ## !! theta, Y replaced by maximum likelihood estimate
            
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
            ##        Approximate posterior of theta given Y
            ##        up to proportionality. This function is 
            ##        optimized with the polytope method (Nelder-Mead)
            ## print("is state maintained?")
            ## cat("\n")
            ## print(X)
            ## cat("\n")
            par.theta   <- exp(par.theta)
            par         <- list(X=X,beta=beta)
            Xbeta       <- rBind(par$X, par$beta) ## starting value
            X           <- Matrix(par$X, ncol=1)
            n.X         <- X %>% length
            betahat     <- NULL
            ## n.beta      <- beta %>% length
            ## 
            ## ---------------------------------------
            ## Newton optimization
            ## ---------------------------------------
            ## starting values
            grad     <- Matrix(gradX.objective(par, theta=par.theta, spde = spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs), ncol=1)
            hess     <- Matrix(hessianX.objective(par, theta=par.theta, spde=spde,
                                                 data=data, mesh=mesh, L=L, A=A, Aobs=Aobs), sparse=TRUE)
            fval.p   <- fval.c <- objective(par=par, theta=par.theta, spde=spde,
                                            data=data, mesh=mesh, L=L, A=A, Aobs=Aobs) %>% as.numeric
            scaled.acc  <- (fval.c %>% abs %>% max(1)) * acc     
            tol.test    <- Inf
            ## uncomment to track progress of Newton algo
            ## d <- -Matrix::solve(hess, grad)
            ## d2 <- -grad
            ## !!!        
            ## 
            print(paste("tol.test is :", tol.test,
                        "scaled.acc is:", scaled.acc))
            ## fdiff <- TRUE
            while(TRUE)    #!!!
            {
                X.p <- X
                alpha   <- 1
                fval.p  <- fval.c
                betahat <- log(n)-log(t(L)%*%exp(A%*%X.p))
                par     <- list(X=X.p, beta=betahat)
                ## 
                ## Newton move
                ## 
                grad  <- Matrix(gradX.objective(par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs), ncol=1)
                hess  <- Matrix(hessianX.objective(par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs),
                                sparse=TRUE)
                ## 
                d       <- -Matrix::solve(hess, grad)
                d2      <- -grad
                X       <- X.p + alpha * d
                betahat <- log(n)-log(t(L)%*%exp(A%*%X.p))
                par     <- list(X=X, beta=betahat)
                ## 
                fval.c  <- as.numeric(objective(par=par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs))
                ## uncomment to track progress of Newton algo
                print(paste("fval.c is:", fval.c))
                print(paste("fval.c is larger than fval.p", (fval.c > fval.p)))
                while(fval.c > fval.p)
                {                
                    alpha <- alpha/2
                    ## d <- grad
                    X     <- X.p + alpha * d2
                    betahat <- log(n)-log(t(L)%*%exp(A%*%X.p))
                    par     <- list(X=X, beta=betahat)
                    ## 
                    fval.c  <- as.numeric(objective(par=par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs))
                    ## 
                    tmp <- (fval.c - fval.p)
                    ## uncomment to track progress of Newton algo
                    print(paste("alpha is:", alpha,", ",
                                "Fc-Fp is:", round(tmp,2)))
                    print(paste("fval.c is:", fval.c))
                    if(alpha <= 2^(-8))
                    {
                        X       <- X.p
                        betahat <- log(n)-log(t(L)%*%exp(A%*%X.p))
                        par     <- list(X=X.p, beta=betahat)
                        ## 
                        fval.c  <- as.numeric(objective(par=par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs))
                        ## uncomment to track progress of Newton algo
                        print(paste("STOP MOVING"))
                        break
                    }
                }
                scaled.acc <- (fval.c  %>% abs %>% max(1)) * acc
                tol.test   <- (grad^2) %>% sum %>% sqrt            
                print(paste("f eval is:", fval.c))
                print(paste("TT is:", round(tol.test,3),
                            "SA is:", round(scaled.acc,3)))
                print(paste("Fp is:", round(fval.p,3),
                            "Fc is:", round(fval.c),3))
                nmv <- (((X^2) %>% sum %>% sqrt) -
                        ((X.p^2) %>% sum %>% sqrt)) %>% abs
                if(tol.test < scaled.acc | nmv < acc)
                {
                    break
                }
            }
            print("EXIT while loop")
            cat("\n \n \n")
            counter <<- counter+1
            print(paste("counter is equal to:",counter))
            cat("\n \n \n ")               
            ## ---------------------------------------------------------
            ## ---------------------------------------------------------
            expAx      <- exp(A%*%X)
            one.p      <-  Matrix(rep(1, nrow(X)), ncol=1)
            hess.beta  <- Matrix(-as.numeric(exp(betahat))*(t(L)%*%expAx))
            hess.Xbeta <- Matrix(-as.numeric(exp(betahat))*(t((L%*%t(one.p))*A) %*% expAx))            
            Lxbeta.given.thetay    <- rBind(cBind((Matrix::Cholesky(hess) %>% (Matrix::expand))$L, -hess.Xbeta),
                                            -cBind(t(hess.Xbeta), -hess.beta))
            detQxbeta.given.thetay <- 2*sum(diag(Lxbeta.given.thetay))
            lp.G     <- .5*log(detQxbeta.given.thetay)    
            lp.Xbeta <- prior.Xbeta(X=par$X, beta=par$beta, spde=spde, theta=c(par.theta),  mesh=mesh)
            llik     <- pp.llik(data=data, X=par$X, mesh=mesh, L=L, beta=par$beta, A=A,  Aobs=Aobs)
            ## output
            Xest       <<- Matrix(par$X, ncol=1)
            betaest    <<- Matrix(par$beta, ncol=1)
            gradest    <<- grad
            Hessianest <<- hess
            return(as.numeric(-(lp.Xbeta + llik - lp.G)))
        }
    }



    
    ptheta.prop.marg.post.out <- function(par.theta, spde, data, X, beta, mesh, L=L, A, Aobs,  acc=1e-7)
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
        ##        Approximate posterior of theta given Y
        ##        up to proportionality. This function is 
        ##        optimized with the polytope method (Nelder-Mead) 
        par.theta   <- exp(par.theta)
        par         <- list(X=X,beta=beta)
        Xbeta       <- rBind(par$X, par$beta) ## starting value
        n.X         <- X %>% length
        ## n.beta      <- beta %>% length
        ## 
        ## ---------------------------------------
        ## Newton optimization
        ## ---------------------------------------
        ## starting values
        grad     <- grad.objective(par, theta=par.theta, spde = spde,
                                   data=data, mesh=mesh, L=L, A=A, Aobs=Aobs) %>% matrix(ncol=1)
        hess     <- Matrix(hessian.objective(par, theta=par.theta, spde=spde,
                                             data=data, mesh=mesh, L=L, A=A, Aobs=Aobs), sparse=TRUE)
        fval.p   <- fval.c <- objective(par=par, theta=par.theta, spde=spde,
                                        data=data, mesh=mesh, L=L, A=A, Aobs=Aobs) %>% as.numeric
        scaled.acc  <- (fval.c %>% abs %>% max(1)) * acc     
        tol.test    <- Inf
        ## uncomment to track progress of Newton algo
        ## d <- -Matrix::solve(hess, grad)
        ## d2 <- -grad
        ## !!!        
        ## 
        print(paste("tol.test is :", tol.test,
                    "scaled.acc is:", scaled.acc))
        ## fdiff <- TRUE
        while(TRUE)    #!!!
        {
            Xbeta.p <- Xbeta
            alpha   <- 1
            fval.p  <- fval.c            
            par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1): (Xbeta.p %>% length)])
            ## 
            ## Newton move
            ## 
            grad  <- grad.objective(par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs) %>% Matrix(ncol=1)
            hess  <- Matrix(hessian.objective(par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs), sparse=TRUE)
            ## 
            d       <- -Matrix::solve(hess,grad)
            d2      <- -grad
            Xbeta   <- Xbeta.p + alpha * d            
            par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1): (Xbeta %>% length)])
            ## 
            fval.c  <- as.numeric(objective(par=par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs))
            ## uncomment to track progress of Newton algo
            print(paste("fval.c is:", fval.c))
            print(paste("fval.c is larger than fval.p", (fval.c > fval.p)))
            while(fval.c > fval.p)
            {                
                alpha <- alpha/2
                ## d <- grad
                Xbeta   <- Xbeta.p + alpha * d2
                par     <- list(X=Xbeta[1:n.X], beta=Xbeta[(n.X+1):(Xbeta %>% length)])
                ## 
                fval.c  <- as.numeric(objective(par=par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs))
                ## 
                tmp     <- (fval.c - fval.p)
                ## uncomment to track progress of Newton algo
                print(paste("alpha is:", alpha,", ",
                            "Fc-Fp is:", round(tmp,2)))
                print(paste("fval.c is:", fval.c))
                if(alpha <= 2^(-8))
                {
                    Xbeta   <- Xbeta.p
                    par     <- list(X=Xbeta.p[1:n.X], beta=Xbeta.p[(n.X+1):(Xbeta.p %>% length)])
                    ## 
                    fval.c  <- as.numeric(objective(par=par, theta=par.theta, spde=spde, data=data, mesh=mesh, L=L, A=A, Aobs=Aobs))
                    ## uncomment to track progress of Newton algo
                    print(paste("STOP MOVING"))
                    break
                }
            }
            scaled.acc <- (fval.c  %>% abs %>% max(1)) * acc
            tol.test   <- (grad^2) %>% sum %>% sqrt            
            print(paste("f eval is:", fval.c))
            print(paste("TT is:", round(tol.test,3),
                        "SA is:", round(scaled.acc,3)))
            print(paste("Fp is:", round(fval.p,3),
                        "Fc is:", round(fval.c),3))
            nmv <- (((Xbeta^2) %>% sum %>% sqrt) - ((Xbeta.p^2) %>% sum %>% sqrt)) %>% abs
            if(tol.test < scaled.acc | nmv < acc)
            {           
                break
            }
        }
        print("EXIT while loop")
        cat("\n \n \n \n \n")
        X    <- par$X
        beta <- par$beta
        ## 
        Lxbeta.given.thetay    <- (Matrix::Cholesky(hess) %>% (Matrix::expand))$L    
        detQxbeta.given.thetay <- 2*sum(diag(Lxbeta.given.thetay))
        lp.G     <- .5*log(detQxbeta.given.thetay)    
        lp.Xbeta <- prior.Xbeta(X=par$X, beta=par$beta, spde=spde, theta=c(par.theta), mesh=mesh)
        llik     <- pp.llik(data=data, X=par$X, mesh=mesh, L=L, beta=par$beta, A=A, Aobs=Aobs)
        ## output
        negval   <- -as.numeric(lp.Xbeta + llik - lp.G)
        gradient <- grad
        Hessian  <- hess
        theta    <- exp(par.theta)
        out      <- list(X, beta, theta, gradient, Hessian)
        names(out) <- c("mx.thetay", "mbeta.thetay", "theta.y", "grad", "QXbeta.thetay")
        out %>% return
    }


    
}


if(FALSE)
{
    X.scaled <- rank((Xbeta[1:nrow(mesh$loc)]))/
        (mesh$loc %>% nrow)

    pdf(file="~/Desktop/modelfit.pdf")
    ## mesh with spatial effects
    plot(mesh, asp=1, main= "spatial effect - model fit")    
    points(mesh$loc, col=gray(1-X.scaled),
           pch=16, cex=1.2,asp=1)
    dev.off()
    ## Linde accidents
    pdf(file="~/Desktop/modelfit_and_linde_data.pdf")
    plot(mesh, asp=1,main= "Linde data - red points")    
    points(mesh$loc,col=gray(1-X.scaled),
           pch=16, cex=1.2,asp=1)
    points(coo, cex=0.5, col=2, pch=16)
    dev.off()
    ##
    pdf(file="~/Desktop/modelfit_and_locations_to_predict.pdf")
    plot(mesh, asp=1,
         main= "Locations to predict - blue points")    
    points(mesh$loc, col=gray(1-X.scaled),
           pch=16, cex=1.2,asp=1)
    points(new.data$Long,new.data$Lat,
           col=4, pch=16, cex=1)
    dev.off()




    ## intensity

    plot(mesh, asp=1, main= "spatial effect - model fit")    
    points(mesh$loc, col=lambda,
           pch=16, cex=1.2,asp=1)


    ## 
    points(coo, cex=1, col=2, pch=16)
    plot(mesh, asp=1)
    points(mesh$loc, col=gray(1-X.scaled),
           pch=16, cex=2,asp=1)
    points(coo, col=2, pch=16, cex=1)
    points(new.data$Long,new.data$Lat,
           col=4, pch=16, cex=2)
    
}




