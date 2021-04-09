'oscillating.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        rho     <- exp(theta[1L])
        kappa   <- sqrt(8)/rho
        sigma   <- exp(theta[2L])
        phi     <- (1-exp(-theta[3L]))/(1+exp(-theta[3L]))
        sincpth <- sqrt(1-phi^2)/acos(phi)
        tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
        z       <- list(tausq = tausq, kappa  = kappa, phi = phi)
        return(z)
    }
    graph <- function(){
        return(M$M2)
    }
    Q <- function() {
        require(Matrix)
        param <- interpret.theta()
        precision <- param$tausq*(param$kappa^2 * M$M0 + 2*param$phi * param$kappa^2 * M$M1 + M$M2)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        sigmaLN  <- 0.5
        murho    <- 25
        rho      <- sqrt(8)/param$kappa
        sigma    <- sqrt(param$tausq)
        phi      <- param$phi
        lrho.sp    <- dlnorm(rho, log(murho), sigmaLN, log=TRUE)    
        lsigma.sp  <- dexp(sigma, 1/2, log = TRUE)
        lpphi.sp   <- 0 ## prior.phi_osc(phi, a=1, b=20, lg=TRUE)
        res        <- lpphi.sp + lrho.sp + lsigma.sp
        return(res)
    }
    initial <- function()  return(c(0, 0, 0))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}


## oscillating.rgeneric <- inla.rgeneric.define(oscillating.model) #arguments that need to be passed are M0, M1 and M2 

'circular1D.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
                               theta = NULL) {
    interpret.theta <- function() {
        kappa    <- exp(theta[1L])
        sigma   <- exp(theta[2L])        
        tausq   <- (pi + sinh(2 * pi * kappa)/(2*kappa))/((sinh(pi*kappa)^2)*(sigma^2)*(2*kappa.dir)^2)
        z       <- list(kappa  = kappa, tausq = tausq)
        return(z)
    }
    graph <- function() require(Matrix); return(Q())
    Q <- function() {
        require(Matrix)
        param <- interpret.theta()
        precision <- param$tausq*(param$kappa^2 * M0 + 2 * param$kappa^2 * M1 + M2)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {
        param = interpret.theta()
        sigma.hd <- 1/sqrt(param$tausq)
        rho.hd      <- sqrt(8)/param$kappa
        lrho.hd    <- dexp(rho.hd, 1/pi, log=TRUE)        
        lsigma.hd  <- dexp(sigma.hd, 1, log = TRUE)
        if(FALSE){
            rpi.div.r0 <- 2*(pi*param$kappa*cosh(pi*param$kappa)+sinh(pi*param$kappa))/(2*pi*param$kappa + sinh(2*pi*param$kappa))
            lrpi.div.r0.hd <- dbeta(rbi.div.r0, alpha=1, beta=10)
        }
        res <- lrho.hd + lsigma.hd
        return(res)
    }
    initial <- function()  return(c(0, 0))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}

circular1D <- inla.rgeneric.define(circular1D.model, M0=NULL, M1=NULL, M2=NULL) 

