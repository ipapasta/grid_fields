'oscillating.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        rho     <- 5 + exp(theta[1L])
        kappa   <- sqrt(8)/rho
        sigma   <- exp(theta[2L])
        ## phi     <- (1-exp(-theta[3L]))/(1+exp(-theta[3L]))
        phi     <- 1/(1+exp(-theta[3L]))-1
        sincpth <- sqrt(1-phi^2)/acos(phi)
        tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
        z       <- list(tausq = tausq, rho  = rho, phi = phi)
        return(z)
    }
    graph <- function() return(M$M2)
    Q <- function() {
        require(Matrix)
        param     <- interpret.theta()
        kappa     <- sqrt(8)/param$rho
        precision <- param$tausq*(kappa^4 * M$M0 + 2*param$phi * kappa^2 * M$M1 + M$M2)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        prior.phi_osc <- function(phi, a, b, l=(-0.998), u=0, lg=TRUE){
            if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
            if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
        }
        rho        <- sqrt(8)/param$kappa
        sigma      <- sqrt(param$tausq)
        phi        <- param$phi
        ## distribution of hyperparameters         
        sigma.spatial.oscillating <- hyperpar$sigma.spatial.oscillating
        lrho.sp    <- dlnorm(rho-5, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
        lsigma.sp  <- dexp(sigma, sigma.spatial.oscillating, log = TRUE)
        lpphi.sp <- dunif(phi, -1, 0)
        ## lpphi.sp   <- prior.phi_osc(phi,
        ##                             a=hyperpar$a.par.phi.prior.spatial.oscillating,
        ##                             b=hyperpar$b.par.phi.prior.spatial.oscillating,
        ##                             lg=TRUE)
        res        <- lpphi.sp + lrho.sp + lsigma.sp
        return(res)
    }
    initial <- function()  return(c(initial.space$theta1, initial.space$theta2, initial.space$theta3))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}


'space.direction.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        ## space parameters
        rho.space     <- exp(theta[1L])
        kappa.space   <- sqrt(8)/rho.space
        sigma.space   <- exp(theta[2L])
        phi.space     <- (1-exp(-theta[3L]))/(1+exp(-theta[3L]))
        ## tausq.space   <- 1/(4*pi*(sigma.space^2)*(kappa.space^2)*sincpth)
        ## direction parameters
        rho.direction     <- exp(theta[4L])
        kappa.direction   <- sqrt(8*(3/2))/rho.direction
        sigma.direction   <- exp(theta[5L])
        ## 
        z       <- list(sigma.space = sigma.space,
                        rho.space  = rho.space,
                        phi.space = phi.space,
                        sigma.direction = sigma.direction,
                        rho.direction=rho.direction)
        return(z)
    }
    graph <- function() return(kronecker(M$M2.direction, M$M2.space))
    Q <- function() {
        require(Matrix)
        param               <- interpret.theta()
        kappa.space         <- sqrt(8)/param$rho.space
        kappa.direction     <- sqrt(8*(3/2))/param$rho.direction
        sincpth             <- sqrt(1-param$phi.space^2)/acos(param$phi.space)
        tausq.space         <- 1/(4*pi*(param$sigma.space^2)*(kappa.space^2)*sincpth)
        tausq.direction     <- 1/(4*(param$sigma.direction^2)*(kappa.direction^3))
        precision.space     <- tausq.space*(kappa.space^4 * M$M0.space + 2*param$phi.space * kappa.space^2 * M$M1.space + M$M2.space)
        precision.direction <- tausq.direction*(kappa.direction^4 * M$M0.direction + 2 * kappa.direction^2 * M$M1.direction + M$M2.direction)
        precision <- kronecker(precision.direction, precision.space)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        prior.phi_osc <- function(phi, a, b, l=(-0.998), u=1, lg=TRUE){
            if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
            if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
        }
        sigmaLN  <- 3
        murho    <- 30
        sigma.space               <- param$sigma.space
        phi.space                 <- param$phi.space
        sigma.direction           <- param$sigma.direction
        sigmaLN                   <- hyperpar$sigma.range.spatial.oscillating
        murho                     <- hyperpar$mu.range.spatial.oscillating
        sigma.spatial.oscillating <- hyperpar$sigma.spatial.oscillating
        lrho.space                <- dlnorm(rho, log(murho), sigmaLN, log=TRUE)
        lsigma.space              <- dexp(sigma, sigma.spatial.oscillating, log = TRUE)       
        ## lrho.space        <- dlnorm(param$rho.space, log(murho), sigmaLN, log=TRUE)
        ## lsigma.space      <- dexp(sigma.space, 1/2, log = TRUE)
        lpphi.space       <- prior.phi_osc(phi.space,
                                           a=hyperpar$a.par.phi.prior.spatial.oscillating,
                                           b=hyperpar$b.par.phi.prior.spatial.oscillating,
                                           lg=TRUE)
        lrho.direction    <- dexp(param$rho.direction, hyperpar$rho.directional, log=TRUE)   
        lsigma.direction  <- dexp(sigma.direction, hyperpar$sigma.directional, log = TRUE)
        res               <- lpphi.space + lrho.space + lsigma.space + lrho.direction + lsigma.direction
        return(res)
        ## rho0     <- 30
        ## lrho.sp    <- dlnorm(rho, log(murho), sigmaLN, log=TRUE)
        ## alpha      <- 1-10^(-15)
        ## lambda.rho <- -log(alpha) * rho0 
        ## lrho.sp    <- log(lambda.rho) - 2*log(rho) -lambda.rho/rho
    }
    initial <- function()  return(c(log(21), log(0.8), atanh(-.92), log(3.22), log(1.86)))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}

'temporal.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        rho     <- exp(theta[1L])
        kappa   <- sqrt(8*(3/2))/rho
        sigma   <- exp(theta[2L])
        z       <- list(sigma.temporal = sigma, rho.temporal  = rho)
        return(z)
    }
    graph <- function() return(M$M2.temporal)
    Q <- function() {
        require(Matrix)
        param            <- interpret.theta()        
        kappa.temporal   <- sqrt(8*(3/2))/param$rho.temporal
        tausq.temporal   <- 1/(4*(param$sigma.temporal^2)*(kappa.temporal^3))
        precision        <- tausq.temporal*(kappa.temporal^4 * M$M0.temporal + 2 * kappa.temporal^2 * M$M1.temporal + M$M2.temporal)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()        
        lrho.temporal    <- dexp(param$rho.temporal, hyperpar$rho.temporal , log=TRUE)
        lsigma.temporal  <- dexp(param$sigma.temporal, hyperpar$sigma.temporal, log = TRUE)
        res              <- lrho.temporal + lsigma.temporal
        return(res)
    }
    initial <- function()  return(c(log(22), log(0.68)))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}

space <- function(time){
    f <- function(time.scalar){
        dat <<- data$Ypos
        if(time.scalar <= dat$time[1]) return(dat$coords[1,])
        if(time.scalar >= dat$time.lead[length(dat$time.lead)]) return(dat$coords[nrow(dat),])
        else{
            wh.start <- max(which(dat$time <= time.scalar))
            bary <- (dat$time.lead[wh.start] - time.scalar)/(dat$time.lead[wh.start]- dat$time[wh.start])
            return(bary*dat$coords[wh.start,] + (1-bary)*dat$coords.lead[wh.start,])
        }
    }
    t((Vectorize(f, vectorize.args="time.scalar"))(time))
}

direction <- function(time){
    f <- function(time.scalar){
        dat <<- data$Ypos
        wh.start <- max(which(dat$time <= time.scalar))
        bary <- (dat$time.lead[wh.start] - time.scalar)/(dat$time.lead[wh.start]- dat$time[wh.start])
        bary*dat$hd[wh.start] + (1-bary)*dat$hd.lead[wh.start,]
    }
    (Vectorize(f, vectorize.args="time.scalar"))(time)
}


## oscillating.rgeneric <- inla.rgeneric.define(oscillating.model) #arguments that need to be passed are M0, M1 and M2 

'circular1D.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
                               theta = NULL) {
    envir <- parent.env(environment())
    interpret.theta <- function() {
        kappa    <- exp(theta[1L])
        sigma   <- exp(theta[2L])        
        tausq   <- (pi + sinh(2 * pi * kappa)/(2*kappa))/((sinh(pi*kappa)^2)*(sigma^2)*(2*kappa.dir)^2)
        z       <- list(kappa  = kappa, tausq = tausq)
        return(z)
    }
    graph <- function() return(M2)
    Q <- function() {
        require(Matrix)
        param <- interpret.theta()
        precision <- param$tausq*(param$kappa^4 * M0 + 2 * param$kappa^2 * M1 + M2)
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

## circular1D <- inla.rgeneric.define(circular1D.model, M0=NULL, M1=NULL, M2=NULL) 

