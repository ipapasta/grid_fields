'oscillating.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        rho      <- theta.functions$theta.2.rho(theta[1L])
        sigma    <- theta.functions$theta.2.sigma(theta[2L])
        phi      <- theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u)
        kappa    <- sqrt(8)/rho
        sincpth  <- sqrt(1-phi^2)/acos(phi)
        tausq    <- 1/(4*pi*(kappa^2)*(sigma^2)*sincpth)
        ljac.rho   <- attr(rho, "ljacobian")
        ljac.sigma <- attr(sigma, "ljacobian")
        ljac.phi   <- attr(phi, "ljacobian")
        ##
        z        <- list(rho  = rho, tausq = tausq, phi = phi, sigma=sigma, kappa=kappa,
                         ljac.rho    = ljac.rho, ljac.sigma  = ljac.sigma, ljac.phi    = ljac.phi)
        return(z)
    }
    graph <- function() return(M$M2)
    Q <- function() {
        require(Matrix)
        param     <- interpret.theta()
        precision <- param$tausq*(param$kappa^4 * M$M0 + 2*param$phi * param$kappa^2 * M$M1 + M$M2)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        lrho.sp    <- dlnorm(param$rho-5, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
        lsigma.sp  <- dexp(param$sigma, hyperpar$sigma.spatial.oscillating, log = TRUE)
        lpphi.sp   <- prior.functions$prior.phi_osc(param$phi,
                                                    a=hyperpar$a.par.phi.prior.spatial.oscillating,
                                                    b=hyperpar$b.par.phi.prior.spatial.oscillating,
                                                    l=theta.functions$l, u=theta.functions$u,
                                                    lg=TRUE)
        res        <- lpphi.sp + lrho.sp + lsigma.sp + param$ljac.phi + param$ljac.sigma + param$ljac.rho
        return(res)
    }
    initial <- function()  return(c(initial.space$theta1, initial.space$theta2, initial.space$theta3))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}

'oscillating.model_fixed.phi' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        rho      <- theta.functions$theta.2.rho(theta[1L])
        sigma    <- theta.functions$theta.2.sigma(theta[2L])
        kappa    <- sqrt(8)/rho
        sincpth  <- sqrt(1-hyperpar$phi^2)/acos(hyperpar$phi)
        tausq    <- 1/(4*pi*(kappa^2)*(sigma^2)*sincpth)
        ##
        ljac.rho   <- attr(rho, "ljacobian")
        ljac.sigma <- attr(sigma, "ljacobian")
        z        <- list(tausq = tausq, rho  = rho, kappa=kappa, sigma=sigma,
                         ljac.sigma=ljac.sigma, ljac.rho = ljac.rho )
        return(z)
    }
    Q <- function() {
        require(Matrix)
        param     <- interpret.theta()
        precision <- param$tausq*(param$kappa^4 * M$M0 + 2*hyperpar$phi * param$kappa^2 * M$M1 + M$M2)
        return(precision)
    }
    graph <- function() return(M$M2)
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function(){        
        param      = interpret.theta()
        ## sigma      <- sqrt(1/tausq)
        lrho.sp    <- dlnorm(param$rho-5, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
        lsigma.sp  <- dexp(param$sigma, hyperpar$sigma.spatial.oscillating, log = TRUE)
        res        <- lrho.sp + lsigma.sp + param$ljac.sigma + param$ljac.rho
        return(res)
    }
    initial <- function()  return(c(initial.space$theta1, initial.space$theta2))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}

'space.direction.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        ## space parameters
        rho.space            <- theta.functions$theta.2.rho(theta[1L])
        sigma.space          <- theta.functions$theta.2.sigma(theta[2L])
        phi.space            <- theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u)
        kappa.space          <- sqrt(8)/rho.space
        sincpth              <- sqrt(1-phi.space^2)/acos(phi.space)
        tausq.space          <- 1/(4*pi*(kappa.space^2)*(sigma.space^2)*sincpth)
        ljac.rho.space       <- attr(rho.space, "ljacobian")
        ljac.sigma.space     <- attr(sigma.space, "ljacobian")
        ljac.phi.space       <- attr(phi.space, "ljacobian")
        rho.direction        <- theta.functions$theta.2.rho.direction(theta[4L])
        kappa.direction      <- sqrt(8*(3/2))/rho.direction
        sigma.direction      <- theta.functions$theta.2.sigma.direction(theta[5L])
        tausq.direction      <- 1/(4*(sigma.direction^2)*(kappa.direction^3))
        ljac.rho.direction   <- attr(rho.direction, "ljacobian")
        ljac.sigma.direction <- attr(sigma.direction, "ljacobian")
        ## 
        z       <- list(sigma.space           = sigma.space,
                        rho.space             = rho.space,
                        phi.space             = phi.space,
                        kappa.space           = kappa.space,
                        tausq.space           = tausq.space,
                        sigma.direction       = sigma.direction,
                        rho.direction         = rho.direction,
                        kappa.direction       = kappa.direction,
                        tausq.direction       = tausq.direction,
                        ljac.rho.space        = ljac.rho.space,
                        ljac.sigma.space      = ljac.sigma.space,
                        ljac.phi.space        = ljac.phi.space,
                        ljac.rho.direction    = ljac.rho.direction, 
                        ljac.sigma.direction  = ljac.sigma.direction)
        return(z)
    }
    graph <- function() return(kronecker(M$M2.direction, M$M2.space))
    Q <- function() {
        require(Matrix)
        param               <- interpret.theta()
        precision.space     <- param$tausq.space*(param$kappa.space^4 * M$M0.space + 2*param$phi.space * param$kappa.space^2 * M$M1.space + M$M2.space)
        precision.direction <- param$tausq.direction*(param$kappa.direction^4 * M$M0.direction + 2 * param$kappa.direction^2 * M$M1.direction + M$M2.direction)
        precision <- kronecker(precision.direction, precision.space)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        lrho.space                   <- dlnorm(param$rho.space, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
        lsigma.space                 <- dexp(param$sigma.space, hyperpar$sigma.spatial.oscillating, log = TRUE)       
        lpphi.space       <- prior.functions$prior.phi_osc(param$phi.space,
                                                           a=hyperpar$a.par.phi.prior.spatial.oscillating,
                                                           b=hyperpar$b.par.phi.prior.spatial.oscillating,
                                                           lg=TRUE)
        lrho.direction    <- dexp(param$rho.direction, hyperpar$rho.directional, log=TRUE)   
        lsigma.direction  <- dexp(param$sigma.direction, hyperpar$sigma.directional, log = TRUE)
        res               <- lpphi.space + lrho.space + lsigma.space + lrho.direction + lsigma.direction +
            param$ljac.phi.space + param$ljac.sigma.space + param$ljac.rho.space +
            param$ljac.rho.direction + param$ljac.sigma.direction 
        return(res)
    }
    initial <- function()  return(c(initial.space$theta1, initial.space$theta2, initial.space$theta3, initial.direction$theta4, initial.direction$theta5))
    quit    <- function()  return(invisible())
    res     <- do.call(match.arg(cmd), args = list())
    return(res)
}



'space.direction2.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        ## space parameters
        rho.space            <- theta.functions$theta.2.rho(theta[1L])
        sigma.space          <- theta.functions$theta.2.sigma(theta[2L])
        phi.space            <- theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u)
        kappa.space          <- sqrt(8)/rho.space
        sincpth              <- sqrt(1-phi.space^2)/acos(phi.space)
        tausq.space          <- 1/(4*pi*(kappa.space^2)*(sigma.space^2)*sincpth)
        ljac.rho.space       <- attr(rho.space, "ljacobian")
        ljac.sigma.space     <- attr(sigma.space, "ljacobian")
        ljac.phi.space       <- attr(phi.space, "ljacobian")
        rho.direction        <- theta.functions$theta.2.rho.direction(theta[4L])
        kappa.direction      <- sqrt(8*(3/2))/rho.direction
        sigma.direction      <- theta.functions$theta.2.sigma.direction(theta[5L])
        phi.direction        <- theta.functions$theta.2.phi(theta[6L], l=theta.functions$l, u=theta.functions$u)
        ## -----------------------------------
        ## variance of damped circular field
        ## -----------------------------------
        a                    <- -phi.direction; b          <- -(1-phi.direction^2)^(1/2)        
        x                    <- cos(atan2(y=b, x=a)/2); y  <- sin(atan2(y=b, x=a)/2)
        tausq.direction      <- ((y*sin(2*kappa.direction*pi*x)+x*sinh(2*kappa.direction*pi*y)))/
            (2*b*(kappa.direction^3)*sigma.direction*(cosh(2*kappa.direction*pi*y)-cos(2*kappa.direction*pi*x)))         
        ## tausq.direction      <- 1/(4*(sigma.direction^2)*(kappa.direction^3))
        ljac.rho.direction   <- attr(rho.direction, "ljacobian")
        ljac.sigma.direction <- attr(sigma.direction, "ljacobian")
        ljac.phi.direction   <- attr(phi.direction, "ljacobian")
        ## 
        z       <- list(sigma.space           = sigma.space,
                        rho.space             = rho.space,
                        phi.space             = phi.space,
                        kappa.space           = kappa.space,
                        tausq.space           = tausq.space,
                        sigma.direction       = sigma.direction,
                        rho.direction         = rho.direction,
                        phi.direction         = phi.direction,
                        kappa.direction       = kappa.direction,
                        tausq.direction       = tausq.direction,
                        ljac.rho.space        = ljac.rho.space,
                        ljac.sigma.space      = ljac.sigma.space,
                        ljac.phi.space        = ljac.phi.space,
                        ljac.rho.direction    = ljac.rho.direction, 
                        ljac.sigma.direction  = ljac.sigma.direction,
                        ljac.phi.direction    = ljac.phi.direction)
        return(z)
    }
    graph <- function() return(kronecker(M$M2.direction, M$M2.space))
    Q <- function() {
        require(Matrix)
        param               <- interpret.theta()
        precision.space     <- param$tausq.space*(param$kappa.space^4 * M$M0.space +
                                                  2*param$phi.space * param$kappa.space^2 * M$M1.space +
                                                  M$M2.space)
        precision.direction <- param$tausq.direction*(param$kappa.direction^4 * M$M0.direction +
                                                      2 * param$phi.direction * param$kappa.direction^2 * M$M1.direction +
                                                      M$M2.direction)
        precision <- kronecker(precision.direction, precision.space)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        lrho.space          <- dlnorm(param$rho.space, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
        lsigma.space        <- dexp(param$sigma.space, hyperpar$sigma.spatial.oscillating, log = TRUE)       
        lpphi.space       <- prior.functions$prior.phi_osc(param$phi.space,
                                                           a=hyperpar$a.par.phi.prior.spatial.oscillating,
                                                           b=hyperpar$b.par.phi.prior.spatial.oscillating,
                                                           lg=TRUE)
        lrho.direction    <- dexp(param$rho.direction, hyperpar$rho.directional, log=TRUE)   
        lsigma.direction  <- dexp(param$sigma.direction, hyperpar$sigma.directional, log = TRUE)
        lpphi.direction   <- prior.functions$prior.phi_osc(param$phi.direction,
                                                           a=12,
                                                           b=4,
                                                           lg=TRUE)
        res               <- lpphi.space + lrho.space + lsigma.space + lrho.direction + lsigma.direction + lpphi.direction +
            param$ljac.phi.space + param$ljac.sigma.space + param$ljac.rho.space +
            param$ljac.rho.direction + param$ljac.sigma.direction + param$ljac.phi.direction 
        return(res)
    }
    initial <- function()  return(c(initial.space$theta1, initial.space$theta2, initial.space$theta3,
                                    initial.direction$theta4, initial.direction$theta5, initial.direction$theta6))
    quit    <- function()  return(invisible())
    res     <- do.call(match.arg(cmd), args = list())
    return(res)
}


## oscillatory circular field
## cot <- function(x) 1/tan(x)
## phi_test <- -0.5
## a <- -phi_test
## b <- -(1-phi_test^2)^(1/2)

## x <- Re(complex(real=a,imaginary=b)^(1/2))
## y <- Im(complex(real=a,imaginary=b)^(1/2))

## x <- cos(atan2(y=b, x=a)/2)
## y <- sin(atan2(y=b, x=a)/2)

## a_test <- cos(atan2(y=b, x=a))
## b_test <- sin(atan2(y=b, x=a))

## ((cot(pi*sqrt(complex(real=a, imaginary=b)))/sqrt(complex(real=a, imaginary=b)))-
##  (cot(pi*sqrt(complex(real=a, imaginary=-b)))/sqrt(complex(real=a, imaginary=-b))))/(complex(real=0,imaginary=-4*b))
## 2*((y*sin(2*pi*x)+x*sinh(2*pi*y)))/(4*b*(cosh(2*pi*y)-cos(2*pi*x))) 


'temporal.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        theta.pars      <- tail(theta, 2) 
        rho.time        <- theta.functions$theta.2.rho.time(theta.pars[1])
        kappa.time      <- sqrt(8*(3/2))/rho.time        
        sigma.time      <- theta.functions$theta.2.sigma.time(theta.pars[2])
        ljac.rho.time   <- attr(rho.time, "ljacobian")
        ljac.sigma.time <- attr(sigma.time, "ljacobian")
        z          <- list(sigma.time           = sigma.time,
                           rho.time             = rho.time,
                           kappa.time           = kappa.time,
                           ljac.rho.time        = ljac.rho.time,
                           ljac.sigma.time      = ljac.sigma.time)
        return(z)
    }
    graph <- function() return(M$M2.temporal)
    Q <- function() {
        require(Matrix)
        param             <- interpret.theta()        
        tausq.time        <- 1/(4*(param$sigma.time^2)*(param$kappa.time^3))
        precision         <- tausq.time*(param$kappa.time^4 * M$M0.temporal + 2 * param$kappa.time^2 * M$M1.temporal + M$M2.temporal)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()        
        lrho.time    <- dexp(param$rho.time, hyperpar$rho.prior.rate.time , log=TRUE)
        lsigma.time  <- dexp(param$sigma.time, hyperpar$sigma.prior.rate.time, log = TRUE)
        res              <- lrho.time + lsigma.time + param$ljac.rho.time + + param$ljac.sigma.time
        return(res)
    }
    initial <- function()  return(c(initial.time$theta1, initial.time$theta2))
    quit <- function()  return(invisible())
    res <- do.call(match.arg(cmd), args = list())
    return(res)
}


'temporal.model2' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
    envir <- parent.env(environment())
    interpret.theta <- function() {
        theta.pars      <- tail(theta, 3)
        rho.time        <- theta.functions$theta.2.rho.time(theta.pars[1])
        kappa.time      <- sqrt(8*(3/2))/rho.time        
        sigma.time      <- theta.functions$theta.2.sigma.time(theta.pars[2])
        phi.time         <- theta.functions$theta.2.phi.time(theta.pars[3],
                                                             l=theta.functions$l,
                                                             u=theta.functions$u)
        th.osc <- uniroot(f=function(x) {
            sin(x)-x*sqrt(1-phi.time^2)/acos(phi.time)
        }, interval=c(1e-20, pi-1e-20))$root
        tausq.time      <- 1/(4*(sigma.time^2)*(kappa.time^3)*cos(th.osc/2))
        ljac.rho.time   <- attr(rho.time, "ljacobian")
        ljac.sigma.time <- attr(sigma.time, "ljacobian")
        ljac.phi.time   <- attr(phi.time, "ljacobian")
        z               <- list(sigma.time       = sigma.time,
                                rho.time         = rho.time,
                                kappa.time       = kappa.time,
                                tausq.time       = tausq.time,
                                phi.time         = phi.time,
                                ljac.rho.time    = ljac.rho.time,
                                ljac.sigma.time  = ljac.sigma.time,
                                ljac.phi.time    = ljac.phi.time)
        return(z)
    }
    graph <- function() return(M$M2.temporal)
    Q <- function() {
        require(Matrix)
        param             <- interpret.theta()        
        precision         <- param$tausq.time*(param$kappa.time^4 * M$M0.temporal + 2 * param$kappa.time^2 * param$phi.time*M$M1.temporal + M$M2.temporal)
        return(precision)
    }
    mu <- function() return(numeric(0))
    log.norm.const <- function() return(numeric(0))
    log.prior <- function() {        
        param = interpret.theta()
        lpphi.time   <- prior.functions$prior.phi_osc(param$phi.space,
                                                      a=hyperpar$a.par.phi.time,
                                                      b=hyperpar$b.par.phi.time,
                                                      lg=TRUE)
        lrho.time    <- dexp(param$rho.time, hyperpar$rho.prior.rate.time , log=TRUE)
        lsigma.time  <- dexp(param$sigma.time,hyperpar$sigma.prior.rate.time, log = TRUE)
        res          <- lrho.time + lsigma.time +
            param$ljac.rho.time + + param$ljac.sigma.time  +  param$ljac.phi.time
        return(res)
    }
    initial <- function()  return(c(initial.time$theta1, initial.time$theta2, initial.time$theta3))
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



## 'oscillating.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
##     envir <- parent.env(environment())
##     interpret.theta <- function() {
##         rho     <- theta.functions$theta.2.rho(theta=theta[1L])
##         kappa   <- sqrt(8)/rho
##         sigma   <- theta.functions$theta.2.sigma(theta=theta[2L])
##         ## phi     <- (1-exp(-theta[3L]))/(1+exp(-theta[3L]))
##         phi     <- theta.functions$theta.2.phi(theta=theta[3L], l=theta.functions$l, u=theta.functions$u)
##         sincpth <- sqrt(1-phi^2)/acos(phi)
##         tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
##         z       <- list(tausq = tausq, rho  = rho, phi = phi)
##         return(z)
##     }
##     Q <- function() {
##         require(Matrix)
##         param     <- interpret.theta()
##         kappa     <- sqrt(8)/param$rho
##         precision <- param$tausq*(kappa^4 * M$M0 + 2*param$phi * kappa^2 * M$M1 + M$M2)
##         return(precision)
##     }
##     graph <- function() return(M$M2)
##     mu    <- function() return(numeric(0))
##     log.norm.const <- function() return(numeric(0))
##     log.prior <- function(){        
##         param      = interpret.theta()        
##         rho        <- sqrt(8)/param$kappa
##         sigma      <- sqrt(param$tausq)        
##         ## distribution of hyperparameters         
##         sigma.spatial.oscillating <- hyperpar$sigma.spatial.oscillating
##         lrho.sp    <- dlnorm(rho-5, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
##         ## lrho.sp    <- dexp(rho-5, 1/10, log=TRUE)
##         ## lrho.sp    <- dexp(rho-5, 1/5, log=TRUE)
##         lsigma.sp  <- dexp(sigma, sigma.spatial.oscillating, log = TRUE)
##         lpphi.sp   <- prior.functions$prior.phi_osc(param$phi,
##                                                     a=hyperpar$a.par.phi.prior.spatial.oscillating,
##                                                     b=hyperpar$b.par.phi.prior.spatial.oscillating,
##                                                     l=theta.functions$l, u=theta.functions$u,
##                                                     lg=TRUE)
##         ljac.phi   <- log(attr(theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u), "jacobian"))
##         ljac.sigma <- theta[2L]
##         ljac.rho   <- theta[1L]
##         res        <- lpphi.sp + lrho.sp + lsigma.sp + ljac.phi + ljac.sigma + ljac.rho
##         ## print(paste("theta.phi: ", theta[3L], "theta.rho: ", ljac.rho, "theta.sigma: ", ljac.sigma))
##         return(res)
##     }
##     initial <- function()  return(c(initial.space$theta1, initial.space$theta2, initial.space$theta3))
##     quit <- function()  return(invisible())
##     res <- do.call(match.arg(cmd), args = list())
##     return(res)
## }

## 'oscillating.model.fixed.phi' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
##     envir <- parent.env(environment())
##     interpret.theta <- function() {
##         rho     <- theta.functions$theta.2.rho(theta=theta[1L])
##         kappa   <- sqrt(8)/rho
##         sigma   <- theta.functions$theta.2.sigma(theta=theta[2L])
##         ## phi     <- (1-exp(-theta[3L]))/(1+exp(-theta[3L]))
##         ## phi     <- theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u)
##         sincpth <- sqrt(1-hyperpar$phi^2)/acos(hyperpar$phi)
##         tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
##         z       <- list(tausq = tausq, rho  = rho)
##         return(z)
##     }
##     Q <- function() {
##         require(Matrix)
##         param     <- interpret.theta()
##         kappa     <- sqrt(8)/param$rho
##         precision <- param$tausq*(kappa^4 * M$M0 + 2*hyperpar$phi * kappa^2 * M$M1 + M$M2)
##         return(precision)
##     }
##     graph <- function() return(M$M2)
##     mu <- function() return(numeric(0))
##     log.norm.const <- function() return(numeric(0))
##     log.prior <- function(){        
##         param      = interpret.theta()        
##         rho        <- sqrt(8)/param$kappa
##         sigma      <- sqrt(param$tausq)
##         ## distribution of hyperparameters         
##         sigma.spatial.oscillating <- hyperpar$sigma.spatial.oscillating
##         lrho.sp    <- dlnorm(rho-5, log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)
##         ## lrho.sp    <- dexp(rho-5, 1/10, log=TRUE)
##         ## lrho.sp    <- dexp(rho-5, 1/5, log=TRUE)
##         lsigma.sp  <- dexp(sigma, sigma.spatial.oscillating, log = TRUE)
##         ljac.sigma <- theta[2L]
##         ljac.rho   <- theta[1L]
##         res        <- lrho.sp + lsigma.sp + ljac.sigma + ljac.rho
##         ## print(paste("theta.phi: ", theta[3L], "theta.rho: ", ljac.rho, "theta.sigma: ", ljac.sigma))
##         return(res)
##     }
##     initial <- function()  return(c(initial.space$theta1, initial.space$theta2, initial.space$theta3))
##     quit <- function()  return(invisible())
##     res <- do.call(match.arg(cmd), args = list())
##     return(res)
## }


## 'space.direction.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
##     envir <- parent.env(environment())
##     interpret.theta <- function() {
##         ## space parameters
##         rho.space     <- theta.functions$theta.2.rho(theta=theta[1L])
##         sigma.space   <- theta.functions$theta.2.sigma(theta=theta[2L])
##         phi.space     <- theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u)
##         ## direction parameters
##         rho.direction     <- theta.functions$theta.2.rho.direction(theta=theta[4L])
##         sigma.direction   <- theta.functions$theta.2.sigma(theta=theta[5L])
##         ## 
##         z       <- list(sigma.space = sigma.space, rho.space  = rho.space, phi.space = phi.space,
##                         sigma.direction = sigma.direction, rho.direction=rho.direction)
##         return(z)
##         ## rho     <- theta.functions$theta.2.rho(theta=theta[1L])
##         ## kappa   <- sqrt(8)/rho
##         ## sigma   <- theta.functions$theta.2.sigma(theta=theta[2L])
##         ## ## phi     <- (1-exp(-theta[3L]))/(1+exp(-theta[3L]))
##         ## phi     <- theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u)
##         ## sincpth <- sqrt(1-phi^2)/acos(phi)
##         ## tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
##         ## tausq.space   <- 1/(4*pi*(sigma.space^2)*(kappa.space^2)*sincpth)
##     }
##     graph <- function() return(kronecker(M$M2.direction, M$M2.space))
##     Q <- function() {
##         require(Matrix)
##         param               <- interpret.theta()
##         kappa.space         <- sqrt(8)/param$rho.space
##         kappa.direction     <- sqrt(8*(3/2))/param$rho.direction
##         sincpth             <- sqrt(1-param$phi.space^2)/acos(param$phi.space)
##         tausq.space         <- 1/(4*pi*(param$sigma.space^2)*(kappa.space^2)*sincpth)
##         tausq.direction     <- 1/(4*(param$sigma.direction^2)*(kappa.direction^3))
##         precision.space     <- tausq.space*(kappa.space^4 * M$M0.space + 2*param$phi.space * kappa.space^2 * M$M1.space + M$M2.space)
##         precision.direction <- tausq.direction*(kappa.direction^4 * M$M0.direction + 2 * kappa.direction^2 * M$M1.direction + M$M2.direction)
##         precision           <- kronecker(precision.direction, precision.space)
##         return(precision)
##     }
##     mu <- function() return(numeric(0))
##     log.norm.const <- function() return(numeric(0))
##     log.prior <- function() {        
##         param = interpret.theta()        
##         lrho.space          <- dlnorm(param$rho.space-5,
##                                       log(hyperpar$mu.range.spatial.oscillating), hyperpar$sigma.range.spatial.oscillating, log=TRUE)       
##         lsigma.space        <- dexp(param$sigma.space, hyperpar$sigma.spatial.oscillating, log = TRUE)
##         lpphi.space         <- prior.functions$prior.phi_osc(param$phi.space,
##                                                              a=hyperpar$a.par.phi.prior.spatial.oscillating,
##                                                              b=hyperpar$b.par.phi.prior.spatial.oscillating,
##                                                              l=theta.functions$l, u=theta.functions$u,
##                                                              lg=TRUE)
##         ljac.rho.space        <- theta[1L]
##         ljac.sigma.space      <- theta[2L]
##         ljac.phi.space        <- log(attr(theta.functions$theta.2.phi(theta[3L], l=theta.functions$l, u=theta.functions$u), "jacobian"))
##         ## 
##         lrho.direction        <- dexp(param$rho.direction, hyperpar$rho.directional, log=TRUE)
##         lsigma.direction      <- dexp(param$sigma.direction, hyperpar$sigma.directional, log = TRUE)
##         ljac.rho.direction    <- theta[4L]
##         ljac.sigma.direction  <- theta[5L]
##         ## 
##         res                   <- lpphi.space + lrho.space + lsigma.space + ljac.phi.space + ljac.sigma.space + ljac.rho.space +
##             lrho.direction + lsigma.direction + ljac.rho.direction + ljac.sigma.direction 
##         return(res)
##     }
##     initial <- function()  return(c(initial.space.direction$theta1, initial.space.direction$theta2, initial.space.direction$theta3,
##                                     initial.space.direction$theta4, initial.space.direction$theta5))
##     quit <- function()  return(invisible())
##     res <- do.call(match.arg(cmd), args = list())
##     return(res)
## }

## 'temporal.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"), theta = NULL){
##     envir <- parent.env(environment())
##     interpret.theta <- function() {
##         rho     <- exp(theta[1L])
##         kappa   <- sqrt(8*(3/2))/rho
##         sigma   <- exp(theta[2L])
##         z       <- list(sigma.temporal = sigma, rho.temporal  = rho)
##         return(z)
##     }
##     graph <- function() return(M$M2.temporal)
##     Q <- function() {
##         require(Matrix)
##         param            <- interpret.theta()        
##         kappa.temporal   <- sqrt(8*(3/2))/param$rho.temporal
##         tausq.temporal   <- 1/(4*(param$sigma.temporal^2)*(kappa.temporal^3))
##         precision        <- tausq.temporal*(kappa.temporal^4 * M$M0.temporal + 2 * kappa.temporal^2 * M$M1.temporal + M$M2.temporal)
##         return(precision)
##     }
##     mu <- function() return(numeric(0))
##     log.norm.const <- function() return(numeric(0))
##     log.prior <- function() {        
##         param = interpret.theta()        
##         lrho.temporal    <- dexp(param$rho.temporal, hyperpar$rho.temporal , log=TRUE)
##         lsigma.temporal  <- dexp(param$sigma.temporal, hyperpar$sigma.temporal, log = TRUE)
##         res              <- lrho.temporal + lsigma.temporal
##         return(res)
##     }
##     initial <- function()  return(c(log(22), log(0.68)))
##     quit <- function()  return(invisible())
##     res <- do.call(match.arg(cmd), args = list())
##     return(res)
## }

## space <- function(time){
##     f <- function(time.scalar){
##         dat <<- data$Ypos
##         if(time.scalar <= dat$time[1]) return(dat$coords[1,])
##         if(time.scalar >= dat$time.lead[length(dat$time.lead)]) return(dat$coords[nrow(dat),])
##         else{
##             wh.start <- max(which(dat$time <= time.scalar))
##             bary <- (dat$time.lead[wh.start] - time.scalar)/(dat$time.lead[wh.start]- dat$time[wh.start])
##             return(bary*dat$coords[wh.start,] + (1-bary)*dat$coords.lead[wh.start,])
##         }
##     }
##     t((Vectorize(f, vectorize.args="time.scalar"))(time))
## }

## direction <- function(time){
##     f <- function(time.scalar){
##         dat <<- data$Ypos
##         wh.start <- max(which(dat$time <= time.scalar))
##         bary <- (dat$time.lead[wh.start] - time.scalar)/(dat$time.lead[wh.start]- dat$time[wh.start])
##         bary*dat$hd[wh.start] + (1-bary)*dat$hd.lead[wh.start,]
##     }
##     (Vectorize(f, vectorize.args="time.scalar"))(time)
## }


## ## oscillating.rgeneric <- inla.rgeneric.define(oscillating.model) #arguments that need to be passed are M0, M1 and M2 

## 'circular1D.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
##                                theta = NULL) {
##     envir <- parent.env(environment())
##     interpret.theta <- function() {
##         kappa    <- exp(theta[1L])
##         sigma   <- exp(theta[2L])        
##         tausq   <- (pi + sinh(2 * pi * kappa)/(2*kappa))/((sinh(pi*kappa)^2)*(sigma^2)*(2*kappa.dir)^2)
##         z       <- list(kappa  = kappa, tausq = tausq)
##         return(z)
##     }
##     graph <- function() return(M2)
##     Q <- function() {
##         require(Matrix)
##         param <- interpret.theta()
##         precision <- param$tausq*(param$kappa^4 * M0 + 2 * param$kappa^2 * M1 + M2)
##         return(precision)
##     }
##     mu <- function() return(numeric(0))
##     log.norm.const <- function() return(numeric(0))
##     log.prior <- function() {
##         param = interpret.theta()
##         sigma.hd <- 1/sqrt(param$tausq)
##         rho.hd      <- sqrt(8)/param$kappa
##         lrho.hd    <- dexp(rho.hd, 1/pi, log=TRUE)        
##         lsigma.hd  <- dexp(sigma.hd, 1, log = TRUE)
##         if(FALSE){
##             rpi.div.r0 <- 2*(pi*param$kappa*cosh(pi*param$kappa)+sinh(pi*param$kappa))/(2*pi*param$kappa + sinh(2*pi*param$kappa))
##             lrpi.div.r0.hd <- dbeta(rbi.div.r0, alpha=1, beta=10)
##         }
##         res <- lrho.hd + lsigma.hd
##         return(res)
##     }
##     initial <- function()  return(c(0, 0))
##     quit <- function()  return(invisible())
##     res <- do.call(match.arg(cmd), args = list())
##     return(res)
## }

## ## circular1D <- inla.rgeneric.define(circular1D.model, M0=NULL, M1=NULL, M2=NULL) 


## kappa.space         <- sqrt(8)/param$rho.space
## kappa.direction     <- sqrt(8*(3/2))/param$rho.direction
## tausq.space         <- 1/(4*pi*(param$sigma.space^2)*(kappa.space^2)*sincpth)
## sincpth             <- sqrt(1-param$phi.space^2)/acos(param$phi.space)
## tausq.direction     <- 1/(4*(param$sigma.direction^2)*(param$kappa.direction^3))
## sigmaLN  <- 3
## murho    <- 30
## sigma.space               <- param$sigma.space
## phi.space                 <- param$phi.space
## sigma.direction           <- param$sigma.direction
## sigmaLN                   <- hyperpar$sigma.range.spatial.oscillating
## murho                     <- hyperpar$mu.range.spatial.oscillating
## sigma.spatial.oscillating <- hyperpar$sigma.spatial.oscillating
## lrho.space        <- dlnorm(param$rho.space, log(murho), sigmaLN, log=TRUE)
## lsigma.space      <- dexp(sigma.space, 1/2, log = TRUE)
