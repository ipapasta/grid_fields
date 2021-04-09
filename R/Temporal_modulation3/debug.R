library(reprex)
reprex({    
    library(INLA)
    loc  <- cbind(runif(100,0,100), runif(100,0,100))
    mesh <- inla.mesh.create(loc = loc, refine = list(max.edge = 10))
    fem  <- inla.mesh.fem(mesh, order = 2)
    M0 = fem$c0
    M1 = fem$g1
    M2 = fem$g2
    B0 = matrix(c(0,1,0,0), nrow=1)
    B1 = matrix(c(0,0,1,0), nrow=1)
    B2 = matrix(c(0,0,0,1), nrow=1)
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
    spde = inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2, 
                              B0 = B0, B1 = B1, B2 = B1,
                              theta.mu = c(0,0,0), 
                              theta.Q = diag(c(1,1,1)),
                              transform = "log")

    n = 100
    field.fcn = function(loc) (10*cos(2*pi*2*(loc[,1]+loc[,2])))
    loc = matrix(runif(n*2),n,2)
    idx.y = rep(1:n,2)
    y = field.fcn(loc[idx.y,]) + rnorm(length(idx.y))
    mesh = inla.mesh.create(loc, refine=list(max.edge=0.05))
    data = list(y=y, field=mesh$idx$loc[idx.y])
    formula1 = y ~ -1 + f(field, model=spde)
    ## result1 = inla(formula1, data=data, family="normal", verbose=TRUE)
    oscillating.rgeneric <- inla.rgeneric.define(oscillating.model, M = list(M0=M0, M1=M1, M2=M2)) #arguments that need to be passed are M0, M1 and M2
    formula2 = y ~ -1 + f(field, model = oscillating.rgeneric) 
    result2 = inla(formula2, data=data, family="normal", verbose=TRUE)
    ## cmp <- coordinates ~ mySmooth(coordinates, model = ,
    ##                               mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
},style=TRUE,session_info=TRUE, advertise = FALSE,
html_preview    = TRUE,
comment         = "#;-)",
tidyverse_quiet = FALSE,
std_out_err     = TRUE
)




    ## 
    ## tau    <- 1
    ## kappa  <- sqrt(8)/10
    ## phi    <- -.9
    ## Q <- (tau^2)*(kappa^4*M0 + 2*(kappa^2)*phi*M1 + M2)
    ## x <- inla.qsample(1, Q)
