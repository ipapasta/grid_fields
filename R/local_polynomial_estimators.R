## ----------------------------------
## Local linear regression estimator
## ----------------------------------
## objective.path(par, s1, s2, si1)
## lambda(s)
objective.path <- function(par, s1, s2, si1, si2, L, z.mid, h){
    beta0 <- par[1]
    beta1 <- par[2]
    beta2 <- par[3]
    int.ell <- -exp(beta0) * sum(L*exp(beta1*(z.mid[,1]-s1)+beta2*(z.mid[,2]-s2))*(dnorm((z.mid[,1]-s1)/h, 0, 1)/h)*(dnorm((z.mid[,2]-s2)/h, 0, 1)/h))
    ## if(is.nan(int.ell)) int.ell <- 0
    sum.ell <- sum((beta0+beta1*(si1-s1)+beta2*(si2-s2) )*(dnorm((si1-s1)/h, 0, 1)/h)*dnorm((si2-s2)/h, 0, 1)/h)
    o <- -(int.ell + sum.ell)
    return(o)
}

local.linear.nhpp.single.location.path <- function(init, s1, s2, si1, si2, x1, x2, y1, y2, L, z.mid, h){
    res <- optim(par=c(init,0,0), fn=objective.path, s1=s1, s2=s2, si1=si1, si2=si2, L=L, z.mid=z.mid, h=h,
                 control=list(maxit=1000))
    res <- optim(par=res$par, fn=objective.path, s1=s1, s2=s2, si1=si1, si2=si2, L=L, z.mid=z.mid, h=h,
                 method="BFGS")
    return(res$par)
}
local.linear.nhpp.path <- Vectorize(local.linear.nhpp.single.location.path, vectorize.args=c("init", "s1", "s2"))


local.linear.nhpp.path
local.constant.nhpp.single.location.path <- function(s1, s2, si1, si2, L, z.mid, h){
    numer <- sum((dnorm((si1-s1)/h, 0, 1)/h)*dnorm((si2-s2)/h, 0, 1)/h)
    denom <- sum(L*dnorm((z.mid[,1]-s1)/h)*dnorm((z.mid[,2]-s2)/h)/(h^2))
    log(numer/denom)
}
local.constant.nhpp.path <- Vectorize(local.constant.nhpp.single.location.path, vectorize.args=c("s1", "s2"))

local.constant.nhpp.single.location.path.sargolini <- function(s1, s2, si1, si2, Delta.t, z.mid, h){
    numer <- sum((dnorm((si1-s1)/h, 0, 1)/h)*dnorm((si2-s2)/h, 0, 1)/h)
    denom <- sum(Delta.t*dnorm((z.mid[,1]-s1)/h)*dnorm((z.mid[,2]-s2)/h)/(h^2))
    log(numer/denom)
}
local.constant.nhpp.path.sargolini <- Vectorize(local.constant.nhpp.single.location.path.sargolini, vectorize.args=c("s1", "s2"))

## ------------------------------------------
## Simulation of NHPP with Lewis' method, Ogata's method (marked point process X)
## ------------------------------------------
## simulation performed only for model with spatial effect
##



sim.nhpp <- function(data, lambda_star, h){
    df         <- data$Ypos %>% mutate(dist = map(Li, function(x) sum(x)))
    df.firings <- data$Y
    absGamma   <- sum(unlist(df$dist))    
    n.hpp      <- floor(absGamma*lambda_star)
    p.hpp      <- cumsum(rexp(n.hpp, lambda_star))
    p.hpp      <- p.hpp[p.hpp < absGamma & p.hpp]
    Ysim.hpp   <- t(sapply(p.hpp, function(d, d.vec=c(cumsum(unlist(df$dist))),
                                           coords=df$coords,
                                           lead=df$coords.lead){
        ind         <- (max(which(d.vec <= d))) 
        start.point <- coords[ind,]
        dist.diff   <- (d - d.vec[ind])
        unit.vector <- (lead[ind,]-coords[ind,])/
            sqrt(sum((lead[ind,]-coords[ind,])^2))
        start.point + dist.diff*unit.vector
    }))
    ## compute the intensity function at the simulated HPP points
    intensity.const <- exp(local.constant.nhpp.path(s1=Ysim.hpp[,1], s2=Ysim.hpp[,2],
                                                    si1=df.firings$position_x,
                                                    si2=df.firings$position_y,
                                                    x1=-2,
                                                    x2= 103,
                                                    y1=-2,
                                                    y2=103,
                                                    L = unlist(df$dist),
                                                    z.mid = (df$coords + df$coords.lead)/2,
                                                    h=h))
    U                <- runif(nrow(Ysim.hpp))
    thinning.vector  <- U <= intensity.const/lambda_star
    Ysim <- Ysim.hpp[which(thinning.vector), ]
    return(Ysim)
}

## checking integral of estimated intensity matches number of observed events
if(FALSE){
    df <- data$Ypos %>% mutate(dist = map(Li, function(x) sum(x)))
    coords.mid <- (df$coords + df$coords.lead)/2
    lambda.coords.mid <- exp(local.constant.nhpp.2d(s1=coords.mid[,1], s2=coords.mid[,2],
                                                  si1=data$Y$position_x,
                                                  si2=data$Y$position_y,
                                                  x1=-2,
                                                  x2= 103,
                                                  y1=-2,
                                                  y2 =103,
                                                  h=1))

    o <- sum(unlist(df$dist) * lambda.coords.mid) #this gives wrong estimates.g 
}




objective.linear.1d <- function(par, time, dat, T, h){
    beta0 <- par[1]
    beta1 <- par[2]
    int.ell <- -exp(beta0)*exp((h^2 * beta1^2)/2)*
        (pnorm((time+(h^2) * beta1)/h)-pnorm((time-T+(h^2) * beta1)/h))
    sum.ell <- sum((beta0+beta1*(dat-time))*(dnorm((dat-time)/h, 0, 1)/h))
    o <- -(sum.ell + int.ell)
    return(o)
}
local.linear.nhpp.single.location.1d <- function(time, dat, T, h){
    init <- c(0,0)
    res <- optim(par=init, fn=objective.linear.1d, time=time,
                 dat=dat, T=T, h=h, control=list(maxit=10000))
    o <- res$par
    attr(o, "conv") <- res$convergence
    return(o)
}
local.linear.nhpp.1d <- Vectorize(local.linear.nhpp.single.location.1d, vectorize.args=c("time"))


local.constant.nhpp.single.location.1d <- function(time, dat, T, h){
    numer <- sum((dnorm((dat-time)/h, 0, 1)/h))
    denom <- (pnorm((T-time)/h)-pnorm((-time)/h))
    log(numer/denom)
    ## numer
}
local.constant.nhpp.1d <- Vectorize(local.constant.nhpp.single.location.1d, vectorize.args=c("time"))


## do not run 
if(FALSE){
    objective.2d <- function(par, s1, s2, si1, si2, x, y, h){
        beta0 <- par[1]
        beta1 <- par[2]
        beta2 <- par[3]
        int.ell <- -exp(beta0 + (h^2/2)*(beta1^2+beta2^2)) * (pnorm((s1 + beta1*h^2)/h)-pnorm((s1 - x + beta1*h^2)/h)) *
            (pnorm((s2 + beta2*h^2)/h)-pnorm((s2 - y + beta2*h^2)/h))
        if(is.nan(int.ell)) int.ell <- 0
        sum.ell <- sum((beta0+beta1*(si1-s1)+beta2*(si2-s2) )*(dnorm((si1-s1)/h, 0, 1)/h)*dnorm((si2-s2)/h, 0, 1)/h)
        o <- -(int.ell + sum.ell)
        return(o)
    }
    local.linear.nhpp.single.location.2d <- function(s1, s2, s.dat1, s.dat2, h){
        while(TRUE){
            init <- c(runif(1, -.5, .5), runif(1, -.1, .1), runif(1, -.1, .1))
            ## init <- c(-1,-1,-1)
            res <- optim(par=init, fn=objective.2d, s1=s1, s2=s2, si1=s.dat1, si2=s.dat2, x=x, y=y, h=h,
                         control=list(maxit=10000))
            res <- optim(par=res$par, fn=objective.2d, s1=s1, s2=s2, si1=s.dat1, si2=s.dat2, x=x, y=y, h=h,
                         method="BFGS")
            if(res$convergence==0) break
        }
        return(res$par)
    }
    local.linear.nhpp.2d <- Vectorize(local.linear.nhpp.single.location.2d, vectorize.args=c("s1", "s2"))

    local.constant.nhpp.single.location.2d <- function(s1, s2, si1, si2, x1, x2, y1, y2, h){
        numer <- sum((dnorm((si1-s1)/h, 0, 1)/h)*dnorm((si2-s2)/h, 0, 1)/h)
        denom <- (pnorm((s1-x1)/h)-pnorm((s1-x2)/h)) *
            (pnorm((s2-y1)/h)-pnorm((s2-y2)/h))
        log(numer/denom)
    }
    local.constant.nhpp.2d <- Vectorize(local.constant.nhpp.single.location.2d, vectorize.args=c("s1", "s2"))

    ## Local polynomial regression estimators
    local.constant.nhpp.single.location.2d <- function(s1, s2, si1, si2, x, y, h){
        numer <- sum((dnorm((si1-s1)/h, 0, 1)/h)*dnorm((si2-s2)/h, 0, 1)/h)
        denom <- (pnorm(s1/h)-pnorm((s1 - x)/h)) *
            (pnorm((s2 )/h)-pnorm((s2 - y)/h))
        log(numer/denom)
    }
    local.constant.nhpp.2d <- Vectorize(local.constant.nhpp.single.location.2d, vectorize.args=c("s1", "s2"))
}

