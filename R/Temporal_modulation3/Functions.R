filled.contour2 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2]) * par("csi") * 2.54
  par(las = las)
  mar <- mar.orig
  plot.new()
  par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}


ldmvnorm <- function(x, m, Q){
    ## ----------------------------------------------------------------
    ## input: x (evaluation point), 
    ##        m (mean vector), 
    ##        Q (precision matrix/inverse covariance Sigma^{-1})
    ## output: log density of multivariate normal distribution
    ## ----------------------------------------------------------------
    k      <- length(x)
    xminmu <- Matrix(x-m,ncol=1)
    L      <- (Matrix::Cholesky(Q, perm=TRUE) %>% (Matrix::expand))$L
    halflogdetQ   <- sum(log(diag(L)))     #!!!
    ## print(paste("halflogdetQ is:", halflogdetQ))
    out    <- halflogdetQ - (k/2)*log(2*pi) - .5*sum((xminmu)* (Q %*% xminmu)) # remove matrix multiplication
    out <- as.numeric(out)
    return(out)
}

## theta <- xpar.theta  <- c(3, -1/2, -3.8, log(7), log(10))
## theta.X <- theta[1:3]
## theta.Z <- theta[4:5]    
## QX     <- osc.precision(theta=theta.X, mesh=mesh)
## QZ     <- temp.precision(theta=theta.Z, mesh=mesh1d) 

## for(i in 1:1000){
##     Sys.sleep(1)
##     ## par(mfrow=c(1,2))
##     ## plot(i1)
##     ## plot(i2)
##     ## i1 <- image(QZ[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2)
##     ## i2 <- image(QX[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2)    
##     plot(image(QZ[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2))
##     Sys.sleep(2)
##     plot(image(QX[((i-1)*20 + 1):(i*20),((i-1)*20 + 1):(i*20)], asp=1, lwd=2), add=TRUE)
##     Sys.sleep(4)
## }

## x<-Z
## m<-matrix(rep(0, length(Z)), ncol=1)
## Q<-QZ

interpolate <- function(x, y, z){
    interp <- (x) + (cumsum((z)/sum((z))))*((y) - (x))
    Delta <- diff(c(x, interp))
    attr(Delta, "data") <- c(x, interp[-length(interp)]) # note that
                                                         # this way
                                                         # the last
                                                         # value from
                                                         # the
                                                         # combined
                                                         # data vector
                                                         # will be
                                                         # missing so
                                                         # needs to be
                                                         # appended later
    return(Delta)
}

interpolate2 <- function(x, y, z){
    interp <- matrix(nrow=length(z), ncol=2)
    for(i in 1:length(z)){
        interp[i,] <- (x) + (sum(z[1:i])/sum((z)))*((y) - (x))
    }
    o <- rbind(x, interp)
    oo <- cbind(head(o, -1), tail(o, -1))
    rownames(oo) <- NULL
    colnames(oo) <- NULL
    return(oo)
}

split.lines <- function(mesh, sp, ep, filter.zero.length = TRUE, tol= 1e-8) {
                                        # locations for splitting
        loc = as.matrix(rbind(sp,ep))
        idx = 1:dim(sp)[1]
                                        # Filter out segments not on the mesh
        t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
        t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
                                        # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
        sp = matrix(sp[!((t1==0) | (t2==0)),], ncol=2)
        ep = matrix(ep[!((t1==0) | (t2==0)),], ncol=2)
        idx = idx[!((t1==0) | (t2==0))]
        loc = as.matrix(rbind(sp,ep))
                                        # Split them segments into parts
        if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
        np = dim(sp)[1]
        sp.idx = t(rbind(1:np,np+1:np))
        splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
                                        #plot(data$mesh)
                                        #points(loc)
                                        #points(splt$split.loc,col="blue)
        sp = matrix(splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]], ncol=2) # Start point of new segments
        ep = matrix(splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]], ncol=2) # End points of new segments
        idx = idx[splt$split.idx[,1]]
        origin = splt$split.origin
                                        # Filter out zero length segments
        if ( filter.zero.length ) {
            sl = apply((ep-sp)^2,MARGIN=1,sum)
            sp = sp[!(sl<tol^2),]
            ep = ep[!(sl<tol^2),]
            origin = origin[!(sl<tol^2)]
            idx = idx[!(sl<tol^2)]
        }
        return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
    }




split.lines <- function(mesh, sp, ep, filter.zero.length = TRUE, tol= 1e-8, return.filter.index=TRUE) {
    ## locations for splitting
    loc = as.matrix(rbind(sp,ep))
    idx = 1:dim(sp)[1]
    ## Filter out segments not on the mesh
    t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
    t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
    ## if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
    sp = matrix(sp[!((t1==0) | (t2==0)),], ncol=2)
    ep = matrix(ep[!((t1==0) | (t2==0)),], ncol=2)
    idx = idx[!((t1==0) | (t2==0))]
    loc = as.matrix(rbind(sp,ep))
    ## Split them segments into parts
    if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
    np = dim(sp)[1]
    sp.idx = t(rbind(1:np,np+1:np))
    splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
    ##plot(data$mesh)
    ##points(loc)
    ##points(splt$split.loc,col="blue)
    sp = matrix(splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]], ncol=2) # Start point of new segments
    ep = matrix(splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]], ncol=2) # End points of new segments
    idx = idx[splt$split.idx[,1]]
    origin = splt$split.origin
    sl = apply((ep-sp)^2,MARGIN=1,sum)
    filter.index <- !(sl < tol^2)
    ## Filter out zero length segments
    if ( filter.zero.length ) {
        ## sl = apply((ep-sp)^2,MARGIN=1,sum)
        sp = sp[!(sl<tol^2),]
        ep = ep[!(sl<tol^2),]
        origin = origin[!(sl<tol^2)]
        idx = idx[!(sl<tol^2)]
        filter.index = filter.index[!(sl<tol^2)]
        ## drop filter.index
        return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
    }
    return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc, filter.index=filter.index))
}


split.arcs <- function(hd, hd.lead){
    ## print(paste("hd: ",hd, "hd.lead: ",hd))
    hd.int <- which(apply(intervals, 1, function(x) (x[1] <= hd) & (hd <= x[2])))
    hd.lead.int <- which(apply(intervals, 1, function(x) (x[1] <= hd.lead) & (hd.lead <= x[2])))
    if(hd.int==hd.lead.int){
        return(matrix(c(hd, hd.lead),nrow=1))
    }else{
        if(hd < hd.lead){
            if(hd.lead.int-hd.int == 1){
                rbind(c(hd, intervals[hd.int,2]),
                      c(intervals[hd.lead.int,1], hd.lead))
            }else{
                return(rbind(c(hd, intervals[(hd.int+1),1]),
                             intervals[(hd.int+1):(hd.lead.int-1),],
                             c(intervals[(hd.lead.int-1),2], hd.lead)))
            }
        }else{
            if(hd.int-hd.lead.int == 1){
                rbind(c(hd, intervals[hd.int,1]),
                      c(intervals[hd.lead.int,2], hd.lead))
            }else{
            return(rbind(c(hd, intervals[(hd.int-1),2]),
                         intervals[(hd.int-1):(hd.lead.int+1),2:1],
                         c(intervals[(hd.lead.int+1),1], hd.lead)))
            }
            ## rbind(c(hd, intervals[(hd.int-1),1]),
            ##       apply((intervals[(hd.int-1):(hd.lead.int+1),]), 1, function(x)x),
            ##       c(intervals[(hd.lead.int-1),2], hd.lead))
        }
    }
}




## 
## Local linear regression estimator 
##

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
## Simulation of NHPP with Lewis' method
## ------------------------------------------
## simulation performed only for model with spatial effect
##

sim.nhpp <- function(data, lambda_star, h){
    df <- data$Ypos %>% mutate(dist = map(Li, function(x) sum(x)))
    df.firings <- data$Y
    absGamma <- sum(unlist(df$dist))    
    n.hpp <- floor(absGamma*lambda_star)
    p.hpp <- cumsum(rexp(n.hpp, lambda_star))
    p.hpp <- p.hpp[p.hpp < absGamma & p.hpp]
    Ysim.hpp <- t(sapply(p.hpp, function(d, d.vec=c(cumsum(unlist(df$dist))),
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
    o <- -(int.ell + sum.ell)
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
