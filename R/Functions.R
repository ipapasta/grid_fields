## ---------------------------------------
## evaluating multivariate joint density
## ---------------------------------------
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


##
## Input: 
##      - theta1: correlation length parameter. Note that
##                this parameter does not have the same interpretation
##                as in the usual Matern covariance family due to the oscillation
##                of the process. However, it can still be used to describe
##                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
##                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
##      - theta2: scale parameter that controls the variance of the field.
##                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
##                standard deviation of the process but is related to.
##      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
##                -1 < phi < 0: oscillating 
##                 0 < phi < 1: overdampened oscillating 
## Output:
##      - precision matrix for the weights at the mesh vertices 
##        
## 

osc.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    theta3 <- theta[3]
    rho     <- exp(theta1)
    kappa   <- sqrt(8)/rho
    sigma   <- exp(theta2)
    phi     <- (1-exp(-theta3))/(1+exp(-theta3))
    sincpth <- sqrt(1-phi^2)/acos(phi)
    tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
    ## ---------------------------------------
    o       <- inla.mesh.fem(mesh, order = 2)
    ## 
    tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*phi*(o$g1)) + o$g2)
    ## 
}



#' 

#' Input: 
#'      - theta1: correlation length parameter. Note that
#'                this parameter does not have the same interpretation
#'                as in the usual Matern covariance family due to the oscillation
#'                of the process. However, it can still be used to describe
#'                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
#'                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
#'      - theta2: scale parameter that controls the variance of the field.
#'                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
#'                standard deviation of the process but is related to.
#'      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
#'                -1 < phi < 0: oscillating 
#'                 0 < phi < 1: overdampened oscillating 
#' Output:
#'      - precision matrix for the weights at the mesh vertices 
#'        
#' 

temp.precision <- function(theta, mesh, o=2){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    sigma   <- exp(theta2)
    tausq   <- 1/(4*(sigma^2)*(kappa^3))
    ## ---------------------------------------
    o       <- inla.mesh.fem(mesh=mesh, order = o)
    ## 
    tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*(o$g1)) + o$g2)
    ## 
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


split.arcs <- function(hd, hd.lead, mesh.hd){
    ## counter <<- counter + 1
    nodes       <- c(mesh.hd$loc, 2*pi)
    intervals   <- head(cbind(nodes, lead(nodes)), -1)
    ## print(paste("hd: ",hd, "hd.lead: ",hd))
    hd.int <- which(apply(intervals, 1, function(x) (x[1] < hd) & (hd <= x[2])))
    hd.lead.int <- which(apply(intervals, 1, function(x) (x[1] < hd.lead) & (hd.lead <= x[2])))
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


split.segments.wrapper.function <- function(X, mesh, mesh.hd){
    Ypos.tmp <- data.frame(
        hd       = X$hd,
        time     = X$synced_time,
        index.CV = X$index.CV, 
        coords   = I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
        mutate(coords.lead = lead(coords)) %>%
        mutate(time.lead   = lead(X$synced_time)) %>%
        mutate(hd.lead     = lead(X$hd)) %>%
        head(-1)
    Ypos.tmp <- Ypos.tmp %>% mutate(HD.split = map2(hd, hd.lead, function(x, y) split.arcs(x,y, mesh.hd=mesh.hd)),
                                    L.arcs = lapply(HD.split,
                                                    function(x) apply(x, 1,
                                                                      function(y) abs(y[2]-y[1]))),
                                    time.split = pmap(list(time, time.lead, L.arcs), function(x,y,z){
                                        o <- interpolate(x,y,z)
                                        oo <- c(attr(o, "data"), y)
                                        ooo <- head(cbind(oo, lead(oo)), -1)
                                        colnames(ooo) <- NULL
                                        return(ooo)
                                    }),
                                    coords.split=pmap(list(coords, coords.lead, L.arcs), function(x,y,z){
                                        interpolate2(x, y, z)
                                    }),
                                    new.time = lapply(time.split, function(x) x[,1, drop=FALSE]),
                                    new.time.lead= lapply(time.split, function(x) x[,2, drop=FALSE]),
                                    new.hd = lapply(HD.split, function(x) x[,1, drop=FALSE]),
                                    new.hd.lead = lapply(HD.split, function(x) x[,2, drop=FALSE]),
                                    new.coords = lapply(coords.split, function(x) x[,1:2, drop=FALSE]),
                                    new.coords.lead = lapply(coords.split, function(x) x[,3:4, drop=FALSE])
                                    )
    ## Ypos.tmp <- Ypos.tmp %>% dplyr::select(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
    ##     unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))
    ## names(Ypos.tmp) <- c("time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")    
    Ypos.tmp <- Ypos.tmp %>% dplyr::select(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead, index.CV)%>%
        unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead, index.CV))
    names(Ypos.tmp) <- c("time", "time.lead", "hd", "hd.lead", "coords", "coords.lead", "index.CV")    
    line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
                                 filter.zero.length=FALSE,
                                 ep=Ypos.tmp$coords.lead, tol=.0)
    df <- data.frame(origin=line.segments$split.origin,
                     filter.index=line.segments$filter.index,
                     sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                     ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
        group_by(origin) %>%
        summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
        mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
        mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 
    ## 
    ## 
    ## attribute named _data_ stores length of line segments, time differences and arclengths
    Ypos <- inner_join(Ypos.tmp %>%
                       mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
        mutate(Li = map2(sp, ep,
                         function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
        mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
        mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 
    filter.index  <- do.call("c", Ypos$filter.index)
    z <- list()
    z$Ypos <- Ypos
    z$filter.index  <- filter.index
    z$line.segments <- line.segments
    return(z)
}



df.prism.M0.wrapper <- function(Aosc.indices, dGamma, T.data, HD.data, coords.trap) {
    df.prism.M0 <- Aosc.indices %>% group_by(tk) %>% nest() %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma=c(dGamma,0),
            dGamma.lead = lead(dGamma),
            dGamma.lag = lag(dGamma),
            val.M0 = pmap(list(data, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.o[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M0=ooo)
                oooo
            })) %>% 
        dplyr::select(-c("data"))
    return(df.prism.M0)
}

df.prism.M0.wrapper_CV <- function(Aosc.indices, dGamma, T.data, HD.data, coords.trap) {
    df.prism.M0 <- Aosc.indices %>% group_by(tk) %>% nest() %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            index.CV      = index.CV,
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma=c(dGamma,0),
            dGamma.lead = lead(dGamma),
            dGamma.lag = lag(dGamma),
            val.M0 = pmap(list(data, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.o[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M0=ooo)
                oooo
            })) %>% 
        dplyr::select(-c("data"))
    return(df.prism.M0)
}

df.prism.M1.M2.wrapper <- function(At.indices, A.indices, T.data, dGamma, HD.data, coords.trap){
    df.prism.M1_M2 <- full_join(At.indices %>% group_by(tk) %>% nest(),
                                A.indices %>% group_by(tk) %>% nest(), by="tk") %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma=c(dGamma,0),
            dGamma.lead = lead(dGamma),
            dGamma.lag = lag(dGamma),
            val.M1 = pmap(list(data.y, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.ot[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M1=ooo)
                oooo
            }),
            val.M2 = pmap(list(data.x, data.y, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                    x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            })) %>%
        dplyr::select(-c("data.x", "data.y")) 
}

## 
df.prism.M1.M2.wrapper2 <- function(At.indices, Aosc.indices, A.indices, T.data, dGamma, HD.data, coords.trap){
    df.prism.M1_M2 <- full_join(full_join(At.indices %>% group_by(tk) %>% nest(),
                                          A.indices %>% group_by(tk) %>% nest(), by="tk"),
                                Aosc.indices %>% group_by(tk) %>% nest(), by="tk") %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma        = c(dGamma,0),
            dGamma.lead   = lead(dGamma),
            dGamma.lag    = lag(dGamma),
            val.M1 = pmap(list(data.y, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.ot[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M1=ooo)
                oooo
            }),
            val.M2.space.time = pmap(list(data.x, data, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                    x$psi.t[oo[k,1]] * y$psi.o[oo[k,2]]}))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            }),
            val.M2.space.direction.time = pmap(list(data.x, data.y, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                    x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            })) %>%
        dplyr::select(-c("data.x", "data.y", "data")) 
}

df.prism.M1.M2.wrapper2_CV <- function(At.indices, Aosc.indices, A.indices, T.data, dGamma, HD.data, coords.trap){
    df.prism.M1_M2 <- full_join(full_join(At.indices %>% group_by(tk) %>% nest(),
                                          A.indices %>% group_by(tk) %>% nest(), by="tk"),
                                Aosc.indices %>% group_by(tk) %>% nest(), by="tk") %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            index.CV      = index.CV,
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma        = c(dGamma,0),
            dGamma.lead   = lead(dGamma),
            dGamma.lag    = lag(dGamma),
            val.M1 = pmap(list(data.y, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.ot[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M1=ooo)
                oooo
            }),
            val.M2.space.time = pmap(list(data.x, data, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                    x$psi.t[oo[k,1]] * y$psi.o[oo[k,2]]}))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            }),
            val.M2.space.direction.time = pmap(list(data.x, data.y, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                    x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            })) %>%
        dplyr::select(-c("data.x", "data.y", "data")) 
}



## --------------
## link functions
## --------------

## oscillating parameter
theta.2.phi   <- function(theta, l=NULL, u=NULL) {
    if(is.null(l)) l <- get("l", envir = .GlobalEnv)
    if(is.null(u)) u <- get("u", envir = .GlobalEnv)
    res <- l + (u-l)*pnorm(theta, 0, 1)## (1/(1+exp(-theta)))
    attr(res, "ljacobian") <- log(u-l) + dnorm(theta, 0, 1, log=TRUE)
    return(res)
}
theta.2.sigma <- function(theta){
    res <- exp(theta)
    attr(res, "ljacobian") <- theta
    return(res)
}
##
theta.2.rho <- function(theta){
    res <- exp(theta) + 5
    attr(res, "ljacobian") <- theta
    return(res)
}
theta.2.kappa.1d <- function(theta){
    sqrt(8*(3/2))/exp(theta)       
}

##
theta.2.rho.direction <- function(theta){
    res <- exp(theta)
    attr(res, "ljacobian") <- theta
    return(res)
}
theta.2.sigma.direction <- function(theta){
    res <- exp(theta)
    attr(res, "ljacobian") <- theta
    return(res)
}
theta.2.rho.time <- function(theta){
    res <- exp(theta)
    attr(res, "ljacobian") <- theta
    return(res)
}
theta.2.sigma.time <- function(theta){
    res <- exp(theta)
    attr(res, "ljacobian") <- theta
    return(res)
}
theta.2.phi.time   <- function(theta, l=NULL, u=NULL) {
    if(is.null(l)) l <- get("l", envir = .GlobalEnv)
    if(is.null(u)) u <- get("u", envir = .GlobalEnv)
    res <- l + (u-l)*pnorm(theta, 0, 1)## (1/(1+exp(-theta)))
    attr(res, "ljacobian") <- log(u-l) + dnorm(theta, 0, 1, log=TRUE)
    return(res)
}


## phi.seq <- seq(-.99,.99,len=100)
## tmp     <- NULL
## for(i in 1:length(phi.seq))
##     {
##         tmp[i] <- uniroot(f=function(x) {
##             sin(x)-(x*sqrt(1-phi.seq[i]^2)/acos(phi.seq[i]))
##         }, interval=c(1e-20, pi-1e-20))$root
##     }

## helper functions for computing the
## posterior distribution of gridness score
## input is a fitted model with inlabru (e.g. fit.space and fit.space.direction)
## object returns attributes for summaries of the posterior
posterior.spatial.gridness.score <- function(inlabru.fitted.object, theta.mapping){    
    marg                  <- inla.tmarginal(theta.mapping, inlabru.fitted.object$marginals.hyperpar[["Theta3 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.directional.gridness.score <- function(inlabru.fitted.object, theta.mapping){    
    marg                  <- inla.tmarginal(theta.mapping, inlabru.fitted.object$marginals.hyperpar[["Theta6 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}

posterior.spatial.standard.deviation <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta2 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.spatial.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta1 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}

posterior.directional.standard.deviation <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.sigma.direction, inlabru.fitted.object$marginals.hyperpar[["Theta5 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.directional.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho.direction, inlabru.fitted.object$marginals.hyperpar[["Theta4 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.directional.kappa <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.kappa.1d, inlabru.fitted.object$marginals.hyperpar[["Theta4 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}


posterior.temporal.standard.deviation <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta6 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.temporal.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta7 for spde4"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}


## priors of hyperparameters
prior.phi_osc <- function(phi, a, b, l=(-0.998), u=1, lg=TRUE){
    if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
    if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
}
## prior.phi_osc <- function(phi, a, b, l=(-0.998), u=1, lg=TRUE){
##     if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
##     if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
## }






## 
## function below requires all three meshes. These are obtained from
## the training data set and are passed in the function for
## computating integration weights.
## 
weights_line_segments_in_train <- function(X.test, Y.test, mesh, mesh.hd, mesh1d){
    X <- X.test
    Y <- Y.test
    nodes       <- c(mesh.hd$loc, 2*pi)
    ## intervals   <- head(cbind(nodes, lead(nodes)), -1)
    ## df.indices labels the indices of the knots for the directional, the spatial and the spatio-directional basis functions
    ## this data frame is created to create correspondences
    ## between spatio-directional basis knots with spatial basis knots, an
    ## between spatio-directional basis knots with head directional basis knots, respectively.
    df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)), space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
    ## So for example, if the spatio-directional basis knots are labeled as 1, 2, ..., p_Omega * p_Theta
    ## then the function mapindex2space.direction_basis takes as argument the label of spatio-directional basis knot
    ## and returns the coordinates and the head direction associated with the spatial basis function and the
    ## head directional basis function. This function uses mapindex2space.direction_basis which works similarly but
    ## instead of returning coords and angles, it returns the indices of the basis functions.
    mapindex2space.direction_index <- function(index){    
        f<-function(index.single){
            as.numeric(df.indices[which(df.indices$cross==index.single),c("dir","space")])
        }
        t((Vectorize(f, vectorize.args="index.single"))(index))
    }
    mapindex2space.direction_basis <- function(index){    
        f<-function(index.single){
            o <- as.numeric(df.indices[which(df.indices$cross==index.single),c("dir","space")])
            return(c(mesh.hd$loc[o[1]], mesh$loc[o[2],-3]))
        }
        t((Vectorize(f, vectorize.args="index.single"))(index))
    }
    ## 
    Ypos.tmp <- data.frame(
        index.CV = X$index.CV,
        hd=X$hd,
        time=X$synced_time,
        coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
        mutate(coords.lead = lead(coords)) %>%
        mutate(time.lead = lead(X$synced_time)) %>%
        mutate(hd.lead = lead(X$hd)) %>%
        head(-1)
    ## 
    Ypos.tmp <- Ypos.tmp %>% mutate(HD.split = map2(hd, hd.lead, split.arcs, mesh.hd=mesh.hd),
                                    L.arcs = lapply(HD.split,
                                                    function(x) apply(x, 1,
                                                                      function(y) abs(y[2]-y[1]))),
                                    time.split = pmap(list(time, time.lead, L.arcs), function(x,y,z){
                                        o <- interpolate(x,y,z)
                                        oo <- c(attr(o, "data"), y)
                                        ooo <- head(cbind(oo, lead(oo)), -1)
                                        colnames(ooo) <- NULL
                                        return(ooo)
                                    }),
                                    coords.split=pmap(list(coords, coords.lead, L.arcs), function(x,y,z){
                                        interpolate2(x, y, z)
                                    }),
                                    new.time = lapply(time.split, function(x) x[,1, drop=FALSE]),
                                    new.time.lead= lapply(time.split, function(x) x[,2, drop=FALSE]),
                                    new.hd = lapply(HD.split, function(x) x[,1, drop=FALSE]),
                                    new.hd.lead = lapply(HD.split, function(x) x[,2, drop=FALSE]),
                                    new.coords = lapply(coords.split, function(x) x[,1:2, drop=FALSE]),
                                    new.coords.lead = lapply(coords.split, function(x) x[,3:4, drop=FALSE])
                                    ) 
    ## 
    Ypos.tmp <- Ypos.tmp %>% dplyr::select(index.CV, new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
        unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))
    names(Ypos.tmp) <- c("index.CV", "time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")
    ## 
    line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
                                 filter.zero.length=FALSE,
                                 ep=Ypos.tmp$coords.lead, tol=.0)
    ## 
    df <- data.frame(origin=line.segments$split.origin,
                     filter.index=line.segments$filter.index,
                     sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                     ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
        group_by(origin) %>%
        summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
        mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
        mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 
    ## attribute named _data_ stores length of line segments, time differences and arclengths
    Ypos <- inner_join(Ypos.tmp %>%
                       mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
        mutate(Li = map2(sp, ep,
                         function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
        mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
        mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 
    filter.index  <- do.call("c", Ypos$filter.index)
    ## ------------------
    ## Integration points
    ## ------------------
    coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
    HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
    T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))
    index.CV     <- c(do.call("c", lapply(as.list(1:nrow(Ypos)), function(k) rep(Ypos$index.CV[k], nrow(Ypos$sp[[k]])) )), tail(Ypos$index.CV, 1))
    ## 
    Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
    Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
    Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
    A      <- inla.row.kron(Ahd, Aosc)
    ## 
    dGamma <- c(do.call("c", Ypos$Li))
    dT  <- diff(T.data)
    ## spatial basis functions
    Aosctmp <- as(Aosc, "dgTMatrix")
    Aosc.indices <- cbind(cbind(Aosctmp@i+1, Aosctmp@j+1), Aosctmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
    Aosc.indices <- Aosc.indices[order(Aosc.indices[,1]),] %>% as.data.frame %>% mutate(index.CV = sort(rep(index.CV,3))) # 
    ## spatial-directional basis functions
    Atmp <- as(A, "dgTMatrix")
    A.indices <- cbind(cbind(Atmp@i+1, Atmp@j+1), Atmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
    A.indices <- A.indices[order(A.indices[,1]),] %>% as.data.frame %>% mutate(index.CV = sort(rep(index.CV,6)))#
    ## temporal basis functions
    Attmp <- as(Atilde, "dgTMatrix")
    At.indices <- cbind(cbind(Attmp@i+1, Attmp@j+1), Attmp@x) # (i,j, A[i,j]) for which Atilde[i,j] is non-zero (Time)
    At.indices <- At.indices[order(At.indices[,1]),] %>% as.data.frame %>% mutate(index.CV = sort(rep(index.CV,2)))
    names(Aosc.indices) <- c("tk", "i", "psi.o", "index.CV") #ot: omega 
    names(A.indices) <- c("tk", "i", "psi.ot", "index.CV") #ot: omega x theta
    names(At.indices) <- c("tk", "l", "psi.t", "index.CV")
    ## 
    Aosc.indices.group.segments <- Aosc.indices
    while(TRUE){
        if(length(which(diff(unique(Aosc.indices.group.segments$index.CV)) == 1)) == 0){
            break
        }
        CV.groups <- unique(Aosc.indices.group.segments$index.CV)
        wh.group.to.pool.with.previous.group <- which(diff(unique(Aosc.indices.group.segments$index.CV)) == 1) + 1
        ## n.CV.groups <- length(unique(Aosc.indices.group.segments$index.CV))
        ## n.lines.in.groups <- table(Aosc.indices.group.segments$index.CV) %>% as.numeric
        Aosc.indices.group.segments <- Aosc.indices.group.segments %>%
            mutate(index.CV = case_when(
                       index.CV %in% CV.groups[wh.group.to.pool.with.previous.group] ~ (index.CV-1),
                       TRUE ~ index.CV
                   ))
    }
    ## 
    A.indices.group.segments <- A.indices
    while(TRUE){
        if(length(which(diff(unique(A.indices.group.segments$index.CV)) == 1)) == 0){
            break
        }
        CV.groups <- unique(A.indices.group.segments$index.CV)
        wh.group.to.pool.with.previous.group <- which(diff(unique(A.indices.group.segments$index.CV)) == 1) + 1
        ## n.CV.groups <- length(unique(A.indices.group.segments$index.CV))
        ## n.lines.in.groups <- table(A.indices.group.segments$index.CV) %>% as.numeric
        A.indices.group.segments <- A.indices.group.segments %>%
            mutate(index.CV = case_when(
                       index.CV %in% CV.groups[wh.group.to.pool.with.previous.group] ~ (index.CV-1),
                       TRUE ~ index.CV
                   ))
    }
    ## 
    At.indices.group.segments <- At.indices
    while(TRUE){
        if(length(which(diff(unique(At.indices.group.segments$index.CV)) == 1)) == 0){
            break
        }
        CV.groups <- unique(At.indices.group.segments$index.CV)
        wh.group.to.pool.with.previous.group <- which(diff(unique(At.indices.group.segments$index.CV)) == 1) + 1
        ## n.CV.groups <- length(unique(At.indices.group.segments$index.CV))
        ## n.lines.in.groups <- table(At.indices.group.segments$index.CV) %>% as.numeric
        At.indices.group.segments <- At.indices.group.segments %>%
            mutate(index.CV = case_when(
                       index.CV %in% CV.groups[wh.group.to.pool.with.previous.group] ~ (index.CV-1),
                       TRUE ~ index.CV
                   ))
    }
    df.prism.M0 <- Aosc.indices.group.segments %>% group_by(tk) %>%  nest %>% 
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma=c(dGamma,0),
            dGamma.lead = lead(dGamma),
            dGamma.lag = lag(dGamma),
            index.CV   = map(data, function(x){
                x$index.CV
            }),
            val.M0 = pmap(list(data, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.o[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M0=ooo)
                oooo
            })) %>% 
        dplyr::select(-c("data"))
    df.prism.M1_M2 <- full_join(At.indices.group.segments %>% group_by(tk) %>% nest(),
                                A.indices.group.segments %>% group_by(tk) %>% nest(), by=c("tk"="tk")) %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma=c(dGamma,0),
            dGamma.lead   = lead(dGamma),
            dGamma.lag    = lag(dGamma),
            index.M1.CV    = map(data.y, function(x){
                data.frame(index.CV=rep(unique(x$index.CV),6))
            }),
            index.M2.CV    = map(data.y, function(x){
                data.frame(index.CV=rep(unique(x$index.CV),12))
            }),
            val.M1 = pmap(list(data.y, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist(lapply(1:nrow(y), function(k) {
                    y$psi.ot[oo[k]]}))
                oooo <- data.frame(i=y$i[oo],val.M1=ooo)
                oooo
            }),
            val.M2 = pmap(list(data.x, data.y, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
                    x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            })) %>%
        dplyr::select(-c("data.x", "data.y"))
    ## !!
    df.prism.M1_M2.2 <- full_join(full_join(At.indices.group.segments %>% group_by(tk) %>% nest(),
                                            A.indices.group.segments %>% group_by(tk) %>% nest(), by="tk"),
                                  Aosc.indices.group.segments %>% group_by(tk) %>% nest(), by="tk") %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
            direction     = HD.data,
            direction.lag = c(0, HD.data[-length(direction)]),
            coords        = I(coords.trap),
            coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
            dGamma=c(dGamma,0),
            dGamma.lead   = lead(dGamma),
            dGamma.lag    = lag(dGamma),
            index.M1.CV    = map(data.y, function(x){
                data.frame(index.CV=rep(unique(x$index.CV), 6))
            }),
            index.M2.space.time.CV    = map(data, function(x){
                data.frame(index.CV=rep(unique(x$index.CV), 6))
            }),
            index.M2.CV    = map(data.y, function(x){
                data.frame(index.CV=rep(unique(x$index.CV),12))
            }),
            val.M1 = pmap(list(data.y, dGamma), function(y, z) {
                oo  <- 1:nrow(y)
                ooo <- unlist( lapply(1:nrow(y), function(k) y$psi.ot[oo[k]] ))
                oooo <- data.frame(i=y$i[oo],val.M1=ooo)
                oooo
            }),
            val.M2.space.time = pmap(list(data.x, data, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist( lapply(1:(nrow(x) * nrow(y)), function(k) x$psi.t[oo[k,1]] * y$psi.o[oo[k,2]] ))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            }),
            val.M2 = pmap(list(data.x, data.y, dGamma), function(x, y, z){
                oo  <- expand.grid(1:nrow(x), 1:nrow(y))
                ooo <- unlist( lapply(1:(nrow(x) * nrow(y)), function(k) x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]] ))
                oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
                oooo
            })) %>%
        dplyr::select(-c("data.x", "data.y", "data"))
    ## ## !!
    ## --------------
    df.prism.M0 <- df.prism.M0 %>% unnest(cols=c(val.M0, index.CV))
    ## --------------
    df.W.M0 <- df.prism.M0 %>% group_by(index.CV) %>% nest %>%
        mutate(
            df.W.M0 = map(data, function(x){
                tk.min = min(x$tk)
                rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
                      dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0),
                      x %>% 
                      filter(tk!=tk.min) %>%
                      mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                             group=tk-1,
                             dGamma=0) %>%
                      dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0)) %>%
                    arrange(group) %>%
                    mutate(dGamma.trap = dGamma + dGamma.lag) 
            })) %>% dplyr::select(-c("data")) %>%
        mutate(
            W.ipoints.M0 = map(df.W.M0, function(x){
                tol <- 0
                df.dGamma.sum.k.kplus1.M0 <- x %>% group_by(group, i) %>%
                    summarize(val = sum(pmax(dGamma.trap*val.M0, tol))/2,
                              time = unique(time),
                              direction=unique(direction),
                              coords=unique(coords))  %>%
                    ungroup %>% group_by(i) %>%
                    summarize(val = sum(val))
                ## 
                W.M0 <- sparseVector(i=df.dGamma.sum.k.kplus1.M0$i,
                                     x=df.dGamma.sum.k.kplus1.M0$val,
                                     length=mesh$n)
                z <- list()
                z$W.M0 <- W.M0
                ## 
                W.ipoints.M0 <- as(W.M0, "sparseMatrix")
                W.ipoints.M0 <- data.frame(coords.x1 = mesh$loc[W.ipoints.M0@i+1,1],
                                           coords.x2 = mesh$loc[W.ipoints.M0@i+1,2],
                                           weight=W.ipoints.M0@x)
                z$W.ipoints.M0 <- W.ipoints.M0
                return(z)
            })
        )
    ## ------------------------------
    df.prism.M1 <- df.prism.M1_M2 %>% dplyr::select(-val.M2) %>% 
        unnest(cols=c(val.M1, index.M1.CV))
    ## ------------------------------
        df.W.M1 <- df.prism.M1 %>% group_by(index.CV) %>% nest %>%
        mutate(
            df.W.M1 = map(data, function(x){
                tk.min = min(x$tk)           
                rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
                      dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1),
                      x %>% 
                      filter(tk!=tk.min) %>%
                      mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
                             group=tk-1,
                             dGamma=0) %>%
                      dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1)) %>%
                    arrange(group) %>%
                    mutate(dGamma.trap = dGamma + dGamma.lag)
            })) %>%dplyr::select(-c("data")) %>%
        mutate(
            W.ipoints.M1 = map(df.W.M1, function(x){
                tol <- 0    
                df.dGamma.sum.k.kplus1.M1 <- x %>% group_by(group, i) %>%
                    summarize(val = sum(pmax(dGamma.trap*val.M1, tol))/2,
                              time = unique(time),
                              direction=unique(direction),
                              coords=unique(coords))  %>%
                    ungroup %>% group_by(i) %>%
                    summarize(val = sum(val))
                W.M1 <- sparseVector(i=df.dGamma.sum.k.kplus1.M1$i,
                                     x=df.dGamma.sum.k.kplus1.M1$val,
                                     length=mesh$n * mesh.hd$n)
                z <- list()
                z$W.M1 <- W.M1
                ## 
                W.ipoints.M1 <- as(W.M1, "sparseMatrix")
                W.ipoints.M1 <- data.frame(hd=mapindex2space.direction_basis(W.ipoints.M1@i+1)[,1],
                                           coords.x1 = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,2],
                                           coords.x2 = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,3],
                                           weight=W.ipoints.M1@x)
                z$W.ipoints.M1 <- W.ipoints.M1
                return(z)
            })
        )
    ## !!-------------------------------------
    df.prism.M2.space.time <- df.prism.M1_M2.2 %>% dplyr::select(-c(val.M2,val.M1, index.M2.CV, index.M1.CV)) %>%
        unnest(cols=c(val.M2.space.time, index.M2.space.time.CV))
    ## !!-------------------------------------
        df.W.M2.space.time <- df.prism.M2.space.time %>% group_by(index.CV) %>% nest %>%
        mutate(
            df.W.M2.space.time = map(data, function(x){
                tk.min = min(x$tk)           
                rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
                      dplyr::select(group, time, coords, dGamma, dGamma.lag, l, i, val.M2),
                      x %>% 
                      filter(tk!=tk.min) %>%
                      mutate(time=time.lag, coords=coords.lag,
                             group=tk-1,
                             dGamma=0) %>%
                      dplyr::select(group, time, coords, dGamma, dGamma.lag, l, i, val.M2)) %>%
                    arrange(group) %>%
                    mutate(dGamma.trap = dGamma + dGamma.lag)
            })) %>%dplyr::select(-c("data")) %>%
        mutate(
            W.ipoints.M2.space.time = map(df.W.M2.space.time, function(x){
                tol <- 0    
                df.dGamma.sum.k.kplus1.M2.space.time <- x %>% group_by(group, l, i) %>%
                    summarize(val = sum(pmax(dGamma.trap*val.M2, tol))/2,
                              time = unique(time),
                              coords=unique(coords))
                W <- sparseMatrix(i=df.dGamma.sum.k.kplus1.M2.space.time$l,
                                  j=df.dGamma.sum.k.kplus1.M2.space.time$i,
                                  x=df.dGamma.sum.k.kplus1.M2.space.time$val/2)
                W <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(Aosc)-ncol(W)))
                W <- W %>% rbind(Matrix(0, nrow=ncol(Attmp)-nrow(W), ncol=ncol(W)))                
                W.ipoints.M2.space.time <- as(W, "dgTMatrix")
                W.ipoints.M2.space.time <- data.frame(firing_times=mesh1d$loc[W.ipoints.M2.space.time@i+1],
                                                      coords.x1 = mesh$loc[W.ipoints.M2.space.time@j+1,1],
                                                      ## mapindex2space.direction_basis(W.ipoints.M2.space.time@j+1)[,2]
                                                      coords.x2 =mesh$loc[W.ipoints.M2.space.time@j+1,2],
                                                      ## mapindex2space.direction_basis(W.ipoints.M2.space.time@j+1)[,3]
                                                      weight=W.ipoints.M2.space.time@x) %>% arrange(firing_times)                
                z <- list()
                z$W.space.time <- W
                z$W.ipoints.M2.space.time <- W.ipoints.M2.space.time
                return(z)
            })
        )
    ## !!-------------------------------------
    df.prism.M2 <- df.prism.M1_M2 %>% dplyr::select(-val.M1) %>%
        unnest(cols=c(val.M2, index.M2.CV))
    ## !!-------------------------------------
    ## 
    z <- list()
    z$df.W.M0            <- df.W.M0
    z$df.W.M1            <- df.W.M1
    z$df.W.M2.space.time <- df.W.M2.space.time
    return(z)
}


pred.mean.var <- function(weights.mat, post.sample, sample=FALSE){
                                        # posterior sample
    n <- length(post.sample) 
    Int <- lapply(1:n, function(i) post.sample[[i]]$Intercept) # simulated intercepts
    GRF <- lapply(1:n, function(i) post.sample[[i]]$spde2)     # simulated GRFs 
    cond.mean <- sapply(1:n, function(i) exp(Int[[i]])*sum(exp(GRF[[i]])*weights.mat)) # posterior means conditioned on latent parameters
    post.mean <- mean(cond.mean)
    post.var  <- mean(cond.mean) + var(cond.mean)       
    if(!sample) list(post.mean=post.mean, post.var=post.var)  else cond.mean
}

pred.mean.var.M2 <- function(weights.mat, post.sample, sample=FALSE){
                                        # posterior sample
    n <- length(post.sample) 
    Int <- lapply(1:n, function(i) post.sample[[i]]$Intercept) # simulated intercepts
    GRF1 <- lapply(1:n, function(i) post.sample[[i]]$spde1)    # simulated GRFs
    GRF2 <- lapply(1:n, function(i) post.sample[[i]]$spde2)    # simulated GRFs 
    cond.mean <- sapply(1:n, function(i) {
        as.numeric(exp(Int[[i]])*sum((matrix(exp(GRF1[[i]]),nrow=1) %*% weights.mat) * exp(GRF2[[i]])))
    })
    post.mean <- mean(cond.mean)
    post.var  <- mean(cond.mean) + var(cond.mean)
    if(!sample) list(post.mean=post.mean, post.var=post.var)  else cond.mean
    ## list(post.mean=post.mean, post.var=post.var)
}


predictive.M0 <- function(seq, weights.mat, post.sample){
                                        # posterior sample
    n <- length(post.sample) 
    Int <- lapply(1:n, function(i) post.sample[[i]]$Intercept) # simulated intercepts
    GRF <- lapply(1:n, function(i) post.sample[[i]]$spde2)     # simulated GRFs
    predictive <- do.call("rbind",lapply(1:n, function(i) dpois(seq, lambda=exp(Int[[i]])*sum(exp(GRF[[i]])*weights.mat))))
    list(predictive=predictive)
}

predictive.M2 <- function(seq, weights.mat, post.sample){
                                        # posterior sample
    n <- length(post.sample) 
    Int <- lapply(1:n, function(i) post.sample[[i]]$Intercept) # simulated intercepts
    GRF1 <- lapply(1:n, function(i) post.sample[[i]]$spde1)    # simulated GRFs
    GRF2 <- lapply(1:n, function(i) post.sample[[i]]$spde2)    # simulated GRFs 
    predictive <- do.call("rbind",
                          lapply(1:n,
                                 function(i) {
                                     dpois(seq, lambda=as.numeric(exp(Int[[i]])*sum((matrix(exp(GRF1[[i]]),nrow=1) %*% weights.mat) * exp(GRF2[[i]]))))
                                 }))
    list(predictive=predictive)
}

## precision matrix of oscillating field
##
## Input: 
##      - theta1: correlation length parameter. Note that
##                this parameter does not have the same interpretation
##                as in the usual Matern covariance family due to the oscillation
##                of the process. However, it can still be used to describe
##                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
##                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
##      - theta2: scale parameter that controls the variance of the field.
##                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
##                standard deviation of the process but is related to.
##      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
##                -1 < phi < 0: oscillating 
##                 0 < phi < 1: overdampened oscillating 
## Output:
##      - precision matrix for the weights at the mesh vertices 
##        
## 
osc.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    theta3 <- theta[3]
    rho     <- exp(theta1)
    kappa   <- sqrt(8)/rho
    sigma   <- exp(theta2)
    phi     <- (1-exp(-theta3))/(1+exp(-theta3))
    sincpth <- sqrt(1-phi^2)/acos(phi)
    tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
    ## ---------------------------------------
    o       <- inla.mesh.fem(mesh, order = 2)
    ## 
    tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*phi*(o$g1)) + o$g2)
    ## 
}

## precision matrix of oscillating field with cyclic B-splines
hd.bsp.precision <- function(theta, mesh){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    ## phi     <- (1-exp(-theta3))/(1+exp(-theta3))
    ## sincpth <- sqrt(1-phi^2)/acos(phi)
    sigma   <- exp(theta2)
    ## variance needs adjustment
    tausq   <- (pi+(cosh(pi*kappa)*sinh(pi*kappa))/kappa)/
        ((sigma^2)*(2*kappa)^2*sinh(pi*kappa)^2)
    fem       <- inla.mesh.fem(mesh=mesh, order = 2)
    tausq*((kappa^4)*(fem$c0) + (2*(kappa^2)*(fem$g1)) + fem$g2)
}

## theta <- seq(0, 2*pi, len=10)
## mesh.cyclic <- inla.mesh.1d(theta, boundary="cyclic", degree=1)
## fem <- inla.mesh.1d.fem(mesh.cyclic)


## precision matrix of temporal field
##
## Input: 
##      - theta1: correlation length parameter. Note that
##                this parameter does not have the same interpretation
##                as in the usual Matern covariance family due to the oscillation
##                of the process. However, it can still be used to describe
##                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
##                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
##      - theta2: scale parameter that controls the variance of the field.
##                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
##                standard deviation of the process but is related to.
##      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
##                -1 < phi < 0: oscillating 
##                 0 < phi < 1: overdampened oscillating 
## Output:
##      - precision matrix for the weights at the mesh vertices 
##        
## 

temp.precision <- function(theta, mesh, o=2){
    theta1 <- theta[1]
    theta2 <- theta[2]
    rho     <- exp(theta1)
    kappa   <- sqrt(8*(3/2))/rho
    sigma   <- exp(theta2)
    tausq   <- 1/(4*(sigma^2)*(kappa^3))
    ## ---------------------------------------
    o       <- inla.mesh.fem(mesh=mesh, order = o)
    ## 
    tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*(o$g1)) + o$g2)
    ## 
}


trajectory.length <- function(x,y){
    coords.start      <- cbind(x, y)[-length(x),]
    coords.lead       <- cbind(x, y)[-1,]
    coords.matrix     <- cbind(coords.start, coords.lead)
    coords.and.dist   <- cbind(coords.start, cumsum(apply((coords.lead - coords.start)^2,1, function(x) sqrt(sum(x)))))
    length.trajectory <- coords.and.dist[nrow(coords.and.dist),3]
    length.trajectory
}



## Simulation of homogeneous Poisson point process on trajectories
## arguments: constant rate function lambda
## trajectory: x and y coordinates
path.hpp <- function(lambda, x, y, time, hd){
    x.lead = lead(x)
    y.lead = lead(y)
    coords.start      <- cbind(x, y)[-length(x),]
    coords.lead       <- cbind(x.lead, y.lead)[-length(x),]
    coords.matrix     <- cbind(coords.start, coords.lead)
    coords.and.dist   <- cbind(coords.start, cumsum(apply((coords.lead - coords.start)^2,1, function(x) sqrt(sum(x)))))
    length.trajectory <- coords.and.dist[nrow(coords.and.dist),3]
    n                 <- rpois(1,lambda * length.trajectory)
    p.hpp             <- cumsum(rexp(n, lambda))
    sim.hpp <- lapply(as.list(p.hpp), function(z){
        where.on.path.index <- max(which( coords.and.dist[,3] <= z) )
        start.point         <- coords.start[where.on.path.index,]
        dist.diff           <- (z - coords.and.dist[where.on.path.index,3])
        unit.vector <- (coords.lead[where.on.path.index,]-coords.start[where.on.path.index,])/
            sqrt(sum((coords.lead[where.on.path.index,]-coords.start[where.on.path.index,])^2))
        z <- list(a   = start.point, x = start.point + dist.diff*unit.vector, b=coords.lead[where.on.path.index,],
                  hd.start = hd[where.on.path.index], hd.end = lead(hd)[where.on.path.index],
                  time.start = time[where.on.path.index], time.end = lead(time)[where.on.path.index])
        return(z)
    })
    ## return(sim.hpp)
    sp <- do.call("rbind",lapply(sim.hpp, function(y) y$a))
    ep <- do.call("rbind",lapply(sim.hpp, function(y) y$b))
    ip <- do.call("rbind",lapply(sim.hpp, function(y) y$x))
    hd.start <- do.call("rbind",lapply(sim.hpp, function(y) y$hd.start))
    time.end <- do.call("rbind",lapply(sim.hpp, function(y) y$time.end))
    time.start <- do.call("rbind",lapply(sim.hpp, function(y) y$time.start))
    hd.end <- do.call("rbind",lapply(sim.hpp, function(y) y$hd.end))
    volxc <- apply(ep-ip, 1, function(x) sqrt(sum(x^2)))
    volab <- apply(ep-sp, 1, function(x) sqrt(sum(x^2)))
    theta <- volxc/volab
    data.frame(firing_times = (1-theta)*time.start + (theta)*time.end, position_x=ip[,1], position_y=ip[,2], hd=(1-theta)*hd.start + (theta)*hd.end)
    ## z <- list()
    ## z$sp <- sp
    ## z$ep <- ep
    ## z$ip <- ip
    ## z$hd <- (1-theta)*hd.start + (theta)*hd.end
    ## z$theta <- theta
    ## data.frame(firing_times = p.hpp*(max(time)/length.trajectory), position_x=t(sim.hpp)[,1], position_y=t(sim.hpp)[,2])
    ## return(z)
    ## barycentric weights
}

## sim <- path.hpp(lambda=lambda, x=dat$X$position_x, y=dat$X$position_y, time=dat$X$synced_time, hd=dat$X$hd)
## plot(trajectory)
## points(sim$position_x, sim$position_y, col=2, pch=16, cex=.5)




## o$ip[1,]
## (o$theta[1])*o$sp[1,] + (1-o$theta[1])*o$ep[1,]
## (o$theta[2])*o$sp[2,] + (1-o$theta[2])*o$ep[2,]

## hd[where.on.path.index] + 



## n.hpp <- absGamma*max(lambda)       # n points from 
## p.hpp <- cumsum(rexp(n.hpp, max(lambdapred)))
## p.hpp <- p.hpp[p.hpp < max(Ypos$time)]
## ## next need to map p.hpp on locations (x,y) of given trajectory 

## Ysim.hpp <- t(sapply(p.hpp, function(d, d.vec=df.trajectory$time,
##                                      coords=df.trajectory$coords,
##                                      lead=df.trajectory$coords.lead){
##     ind         <- (max(which(d.vec <= d))) 
##     start.point <- coords[ind,]
##     dist.diff   <- (d - d.vec[ind])
##     unit.vector <- (lead[ind,]-coords[ind,])/sqrt(sum((lead[ind,]-coords[ind,])^2))
##     start.point + dist.diff*unit.vector
## }))








## split.lines <- function(mesh, sp, ep, filter.zero.length = TRUE, tol= 1e-8) {
##                                         # locations for splitting
##         loc = as.matrix(rbind(sp,ep))
##         idx = 1:dim(sp)[1]
##                                         # Filter out segments not on the mesh
##         t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
##         t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
##                                         # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
##         sp = matrix(sp[!((t1==0) | (t2==0)),], ncol=2)
##         ep = matrix(ep[!((t1==0) | (t2==0)),], ncol=2)
##         idx = idx[!((t1==0) | (t2==0))]
##         loc = as.matrix(rbind(sp,ep))
##                                         # Split them segments into parts
##         if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
##         np = dim(sp)[1]
##         sp.idx = t(rbind(1:np,np+1:np))
##         splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
##                                         #plot(data$mesh)
##                                         #points(loc)
##                                         #points(splt$split.loc,col="blue)
##         sp = matrix(splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]], ncol=2) # Start point of new segments
##         ep = matrix(splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]], ncol=2) # End points of new segments
##         idx = idx[splt$split.idx[,1]]
##         origin = splt$split.origin
##                                         # Filter out zero length segments
##         if ( filter.zero.length ) {
##             sl = apply((ep-sp)^2,MARGIN=1,sum)
##             sp = sp[!(sl<tol^2),]
##             ep = ep[!(sl<tol^2),]
##             origin = origin[!(sl<tol^2)]
##             idx = idx[!(sl<tol^2)]
##         }
##         return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
##     }






## ## ---------------------------------------
## ## evaluating multivariate joint density
## ## ---------------------------------------
## ldmvnorm <- function(x, m, Q){
##     ## ----------------------------------------------------------------
##     ## input: x (evaluation point), 
##     ##        m (mean vector), 
##     ##        Q (precision matrix/inverse covariance Sigma^{-1})
##     ## output: log density of multivariate normal distribution
##     ## ----------------------------------------------------------------
##     k      <- length(x)
##     xminmu <- Matrix(x-m,ncol=1)
##     L      <- (Matrix::Cholesky(Q, perm=TRUE) %>% (Matrix::expand))$L
##     halflogdetQ   <- sum(log(diag(L)))     #!!!
##     ## print(paste("halflogdetQ is:", halflogdetQ))
##     out    <- halflogdetQ - (k/2)*log(2*pi) - .5*sum((xminmu)* (Q %*% xminmu)) # remove matrix multiplication
##     out <- as.numeric(out)
##     return(out)
## }

## interpolate <- function(x, y, z){
##     interp <- (x) + (cumsum((z)/sum((z))))*((y) - (x))
##     Delta <- diff(c(x, interp))
##     attr(Delta, "data") <- c(x, interp[-length(interp)]) # note that
##                                                          # this way
##                                                          # the last
##                                                          # value from
##                                                          # the
##                                                          # combined
##                                                          # data vector
##                                                          # will be
##                                                          # missing so
##                                                          # needs to be
##                                                          # appended later
##     return(Delta)
## }

## interpolate2 <- function(x, y, z){
##     interp <- matrix(nrow=length(z), ncol=2)
##     for(i in 1:length(z)){
##         interp[i,] <- (x) + (sum(z[1:i])/sum((z)))*((y) - (x))
##     }
##     o <- rbind(x, interp)
##     oo <- cbind(head(o, -1), tail(o, -1))
##     rownames(oo) <- NULL
##     colnames(oo) <- NULL
##     return(oo)
## }



## split.lines <- function(mesh, sp, ep, filter.zero.length = TRUE, tol= 1e-8, return.filter.index=TRUE) {
##     ## locations for splitting
##     loc = as.matrix(rbind(sp,ep))
##     idx = 1:dim(sp)[1]
##     ## Filter out segments not on the mesh
##     t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
##     t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
##     ## if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
##     sp = matrix(sp[!((t1==0) | (t2==0)),], ncol=2)
##     ep = matrix(ep[!((t1==0) | (t2==0)),], ncol=2)
##     idx = idx[!((t1==0) | (t2==0))]
##     loc = as.matrix(rbind(sp,ep))
##     ## Split them segments into parts
##     if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
##     np = dim(sp)[1]
##     sp.idx = t(rbind(1:np,np+1:np))
##     splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
##     ##plot(data$mesh)
##     ##points(loc)
##     ##points(splt$split.loc,col="blue)
##     sp = matrix(splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]], ncol=2) # Start point of new segments
##     ep = matrix(splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]], ncol=2) # End points of new segments
##     idx = idx[splt$split.idx[,1]]
##     origin = splt$split.origin
##     sl = apply((ep-sp)^2,MARGIN=1,sum)
##     filter.index <- !(sl < tol^2)
##     ## Filter out zero length segments
##     if ( filter.zero.length ) {
##         ## sl = apply((ep-sp)^2,MARGIN=1,sum)
##         sp = sp[!(sl<tol^2),]
##         ep = ep[!(sl<tol^2),]
##         origin = origin[!(sl<tol^2)]
##         idx = idx[!(sl<tol^2)]
##         filter.index = filter.index[!(sl<tol^2)]
##         ## drop filter.index
##         return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
##     }
##     return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc, filter.index=filter.index))
## }


## split.arcs <- function(hd, hd.lead, mesh.hd){
##     nodes       <- c(mesh.hd$loc, 2*pi)
##     intervals   <- head(cbind(nodes, lead(nodes)), -1)
##     ## print(paste("hd: ",hd, "hd.lead: ",hd))
##     hd.int <- which(apply(intervals, 1, function(x) (x[1] <= hd) & (hd <= x[2])))
##     hd.lead.int <- which(apply(intervals, 1, function(x) (x[1] <= hd.lead) & (hd.lead <= x[2])))
##     if(hd.int==hd.lead.int){
##         return(matrix(c(hd, hd.lead),nrow=1))
##     }else{
##         if(hd < hd.lead){
##             if(hd.lead.int-hd.int == 1){
##                 rbind(c(hd, intervals[hd.int,2]),
##                       c(intervals[hd.lead.int,1], hd.lead))
##             }else{
##                 return(rbind(c(hd, intervals[(hd.int+1),1]),
##                              intervals[(hd.int+1):(hd.lead.int-1),],
##                              c(intervals[(hd.lead.int-1),2], hd.lead)))
##             }
##         }else{
##             if(hd.int-hd.lead.int == 1){
##                 rbind(c(hd, intervals[hd.int,1]),
##                       c(intervals[hd.lead.int,2], hd.lead))
##             }else{
##             return(rbind(c(hd, intervals[(hd.int-1),2]),
##                          intervals[(hd.int-1):(hd.lead.int+1),2:1],
##                          c(intervals[(hd.lead.int+1),1], hd.lead)))
##             }
##             ## rbind(c(hd, intervals[(hd.int-1),1]),
##             ##       apply((intervals[(hd.int-1):(hd.lead.int+1),]), 1, function(x)x),
##             ##       c(intervals[(hd.lead.int-1),2], hd.lead))
##         }
##     }
## }


## split.segments.wrapper.function <- function(X, mesh, mesh.hd){
##     Ypos.tmp <- data.frame(
##     hd=X$hd, time=X$synced_time,
##     coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
##         mutate(coords.lead = lead(coords)) %>%
##         mutate(time.lead = lead(X$synced_time)) %>%
##         mutate(hd.lead = lead(X$hd)) %>%
##         head(-1)
##     Ypos.tmp <- Ypos.tmp %>% mutate(HD.split = map2(hd, hd.lead, function(x, y) split.arcs(x,y, mesh.hd=mesh.hd)),
##                                     L.arcs = lapply(HD.split,
##                                                     function(x) apply(x, 1,
##                                                                       function(y) abs(y[2]-y[1]))),
##                                     time.split = pmap(list(time, time.lead, L.arcs), function(x,y,z){
##                                         o <- interpolate(x,y,z)
##                                         oo <- c(attr(o, "data"), y)
##                                         ooo <- head(cbind(oo, lead(oo)), -1)
##                                         colnames(ooo) <- NULL
##                                         return(ooo)
##                                     }),
##                                     coords.split=pmap(list(coords, coords.lead, L.arcs), function(x,y,z){
##                                         interpolate2(x, y, z)
##                                     }),
##                                     new.time = lapply(time.split, function(x) x[,1, drop=FALSE]),
##                                     new.time.lead= lapply(time.split, function(x) x[,2, drop=FALSE]),
##                                     new.hd = lapply(HD.split, function(x) x[,1, drop=FALSE]),
##                                     new.hd.lead = lapply(HD.split, function(x) x[,2, drop=FALSE]),
##                                     new.coords = lapply(coords.split, function(x) x[,1:2, drop=FALSE]),
##                                     new.coords.lead = lapply(coords.split, function(x) x[,3:4, drop=FALSE])
##                                     )
##     Ypos.tmp <- Ypos.tmp %>% dplyr::select(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
##         unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))
##     names(Ypos.tmp) <- c("time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")    
##     line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
##                                  filter.zero.length=FALSE,
##                                  ep=Ypos.tmp$coords.lead, tol=.0)
##     df <- data.frame(origin=line.segments$split.origin,
##                      filter.index=line.segments$filter.index,
##                      sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
##                      ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
##         group_by(origin) %>%
##         summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
##         mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
##         mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 
##     ## 
##     ## 
##     ## attribute named _data_ stores length of line segments, time differences and arclengths
##     Ypos <- inner_join(Ypos.tmp %>% mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
##         mutate(Li = map2(sp, ep, function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
##         mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
##         mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 
##     filter.index  <- do.call("c", Ypos$filter.index)
##     z <- list()
##     z$Ypos <- Ypos
##     z$filter.index  <- filter.index
##     z$line.segments <- line.segments
##     return(z)
## }


## ## df.prism.M0    <- df.prism.M0.wrapper(Aosc.indices = Aosc.indices, dGamma=dGamma, T.data=T.data, HD.data=HD.data,
## ##                                       coords.trap=coords.trap) %>% unnest(cols=c(val.M0))
## df.prism.M0.wrapper <- function(Aosc.indices, dGamma, T.data, HD.data, coords.trap) {
##     df.prism.M0 <- Aosc.indices %>% group_by(tk) %>% nest() %>%
##         arrange(tk) %>%
##         ungroup %>% 
##         mutate(
##             time          = T.data,
##             time.lag      = c(0, time[-length(time)]),
##             coords        = I(coords.trap),
##             coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
##             dGamma=c(dGamma,0),
##             dGamma.lead = lead(dGamma),
##             dGamma.lag = lag(dGamma),
##             val.M0 = pmap(list(data, dGamma), function(y, z) {
##                 oo  <- 1:nrow(y)
##                 ooo <- unlist(lapply(1:nrow(y), function(k) {
##                     y$psi.o[oo[k]]}))
##                 oooo <- data.frame(i=y$i[oo],val.M0=ooo)
##                 oooo
##             })) %>% 
##         dplyr::select(-c("data"))
##     return(df.prism.M0)
## }

## df.prism.M1.M2.wrapper <- function(At.indices, A.indices, T.data, dGamma, HD.data, coords.trap){
##     df.prism.M1_M2 <- full_join(At.indices %>% group_by(tk) %>% nest(),
##                                 A.indices %>% group_by(tk) %>% nest(), by="tk") %>%
##         arrange(tk) %>%
##         ungroup %>% 
##         mutate(
##             time          = T.data,
##             time.lag      = c(0, time[-length(time)]),
##             direction     = HD.data,
##             direction.lag = c(0, HD.data[-length(direction)]),
##             coords        = I(coords.trap),
##             coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
##             dGamma=c(dGamma,0),
##             dGamma.lead = lead(dGamma),
##             dGamma.lag = lag(dGamma),
##             val.M1 = pmap(list(data.y, dGamma), function(y, z) {
##                 oo  <- 1:nrow(y)
##                 ooo <- unlist(lapply(1:nrow(y), function(k) {
##                     y$psi.ot[oo[k]]}))
##                 oooo <- data.frame(i=y$i[oo],val.M1=ooo)
##                 oooo
##             }),
##             val.M2 = pmap(list(data.x, data.y, dGamma), function(x, y, z){
##                 oo  <- expand.grid(1:nrow(x), 1:nrow(y))
##                 ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
##                     x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
##                 oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
##                 oooo
##             })) %>%
##         dplyr::select(-c("data.x", "data.y")) 
## }

## ## --------------
## ## link functions
## ## --------------



## theta.2.phi   <- function(theta, l=NULL, u=NULL) {
##     if(is.null(l)) l <- get("l", envir = .GlobalEnv)
##     if(is.null(u)) u <- get("u", envir = .GlobalEnv)
##     res <- l + (u-l)*pcauchy(theta, location=0, scale=1)## (1/(1+exp(-theta)))
##     attr(res, "jacobian") <- (u-l)*dcauchy(theta, location=0, scale=1)
##     return(res)
## }

## ## theta.2.phi <- function(theta){
## ##     ## (1-exp(-theta))/(1+exp(-theta))
## ##     (1/(1+exp(-theta)))-1
## ## }
## ## 
## theta.2.sigma <- function(theta){
##     exp(theta)
## }
## ##
## theta.2.rho <- function(theta){
##     exp(theta) + 5
## }
## theta.2.rho.direction <- function(theta){
##     exp(theta) 
## }
## ## 
## theta.2.kappa.1d <- function(theta){
##     sqrt(8*(3/2))/exp(theta)       
## }

## ## helper functions for computing the
## ## posterior distribution of gridness score
## ## input is a fitted model with inlabru (e.g. fit.space and fit.space.direction)
## ## object returns attributes for summaries of the posterior
## posterior.spatial.gridness.score <- function(inlabru.fitted.object, theta.mapping){    
##     marg                  <- inla.tmarginal(theta.mapping, inlabru.fitted.object$marginals.hyperpar[["Theta3 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }
## posterior.spatial.standard.deviation <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta2 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }
## posterior.spatial.range <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta1 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }



## posterior.directional.standard.deviation <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta5 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }
## posterior.directional.range <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.rho.direction, inlabru.fitted.object$marginals.hyperpar[["Theta4 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }
## posterior.directional.kappa <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.kappa.1d, inlabru.fitted.object$marginals.hyperpar[["Theta4 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }


## posterior.temporal.standard.deviation <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta6 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }
## posterior.temporal.range <- function(inlabru.fitted.object){
##     marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta7 for f"]])
##     summaries             <- inla.zmarginal(marg, silent=TRUE)
##     hpd.interval          <- inla.hpdmarginal(0.95, marg)
##     attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
##     return(marg)
## }


## ## priors of hyperparameters
## ## prior.phi_osc <- function(phi, a, b, l, u, lg=TRUE){
## ##     if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
## ##     if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
## ## }

## prior.phi_osc <- function(phi, a, b, l=(-0.998), u=1, lg=TRUE){
##     if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
##     if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
## }

## circulant <-function(x) {
##     n <- length(x)
##     suppressWarnings( 
##         Matrix(matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n), sparse=TRUE)
##     )
## }


## ## 
## ## function below requires all three meshes. These are obtained from
## ## the training data set and are passed in the function computation
## ## of the integration weights.
## ## 
## weights_line_segments_in_train <- function(X.test, Y.test, mesh, mesh.hd, mesh1d){
##     X <- X.test
##     Y <- Y.test
##     nodes       <- c(mesh.hd$loc, 2*pi)
##     intervals   <- head(cbind(nodes, lead(nodes)), -1)
##     df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)), space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
##     mapindex2space.direction_index <- function(index){    
##         f<-function(index.single){
##             as.numeric(df.indices[which(df.indices$cross==index.single),c("dir","space")])
##         }
##         t((Vectorize(f, vectorize.args="index.single"))(index))
##     }
##     mapindex2space.direction_basis <- function(index){    
##         f<-function(index.single){
##             o <- as.numeric(df.indices[which(df.indices$cross==index.single),c("dir","space")])
##             return(c(mesh.hd$loc[o[1]], mesh$loc[o[2],-3]))
##         }
##         t((Vectorize(f, vectorize.args="index.single"))(index))
##     }
##     ## 
##     Ypos.tmp <- data.frame(
##         index.CV = X$index.CV,
##         hd=X$hd, time=X$synced_time,
##         coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
##         mutate(coords.lead = lead(coords)) %>%
##         mutate(time.lead = lead(X$synced_time)) %>%
##         mutate(hd.lead = lead(X$hd)) %>%
##         head(-1)
##     ## 
##     Ypos.tmp <- Ypos.tmp %>% mutate(HD.split = map2(hd, hd.lead, split.arcs, mesh.hd=mesh.hd),
##                                     L.arcs = lapply(HD.split,
##                                                     function(x) apply(x, 1,
##                                                                       function(y) abs(y[2]-y[1]))),
##                                     time.split = pmap(list(time, time.lead, L.arcs), function(x,y,z){
##                                         o <- interpolate(x,y,z)
##                                         oo <- c(attr(o, "data"), y)
##                                         ooo <- head(cbind(oo, lead(oo)), -1)
##                                         colnames(ooo) <- NULL
##                                         return(ooo)
##                                     }),
##                                     coords.split=pmap(list(coords, coords.lead, L.arcs), function(x,y,z){
##                                         interpolate2(x, y, z)
##                                     }),
##                                     new.time = lapply(time.split, function(x) x[,1, drop=FALSE]),
##                                     new.time.lead= lapply(time.split, function(x) x[,2, drop=FALSE]),
##                                     new.hd = lapply(HD.split, function(x) x[,1, drop=FALSE]),
##                                     new.hd.lead = lapply(HD.split, function(x) x[,2, drop=FALSE]),
##                                     new.coords = lapply(coords.split, function(x) x[,1:2, drop=FALSE]),
##                                     new.coords.lead = lapply(coords.split, function(x) x[,3:4, drop=FALSE])
##                                     ) 
##     ## 
##     Ypos.tmp <- Ypos.tmp %>% dplyr::select(index.CV, new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
##         unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))
##     names(Ypos.tmp) <- c("index.CV", "time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")
##     ## 
##     line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
##                                  filter.zero.length=FALSE,
##                                  ep=Ypos.tmp$coords.lead, tol=.0)
##     ## 
##     df <- data.frame(origin=line.segments$split.origin,
##                      filter.index=line.segments$filter.index,
##                      sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
##                      ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
##         group_by(origin) %>%
##         summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
##         mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
##         mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 
##     ## attribute named _data_ stores length of line segments, time differences and arclengths
##     Ypos <- inner_join(Ypos.tmp %>%
##                        mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
##         mutate(Li = map2(sp, ep,
##                          function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
##         mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
##         mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 
##     filter.index  <- do.call("c", Ypos$filter.index)
##     ## ------------------
##     ## Integration points
##     ## ------------------
##     coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
##     HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
##     T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))
##     index.CV     <- c(do.call("c", lapply(as.list(1:nrow(Ypos)), function(k) rep(Ypos$index.CV[k], nrow(Ypos$sp[[k]])) )), tail(Ypos$index.CV, 1))
##     ## 
##     Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
##     Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
##     Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
##     A      <- inla.row.kron(Ahd, Aosc)
##     ## 
##     dGamma <- c(do.call("c", Ypos$Li))
##     dT  <- diff(T.data)
##     ## spatial basis functions
##     Aosctmp <- as(Aosc, "dgTMatrix")
##     Aosc.indices <- cbind(cbind(Aosctmp@i+1, Aosctmp@j+1), Aosctmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
##     Aosc.indices <- Aosc.indices[order(Aosc.indices[,1]),] %>% as.data.frame %>% mutate(index.CV = sort(rep(index.CV,3))) # 
##     ## spatial-directional basis functions
##     Atmp <- as(A, "dgTMatrix")
##     A.indices <- cbind(cbind(Atmp@i+1, Atmp@j+1), Atmp@x) # (i,j, A[i,j]) for which A[i,j] is non-zero (Omega x Theta)
##     A.indices <- A.indices[order(A.indices[,1]),] %>% as.data.frame %>% mutate(index.CV = sort(rep(index.CV,6)))#
##     ## temporal basis functions
##     Attmp <- as(Atilde, "dgTMatrix")
##     At.indices <- cbind(cbind(Attmp@i+1, Attmp@j+1), Attmp@x) # (i,j, A[i,j]) for which Atilde[i,j] is non-zero (Time)
##     At.indices <- At.indices[order(At.indices[,1]),] %>% as.data.frame %>% mutate(index.CV = sort(rep(index.CV,2)))
##     names(Aosc.indices) <- c("tk", "i", "psi.o", "index.CV") #ot: omega 
##     names(A.indices) <- c("tk", "i", "psi.ot", "index.CV") #ot: omega x theta
##     names(At.indices) <- c("tk", "l", "psi.t", "index.CV")                                        #
##     Aosc.indices.group.segments <- Aosc.indices
##     while(TRUE){
##         if(length(which(diff(unique(Aosc.indices.group.segments$index.CV)) == 1)) == 0){
##             break
##         }
##         CV.groups <- unique(Aosc.indices.group.segments$index.CV)
##         wh.group.to.pool.with.previous.group <- which(diff(unique(Aosc.indices.group.segments$index.CV)) == 1) + 1
##         ## n.CV.groups <- length(unique(Aosc.indices.group.segments$index.CV))
##         ## n.lines.in.groups <- table(Aosc.indices.group.segments$index.CV) %>% as.numeric
##         Aosc.indices.group.segments <- Aosc.indices.group.segments %>%
##             mutate(index.CV = case_when(
##                        index.CV %in% CV.groups[wh.group.to.pool.with.previous.group] ~ (index.CV-1),
##                        TRUE ~ index.CV
##                    ))
##     }
##     ## 
##     A.indices.group.segments <- A.indices
##     while(TRUE){
##         if(length(which(diff(unique(A.indices.group.segments$index.CV)) == 1)) == 0){
##             break
##         }
##         CV.groups <- unique(A.indices.group.segments$index.CV)
##         wh.group.to.pool.with.previous.group <- which(diff(unique(A.indices.group.segments$index.CV)) == 1) + 1
##         ## n.CV.groups <- length(unique(A.indices.group.segments$index.CV))
##         ## n.lines.in.groups <- table(A.indices.group.segments$index.CV) %>% as.numeric
##         A.indices.group.segments <- A.indices.group.segments %>%
##             mutate(index.CV = case_when(
##                        index.CV %in% CV.groups[wh.group.to.pool.with.previous.group] ~ (index.CV-1),
##                        TRUE ~ index.CV
##                    ))
##     }
##     ## 
##     At.indices.group.segments <- At.indices
##     while(TRUE){
##         if(length(which(diff(unique(At.indices.group.segments$index.CV)) == 1)) == 0){
##             break
##         }
##         CV.groups <- unique(At.indices.group.segments$index.CV)
##         wh.group.to.pool.with.previous.group <- which(diff(unique(At.indices.group.segments$index.CV)) == 1) + 1
##         ## n.CV.groups <- length(unique(At.indices.group.segments$index.CV))
##         ## n.lines.in.groups <- table(At.indices.group.segments$index.CV) %>% as.numeric
##         At.indices.group.segments <- At.indices.group.segments %>%
##             mutate(index.CV = case_when(
##                        index.CV %in% CV.groups[wh.group.to.pool.with.previous.group] ~ (index.CV-1),
##                        TRUE ~ index.CV
##                    ))
##     }
##     df.prism.M0 <- Aosc.indices.group.segments %>% group_by(tk) %>%  nest %>% 
##         arrange(tk) %>%
##         ungroup %>% 
##         mutate(
##             time          = T.data,
##             time.lag      = c(0, time[-length(time)]),
##             direction     = HD.data,
##             direction.lag = c(0, HD.data[-length(direction)]),
##             coords        = I(coords.trap),
##             coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
##             dGamma=c(dGamma,0),
##             dGamma.lead = lead(dGamma),
##             dGamma.lag = lag(dGamma),
##             index.CV   = map(data, function(x){
##                 x$index.CV
##             }),
##             val.M0 = pmap(list(data, dGamma), function(y, z) {
##                 oo  <- 1:nrow(y)
##                 ooo <- unlist(lapply(1:nrow(y), function(k) {
##                     y$psi.o[oo[k]]}))
##                 oooo <- data.frame(i=y$i[oo],val.M0=ooo)
##                 oooo
##             })) %>% 
##         dplyr::select(-c("data"))
##     df.prism.M1_M2 <- full_join(At.indices.group.segments %>% group_by(tk) %>% nest(),
##                                 A.indices.group.segments %>% group_by(tk) %>% nest(), by=c("tk"="tk")) %>%
##         arrange(tk) %>%
##         ungroup %>% 
##         mutate(
##             time          = T.data,
##             time.lag      = c(0, time[-length(time)]),
##             direction     = HD.data,
##             direction.lag = c(0, HD.data[-length(direction)]),
##             coords        = I(coords.trap),
##             coords.lag    = I(rbind(c(0, 0), coords.trap[-nrow(coords.trap),])),
##             dGamma=c(dGamma,0),
##             dGamma.lead   = lead(dGamma),
##             dGamma.lag    = lag(dGamma),
##             index.M1.CV    = map(data.y, function(x){
##                 data.frame(index.CV=rep(unique(x$index.CV),6))
##             }),
##             index.M2.CV    = map(data.y, function(x){
##                 data.frame(index.CV=rep(unique(x$index.CV),12))
##             }),
##             val.M1 = pmap(list(data.y, dGamma), function(y, z) {
##                 oo  <- 1:nrow(y)
##                 ooo <- unlist(lapply(1:nrow(y), function(k) {
##                     y$psi.ot[oo[k]]}))
##                 oooo <- data.frame(i=y$i[oo],val.M1=ooo)
##                 oooo
##             }),
##             val.M2 = pmap(list(data.x, data.y, dGamma), function(x, y, z){
##                 oo  <- expand.grid(1:nrow(x), 1:nrow(y))
##                 ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
##                     x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
##                 oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val.M2=ooo)
##                 oooo
##             })) %>%
##         dplyr::select(-c("data.x", "data.y"))
##     ## 
##     df.prism.M0 <- df.prism.M0 %>% unnest(cols=c(val.M0, index.CV))
##     df.prism.M1 <- df.prism.M1_M2 %>% dplyr::select(-val.M2) %>% 
##         unnest(cols=c(val.M1, index.M1.CV))
##     df.prism.M2 <- df.prism.M1_M2 %>% dplyr::select(-val.M1) %>%
##         unnest(cols=c(val.M2, index.M2.CV))
##     df.W.M0 <- df.prism.M0 %>% group_by(index.CV) %>% nest %>%
##         mutate(
##             df.W.M0 = map(data, function(x){
##                 tk.min = min(x$tk)
##                 rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
##                       dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0),
##                       x %>% 
##                       filter(tk!=tk.min) %>%
##                       mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
##                              group=tk-1,
##                              dGamma=0) %>%
##                       dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M0)) %>%
##                     arrange(group) %>%
##                     mutate(dGamma.trap = dGamma + dGamma.lag) 
##             })) %>% dplyr::select(-c("data")) %>%
##         mutate(
##             W.ipoints.M0 = map(df.W.M0, function(x){
##                 tol <- 0
##                 df.dGamma.sum.k.kplus1.M0 <- x %>% group_by(group, i) %>%
##                     summarize(val = sum(pmax(dGamma.trap*val.M0, tol))/2,
##                               time = unique(time),
##                               direction=unique(direction),
##                               coords=unique(coords))  %>%
##                     ungroup %>% group_by(i) %>%
##                     summarize(val = sum(val))
##                 ## 
##                 W.M0 <- sparseVector(i=df.dGamma.sum.k.kplus1.M0$i,
##                                      x=df.dGamma.sum.k.kplus1.M0$val,
##                                      length=mesh$n)
##                 z <- list()
##                 z$W.M0 <- W.M0
##                 ## 
##                 W.ipoints.M0 <- as(W.M0, "sparseMatrix")
##                 W.ipoints.M0 <- data.frame(coords.x1 = mesh$loc[W.ipoints.M0@i+1,1],
##                                            coords.x2 = mesh$loc[W.ipoints.M0@i+1,2],
##                                            weight=W.ipoints.M0@x)
##                 z$W.ipoints.M0 <- W.ipoints.M0
##                 return(z)
##             })
##         )
##     ## 
##     df.W.M1 <- df.prism.M1 %>% group_by(index.CV) %>% nest %>%
##         mutate(
##             df.W.M1 = map(data, function(x){
##                 tk.min = min(x$tk)           
##                 rbind(x %>% mutate(group=tk, dGamma.lag=0) %>%
##                       dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1),
##                       x %>% 
##                       filter(tk!=tk.min) %>%
##                       mutate(time=time.lag, direction=direction.lag, coords=coords.lag,
##                              group=tk-1,
##                              dGamma=0) %>%
##                       dplyr::select(group, time, direction, coords, dGamma, dGamma.lag, i, val.M1)) %>%
##                     arrange(group) %>%
##                     mutate(dGamma.trap = dGamma + dGamma.lag)
##             })) %>%dplyr::select(-c("data")) %>%
##         mutate(
##             W.ipoints.M1 = map(df.W.M1, function(x){
##                 tol <- 0    
##                 df.dGamma.sum.k.kplus1.M1 <- x %>% group_by(group, i) %>%
##                     summarize(val = sum(pmax(dGamma.trap*val.M1, tol))/2,
##                               time = unique(time),
##                               direction=unique(direction),
##                               coords=unique(coords))  %>%
##                     ungroup %>% group_by(i) %>%
##                     summarize(val = sum(val))
##                 W.M1 <- sparseVector(i=df.dGamma.sum.k.kplus1.M1$i,
##                                      x=df.dGamma.sum.k.kplus1.M1$val,
##                                      length=mesh$n * mesh.hd$n)
##                 z <- list()
##                 z$W.M1 <- W.M1
##                 ## 
##                 W.ipoints.M1 <- as(W.M1, "sparseMatrix")
##                 W.ipoints.M1 <- data.frame(hd=mapindex2space.direction_basis(W.ipoints.M1@i+1)[,1],
##                                            coords.x1 = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,2],
##                                            coords.x2 = mapindex2space.direction_basis(W.ipoints.M1@i+1)[,3],
##                                            weight=W.ipoints.M1@x)
##                 z$W.ipoints.M1 <- W.ipoints.M1
##                 return(z)
##             })
##         )
##     z <- list()
##     z$df.W.M0 <- df.W.M0
##     z$df.W.M1 <- df.W.M1
##     return(z)
## }


## ## precision matrix of oscillating field
## ##
## ## Input: 
## ##      - theta1: correlation length parameter. Note that
## ##                this parameter does not have the same interpretation
## ##                as in the usual Matern covariance family due to the oscillation
## ##                of the process. However, it can still be used to describe
## ##                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
## ##                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
## ##      - theta2: scale parameter that controls the variance of the field.
## ##                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
## ##                standard deviation of the process but is related to.
## ##      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
## ##                -1 < phi < 0: oscillating 
## ##                 0 < phi < 1: overdampened oscillating 
## ## Output:
## ##      - precision matrix for the weights at the mesh vertices 
## ##        
## ## 
## osc.precision <- function(theta, mesh){
##     theta1 <- theta[1]
##     theta2 <- theta[2]
##     theta3 <- theta[3]
##     rho     <- exp(theta1)
##     kappa   <- sqrt(8)/rho
##     sigma   <- exp(theta2)
##     phi     <- (1-exp(-theta3))/(1+exp(-theta3))
##     sincpth <- sqrt(1-phi^2)/acos(phi)
##     tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
##     ## ---------------------------------------
##     o       <- inla.mesh.fem(mesh, order = 2)
##     ## 
##     tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*phi*(o$g1)) + o$g2)
##     ## 
## }

## ## precision matrix of oscillating field with cyclic B-splines
## hd.bsp.precision <- function(theta, mesh){
##     theta1 <- theta[1]
##     theta2 <- theta[2]
##     rho     <- exp(theta1)
##     kappa   <- sqrt(8*(3/2))/rho
##     ## phi     <- (1-exp(-theta3))/(1+exp(-theta3))
##     ## sincpth <- sqrt(1-phi^2)/acos(phi)
##     sigma   <- exp(theta2)
##     ## variance needs adjustment
##     tausq   <- (pi+(cosh(pi*kappa)*sinh(pi*kappa))/kappa)/
##         ((sigma^2)*(2*kappa)^2*sinh(pi*kappa)^2)
##     fem       <- inla.mesh.fem(mesh=mesh, order = 2)
##     tausq*((kappa^4)*(fem$c0) + (2*(kappa^2)*(fem$g1)) + fem$g2)
## }

## ## precision matrix of temporal field
## ##
## ## Input: 
## ##      - theta1: correlation length parameter. Note that
## ##                this parameter does not have the same interpretation
## ##                as in the usual Matern covariance family due to the oscillation
## ##                of the process. However, it can still be used to describe
## ##                the length scale beyond which correlation decays (e.g., as with dampened sinusoids)
## ##                exp(theta1) = rho = sqrt(8 nu)/kappa, nu = alpha - d/2 = 2-(2/2) = 1
## ##      - theta2: scale parameter that controls the variance of the field.
## ##                exp(theta2) = sigma. Similar to theta1, sigma here does not equal the marginal
## ##                standard deviation of the process but is related to.
## ##      - theta3: oscillation parameter phi = 1-exp(-theta3)/(1+exp(-theta3))
## ##                -1 < phi < 0: oscillating 
## ##                 0 < phi < 1: overdampened oscillating 
## ## Output:
## ##      - precision matrix for the weights at the mesh vertices 
## ##        
## ## 

## temp.precision <- function(theta, mesh, o=2){
##     theta1 <- theta[1]
##     theta2 <- theta[2]
##     rho     <- exp(theta1)
##     kappa   <- sqrt(8*(3/2))/rho
##     sigma   <- exp(theta2)
##     tausq   <- 1/(4*(sigma^2)*(kappa^3))
##     ## ---------------------------------------
##     o       <- inla.mesh.fem(mesh=mesh, order = o)
##     ## 
##     tausq*((kappa^4)*(o$c0) + (2*(kappa^2)*(o$g1)) + o$g2)
##     ## 
## }



