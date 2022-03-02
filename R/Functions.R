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
    nodes       <- c(mesh.hd$loc, 2*pi)
    intervals   <- head(cbind(nodes, lead(nodes)), -1)
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


split.segments.wrapper.function <- function(X, mesh, mesh.hd){
    Ypos.tmp <- data.frame(
    hd=X$hd, time=X$synced_time,
    coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
        mutate(coords.lead = lead(coords)) %>%
        mutate(time.lead = lead(X$synced_time)) %>%
        mutate(hd.lead = lead(X$hd)) %>%
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
    Ypos.tmp <- Ypos.tmp %>% dplyr::select(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
        unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))
    names(Ypos.tmp) <- c("time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")    
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
    Ypos <- inner_join(Ypos.tmp %>% mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
        mutate(Li = map2(sp, ep, function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
        mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
        mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 
    filter.index  <- do.call("c", Ypos$filter.index)
    z <- list()
    z$Ypos <- Ypos
    z$filter.index  <- filter.index
    z$line.segments <- line.segments
    return(z)
}


## df.prism.M0    <- df.prism.M0.wrapper(Aosc.indices = Aosc.indices, dGamma=dGamma, T.data=T.data, HD.data=HD.data,
##                                       coords.trap=coords.trap) %>% unnest(cols=c(val.M0))
df.prism.M0.wrapper <- function(Aosc.indices, dGamma, T.data, HD.data, coords.trap) {
    df.prism.M0 <- Aosc.indices %>% group_by(tk) %>% nest() %>%
        arrange(tk) %>%
        ungroup %>% 
        mutate(
            time          = T.data,
            time.lag      = c(0, time[-length(time)]),
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

## --------------
## link functions
## --------------



theta.2.phi   <- function(theta, l=NULL, u=NULL) {
    if(is.null(l)) l <- get("l", envir = .GlobalEnv)
    if(is.null(u)) u <- get("u", envir = .GlobalEnv)
    res <- l + (u-l)*pcauchy(theta, location=0, scale=1)## (1/(1+exp(-theta)))
    attr(res, "jacobian") <- (u-l)*dcauchy(theta, location=0, scale=1)
    return(res)
}

## theta.2.phi <- function(theta){
##     ## (1-exp(-theta))/(1+exp(-theta))
##     (1/(1+exp(-theta)))-1
## }
## 
theta.2.sigma <- function(theta){
    exp(theta)
}
##
theta.2.rho <- function(theta){
    exp(theta) + 5
}
theta.2.rho.direction <- function(theta){
    exp(theta) 
}
## 
theta.2.kappa.1d <- function(theta){
    sqrt(8*(3/2))/exp(theta)       
}

## helper functions for computing the
## posterior distribution of gridness score
## input is a fitted model with inlabru (e.g. fit.space and fit.space.direction)
## object returns attributes for summaries of the posterior
posterior.spatial.gridness.score <- function(inlabru.fitted.object, theta.mapping){    
    marg                  <- inla.tmarginal(theta.mapping, inlabru.fitted.object$marginals.hyperpar[["Theta3 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.spatial.standard.deviation <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta2 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.spatial.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta1 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}



posterior.directional.standard.deviation <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta5 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.directional.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho.direction, inlabru.fitted.object$marginals.hyperpar[["Theta4 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.directional.kappa <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.kappa.1d, inlabru.fitted.object$marginals.hyperpar[["Theta4 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}


posterior.temporal.standard.deviation <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta6 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.temporal.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta7 for f"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}


## priors of hyperparameters
## prior.phi_osc <- function(phi, a, b, l, u, lg=TRUE){
##     if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
##     if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
## }

prior.phi_osc <- function(phi, a, b, l=(-0.998), u=1, lg=TRUE){
    if(lg)  return(-log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE))
    if(!lg)  return((1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b))
}

circulant <-function(x) {
    n <- length(x)
    suppressWarnings( 
        Matrix(matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n), sparse=TRUE)
    )
}


## 
## function below requires all three meshes. These are obtained from
## the training data set and are passed in the function computation
## of the integration weights.
## 
weights_line_segments_in_train <- function(X.test, Y.test, mesh, mesh.hd, mesh1d){
    X <- X.test
    Y <- Y.test
    nodes       <- c(mesh.hd$loc, 2*pi)
    intervals   <- head(cbind(nodes, lead(nodes)), -1)
    df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)), space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
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
        hd=X$hd, time=X$synced_time,
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
    names(At.indices) <- c("tk", "l", "psi.t", "index.CV")                                        #
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
    ## 
    df.prism.M0 <- df.prism.M0 %>% unnest(cols=c(val.M0, index.CV))
    df.prism.M1 <- df.prism.M1_M2 %>% dplyr::select(-val.M2) %>% 
        unnest(cols=c(val.M1, index.M1.CV))
    df.prism.M2 <- df.prism.M1_M2 %>% dplyr::select(-val.M1) %>%
        unnest(cols=c(val.M2, index.M2.CV))
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
    ## 
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
    z <- list()
    z$df.W.M0 <- df.W.M0
    z$df.W.M1 <- df.W.M1
    return(z)
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



## 
## ipoints2 <- function(samplers = NULL, domain = NULL, name = NULL, group = NULL,  int.args = NULL, project = NULL){
##     int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
##     if (is.null(int.args)) {
##         int.args <- list()
##     }
##     missing.args <- setdiff(names(int.args.default), names(int.args))
##     int.args[missing.args] <- int.args.default[missing.args]
##     if (!is.null(int.args[["nsub"]])) {
##         int.args[["nsub1"]] <- int.args[["nsub"]]
##     }
##     if (!is.null(int.args[["nsub"]])) {
##         int.args[["nsub2"]] <- int.args[["nsub"]]
##     }
##     if (!is.null(project)) {
##         if (project && !identical(int.args$method, "stable")) {
##             stop("ipoints(project=TRUE) is deprecated, and int.args$methods != 'stable'")
##         }
##         else if (!project && identical(int.args$method, "stable")) {
##             stop("ipoints(project=FALSE) is deprecated, and int.args$methods == 'stable'")
##         }
##         warning("ipoints(project=", ifelse(project, "TRUE", "FALSE"), 
##             ") is deprecated. Will use int.args$method = '", 
##             int.args[["method"]], "' instead.")
##     }
##     if (is.null(domain) && inherits(samplers, c("inla.mesh.1d", 
##         "inla.mesh"))) {
##         domain <- samplers
##         samplers <- NULL
##     }
##     is_2d <- (!is.null(samplers) && inherits(samplers, c("SpatialPoints", 
##         "SpatialPointsDataFrame", "SpatialPolygons", "SpatialPolygonsDataFrame", 
##         "SpatialLines", "SpatialLinesDataFrame"))) || inherits(domain, 
##         "inla.mesh")
##     is_1d <- !is_2d && ((!is.null(samplers) && is.numeric(samplers)) || 
##         (!is.null(domain) && (is.numeric(domain) || inherits(domain, 
##             "inla.mesh.1d"))))
##     if (!is_1d && !is_2d) {
##         stop("Unable to determine integration domain definition")
##     }
##     if (is_1d && !is.null(samplers) && !is.null(domain) && is.numeric(domain) && 
##         length(domain) == 1) {
##         int.args[["nsub1"]] <- domain
##         domain <- NULL
##         int.args[["method"]] <- "direct"
##     }
##     if (is_2d && !is.null(samplers) && !is.null(domain) && is.numeric(domain) && 
##         length(domain) == 1) {
##         int.args[["nsub2"]] <- domain
##         domain <- NULL
##         int.args[["method"]] <- "direct"
##     }
##     if (is.null(domain) && inherits(samplers, c("inla.mesh.1d", 
##         "inla.mesh"))) {
##         domain <- samplers
##         samplers <- NULL
##     }
##     if (is_1d && is.null(name)) {
##         name <- "x"
##     }
##     if (is.data.frame(samplers)) {
##         if (!("weight" %in% names(samplers))) {
##             samplers$weight <- 1
##         }
##         ips <- samplers
##     }
##     else if (is_1d && is.null(samplers) && is.numeric(domain)) {
##         ips <- data.frame(x = as.vector(domain), weight = 1)
##         colnames(ips) <- c(name, "weight")
##     }
##     else if (is_1d && is.null(domain) && is.integer(samplers)) {
##         ips <- data.frame(x = as.vector(samplers), weight = 1)
##         colnames(ips) <- c(name, "weight")
##     }
##     else if (is_1d && is.null(samplers) && inherits(domain, "inla.mesh.1d") && 
##         identical(int.args[["method"]], "stable")) {
##         ips <- data.frame(x = domain$loc, weight = Matrix::diag(INLA::inla.mesh.fem(domain)$c0))
##         colnames(ips) <- c(name, "weight")
##     }
##     else if (is_1d) {
##         domain_range <- if (inherits(domain, "inla.mesh.1d")) {
##             domain$interval
##         }
##         else {
##             NULL
##         }
##         if (is.null(samplers)) {
##             samplers <- matrix(domain_range, 1, 2)
##         }
##         else {
##             if (is.null(dim(samplers))) {
##                 samplers <- matrix(samplers, nrow = 1)
##             }
##             if (ncol(samplers) != 2) {
##                 stop("Interval description matrix must have 2 elements or be a 2-column matrix.")
##             }
##             if (is.null(domain)) {
##                 domain <- INLA::inla.mesh.1d(sort(unique(as.vector(samplers))))
##             }
##         }
##         ips <- list()
##         if (domain$degree >= 2) {
##             warning("Integration points projected onto knots may lead to instability for degree >= 2 basis functions.")
##         }
##         nsub <- int.args[["nsub1"]]
##         u <- rep((seq_len(nsub) - 0.5)/nsub, domain$n - 1)
##         int_loc <- domain$loc[rep(seq_len(domain$n - 1), each = nsub)] * 
##             (1 - u) + domain$loc[rep(seq_len(domain$n - 1) + 
##             1, each = nsub)] * u
##         int_w <- (domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] - 
##             domain$loc[rep(seq_len(domain$n - 1), each = nsub)])/nsub
##         for (j in seq_len(nrow(samplers))) {
##             subsamplers <- samplers[j, ]
##             if (identical(int.args[["method"]], "stable")) {
##                 A_w <- INLA::inla.spde.make.A(domain, int_loc, 
##                   weights = int_w * (int_loc >= min(subsamplers)) * 
##                     (int_loc <= max(subsamplers)))
##                 ips[[j]] <- data.frame(loc = domain$loc, weight = Matrix::colSums(A_w))
##             }
##             else {
##                 inside <- (int_loc >= min(subsamplers)) & (int_loc <= 
##                   max(subsamplers))
##                 ips[[j]] <- data.frame(loc = int_loc[inside], 
##                   weight = int_w[inside])
##             }
##             colnames(ips[[j]]) <- c(name, "weight")
##         }
##         ips <- do.call(rbind, ips)
##     }
##     else if (inherits(domain, "inla.mesh") && is.null(samplers) && 
##         identical(int.args[["method"]], "stable")) {
##         coord_names <- c("x", "y", "coordinateZ")
##         if (!is.null(name)) {
##             coord_names[seq_along(name)] <- name
##         }
##         if (!fm_crs_is_null(domain$crs)) {
##             crs <- domain$crs
##             samplers <- stransform(domain, crs = CRS("+proj=cea +units=km"))
##         }
##         ips <- vertices(domain)
##         ips$weight <- INLA::inla.mesh.fem(domain, order = 1)$va
##         if (!fm_crs_is_null(domain$crs)) {
##             ips <- stransform(ips, crs = crs)
##         }
##         coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
##     }
##     else if (class(samplers) == "SpatialPoints") {
##         ips <- samplers
##         ips$weight <- 1
##     }
##     else if (class(samplers) == "SpatialPointsDataFrame") {
##         if (!("weight" %in% names(samplers))) {
##             warning("The integration points provided have no weight column. Setting weights to 1.")
##             samplers$weight <- 1
##         }
##         ips <- samplers
##     }
##     else if (inherits(samplers, "SpatialLines") || inherits(samplers, 
##         "SpatialLinesDataFrame")) {
##         if (inherits(samplers, "SpatialLines") && !inherits(samplers, 
##             "SpatialLinesDataFrame")) {
##             samplers <- SpatialLinesDataFrame(samplers, data = data.frame(weight = rep(1, 
##                 length(samplers))))
##         }
##         if (!("weight" %in% names(samplers))) {
##             samplers$weight <- 1
##         }
##         ips <- int.slines(samplers, domain, group = group, project = identical(int.args[["method"]], 
##             "stable"))
##         coord_names <- c("x", "y", "coordinateZ")
##         if (!is.null(coordnames(samplers))) {
##             coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
##         }
##         else if (!is.null(name)) {
##             coord_names[seq_along(name)] <- name
##         }
##         coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
##     }
##     else if (is_2d && (inherits(samplers, c("SpatialPolygons", 
##         "SpatialPolygonsDataFrame")) || is.null(samplers))) {
##         if (is.null(samplers)) {
##             stop("Direct integration scheme for mesh domain with no samplers is not yet implemented.")
##         }
##         if (class(samplers)[1] == "SpatialPolygons") {
##             samplers <- SpatialPolygonsDataFrame(samplers, data = data.frame(weight = rep(1, 
##                 length(samplers))), match.ID = FALSE)
##         }
##         else if (is.null(samplers@data[["weight"]])) {
##             samplers@data[["weight"]] <- 1
##         }
##         cnames <- coordnames(samplers)
##         samplers_crs <- fm_sp_get_crs(samplers)
##         if (!inlabru:::fm_crs_is_null(domain$crs)) {
##             samplers <- stransform(samplers, crs = sp::CRS("+proj=cea +units=km"))
##         }
##         polyloc <- do.call(rbind, lapply(seq_len(length(samplers)), 
##             function(k) {
##                 cbind(x = rev(coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 
##                   1]), y = rev(coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 
##                   2]), group = k)
##             }))
##         poly_segm <- INLA::inla.sp2segment(samplers, join = FALSE)
##         poly_segm <- lapply(seq_along(poly_segm), function(k) {
##             segm <- poly_segm[[k]]
##             segm[["grp"]] <- rep(k, NROW(segm[["idx"]]))
##             segm[["is.bnd"]] <- TRUE
##             segm
##         })
##         if (is.null(domain)) {
##             warning("Computing integration points from polygon; specify a mesh for better numerical control.")
##             max.edge <- max(diff(range(polyloc[, 1])), diff(range(polyloc[, 
##                 2])))/20
##             domain <- INLA::inla.mesh.2d(boundary = samplers, 
##                 max.edge = max.edge)
##             domain$crs <- fm_sp_get_crs(samplers)
##         }
##         else {
##             if (!inlabru:::fm_crs_is_null(domain$crs)) {
##                 domain <- stransform(domain, crs = CRS("+proj=cea +units=km"))
##             }
##         }
##         domain_crs <- inlabru:::fm_ensure_crs(domain$crs)
##         if (identical(int.args[["poly_method"]], "legacy")) {
##             ips <- inlabru:::int.polygon(domain, loc = polyloc[, 1:2], 
##                 group = polyloc[, 3], method = int.args$method, 
##                 nsub = int.args$nsub2)
##         }
##         else {
##             ips <- inlabru:::bru_int_polygon(domain, poly_segm, method = int.args$method, 
##                 nsub = int.args$nsub2)
##         }
##         df <- data.frame(samplers@data[ips$group, group, drop = FALSE], 
##             weight = ips[, "weight"] * samplers@data[ips$group, 
##                 "weight"])
##         ips <- SpatialPointsDataFrame(ips[, c("x", "y")], data = df, 
##             match.ID = FALSE, proj4string = domain_crs)
##         if (!inlabru:::fm_crs_is_null(domain_crs) && !inlabru:::fm_crs_is_null(samplers_crs)) {
##             ips <- stransform(ips, crs = samplers_crs)
##         }
##         coord_names <- c("x", "y", "coordinateZ")
##         if (!is.null(coordnames(samplers))) {
##             coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
##         }
##         else if (!is.null(name)) {
##             coord_names[seq_along(name)] <- name
##         }
##         coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
##     }
##     else {
##         stop("No integration handling code reached; please notify the package developer.")
##     }
##     ips
## }

## int.polygon <- function (mesh, loc, group = NULL, method = NULL, ...) {
##     if (is.null(group)) {
##         group <- rep(1, nrow(loc))
##     }
##     method <- match.arg(method, c("stable", "direct"))
##     ipsl <- list()
##     for (g in unique(group)) {
##         gloc <- loc[group == g, , drop = FALSE]
##         bnd <- INLA::inla.mesh.segment(loc = gloc, is.bnd = TRUE)
##         integ <- make_stable_integration_points(mesh, bnd, ...)
##         if (method %in% c("stable")) {
##             integ <- integration_weight_aggregation(mesh, integ)
##         }
##         ok <- INLA::inla.mesh.project(mesh, integ$loc)$ok & (integ$weight > 
##             0)
##         ips <- data.frame(integ$loc[ok, 1:2, drop = FALSE])
##         colnames(ips) <- c("x", "y")
##         ips$weight <- integ$weight[ok]
##         ips$group <- rep(g, nrow(ips))
##         ipsl <- c(ipsl, list(ips))
##     }
##     do.call(rbind, ipsl)
## }


## bru_int_polygon <- function (mesh, polylist, method = NULL, ...){
##     method <- match.arg(method, c("stable", "direct"))
##     ipsl <- list()
##     for (g in seq_along(polylist)) {
##         poly <- polylist[[g]]
##         integ <- make_stable_integration_points(mesh, poly, ...)
##         ok <- INLA::inla.mesh.project(mesh, integ$loc)$ok & (integ$weight > 
##             0)
##         if (any(ok)) {
##             integ <- list(loc = integ$loc[ok, 1:2, drop = FALSE], 
##                 weight = integ$weight[ok])
##             if (method %in% c("stable")) {
##                 integ <- integration_weight_aggregation(mesh, 
##                   integ)
##             }
##             ips <- data.frame(x = integ$loc[, 1], y = integ$loc[, 
##                 2], weight = integ$weight, group = g)
##             ipsl <- c(ipsl, list(ips))
##         }
##     }
##     do.call(rbind, ipsl)
## }


## make_stable_integration_points <- function (mesh, bnd, nsub = NULL) 
## {
##     if (is.null(nsub)) {
##         nsub <- 9
##     }
##     stopifnot(nsub >= 0)
##     nB <- (nsub + 1)^2
##     b <- seq(1/3, 1/3 + nsub, length = nsub + 1)/(nsub + 1)
##     bb <- as.matrix(expand.grid(b, b))
##     refl <- rowSums(bb) > 1
##     if (any(refl)) {
##         bb[refl, ] <- cbind(1 - bb[refl, 2], 1 - bb[refl, 1])
##     }
##     barycentric_grid <- cbind(1 - rowSums(bb), bb)
##     nT <- nrow(mesh$graph$tv)
##     loc <- matrix(0, nT * nB, 3)
##     idx_end <- 0
##     for (tri in seq_len(nT)) {
##         idx_start <- idx_end + 1
##         idx_end <- idx_start + nB - 1
##         loc[seq(idx_start, idx_end, length = nB), ] <- as.matrix(barycentric_grid %*% 
##             mesh$loc[mesh$graph$tv[tri, ], , drop = FALSE])
##     }
##     weight <- rep(INLA::inla.mesh.fem(mesh, order = 1)$ta/nB, 
##         each = nB)
##     mesh_bnd <- INLA::inla.mesh.create(boundary = bnd)
##     ok <- INLA::inla.mesh.projector(mesh_bnd, loc = loc)$proj$ok
##     list(loc = loc[ok, , drop = FALSE], weight = weight[ok])
## }
