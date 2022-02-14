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

## --------------
## link functions
## --------------

## oscillating parameter 
theta.2.phi <- function(theta){
    (1-exp(-theta))/(1+exp(-theta))
}
## 
theta.2.sigma <- function(theta){
    exp(theta)
}
##
theta.2.rho <- function(theta){
    exp(theta)
}

## helper functions for computing the
## posterior distribution of gridness score
## input is a fitted model with inlabru (e.g. fit.space and fit.space.direction)
## object returns attributes for summaries of the posterior
posterior.spatial.gridness.score <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.phi, inlabru.fitted.object$marginals.hyperpar[["Theta3 for spde2"]])
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
    marg                  <- inla.tmarginal(theta.2.sigma, inlabru.fitted.object$marginals.hyperpar[["Theta5 for spde2"]])
    summaries             <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval          <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}
posterior.directional.range <- function(inlabru.fitted.object){
    marg                  <- inla.tmarginal(theta.2.rho, inlabru.fitted.object$marginals.hyperpar[["Theta4 for spde2"]])
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
    ## df.indices labels the indices of the knots for the directional, the spatial and the spatio-directional basis functions
    ## this data frame is created to create correspondences
    ## between spatio-directional basis knots with spatial basis knots, an
    ## between spatio-directional basis knots with head directional basis knots, respectively.
    df.indices <- data.frame(dir = sort(rep(1:mesh.hd$n, mesh$n)), space = rep(1:mesh$n, mesh.hd$n), cross = 1:(mesh$n*mesh.hd$n))
    ## So for example, if the spatio-directional basis knots are labeled as 1, 2, ..., p_Omega * p_Theta
    ## then the function mapindex2space.direction_basis takes as argument the label of spatio-diretional basis knot
    ## and returns the coordinates and the head direction associated with the spatial basis function and the
    ## head directional basis function. This function uses mapindex2space.direction_basis which works similarly but
    ## instead of returning coords and angles, it returns the indices of the associated basis functions.
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
    Ypos.tmp <- Ypos.tmp %>% mutate(HD.split = map2(hd, hd.lead, split.arcs),
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
                                        #
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
                    summarize(val = sum(max(dGamma.trap*val.M0, tol))/2,
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
                    summarize(val = sum(max(dGamma.trap*val.M1, tol))/2,
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


