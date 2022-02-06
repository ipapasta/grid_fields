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







if(FALSE){
    ipoints <- 
        function (samplers = NULL, domain = NULL, name = NULL, group = NULL, 
                  int.args = NULL, project = NULL) 
    {
        int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
        if (is.null(int.args)) {
            int.args <- list()
        }
        missing.args <- setdiff(names(int.args.default), names(int.args))
        int.args[missing.args] <- int.args.default[missing.args]
        if (!is.null(int.args[["nsub"]])) {
            int.args[["nsub1"]] <- int.args[["nsub"]]
        }
        if (!is.null(int.args[["nsub"]])) {
            int.args[["nsub2"]] <- int.args[["nsub"]]
        }
        if (!is.null(project)) {
            if (project && !identical(int.args$method, "stable")) {
                stop("ipoints(project=TRUE) is deprecated, and int.args$methods != 'stable'")
            }
            else if (!project && identical(int.args$method, "stable")) {
                stop("ipoints(project=FALSE) is deprecated, and int.args$methods == 'stable'")
            }
            warning("ipoints(project=", ifelse(project, "TRUE", "FALSE"), 
                    ") is deprecated. Will use int.args$method = '", 
                    int.args[["method"]], "' instead.")
        }
        if (is.null(domain) && inherits(samplers, c("inla.mesh.1d", 
                                                    "inla.mesh"))) {
            domain <- samplers
            samplers <- NULL
        }
        is_2d <- (!is.null(samplers) && inherits(samplers, c("SpatialPoints", 
                                                             "SpatialPointsDataFrame", "SpatialPolygons", "SpatialPolygonsDataFrame", 
                                                             "SpatialLines", "SpatialLinesDataFrame"))) || inherits(domain, 
                                                                                                                    "inla.mesh")
        is_1d <- !is_2d && ((!is.null(samplers) && is.numeric(samplers)) || 
                            (!is.null(domain) && (is.numeric(domain) || inherits(domain, 
                                                                                 "inla.mesh.1d"))))
        if (!is_1d && !is_2d) {
            stop("Unable to determine integration domain definition")
        }
        if (is_1d && !is.null(samplers) && !is.null(domain) && is.numeric(domain) && 
            length(domain) == 1) {
            int.args[["nsub1"]] <- domain
            domain <- NULL
            int.args[["method"]] <- "direct"
        }
        if (is_2d && !is.null(samplers) && !is.null(domain) && is.numeric(domain) && 
            length(domain) == 1) {
            int.args[["nsub2"]] <- domain
            domain <- NULL
            int.args[["method"]] <- "direct"
        }
        if (is.null(domain) && inherits(samplers, c("inla.mesh.1d", 
                                                    "inla.mesh"))) {
            domain <- samplers
            samplers <- NULL
        }
        if (is_1d && is.null(name)) {
            name <- "x"
        }
        pregroup <- NULL
        if (is.data.frame(samplers)) {
            if (!("weight" %in% names(samplers))) {
                samplers$weight <- 1
            }
            ips <- samplers
        }
        else if (is_1d && is.null(samplers) && is.numeric(domain)) {
            ips <- data.frame(x = as.vector(domain), weight = 1)
            colnames(ips) <- c(name, "weight")
        }
        else if (is_1d && is.null(domain) && is.integer(samplers)) {
            ips <- data.frame(x = as.vector(samplers), weight = 1)
            colnames(ips) <- c(name, "weight")
        }
        else if (is_1d && is.null(samplers) && inherits(domain, "inla.mesh.1d") && 
                 identical(int.args[["method"]], "stable")) {
            ips <- data.frame(x = domain$loc, weight = Matrix::diag(INLA::inla.mesh.fem(domain)$c0))
            colnames(ips) <- c(name, "weight")
        }
        else if (is_1d) {
            domain_range <- if (inherits(domain, "inla.mesh.1d")) {
                                domain$interval
                            }
                            else {
                                NULL
                            }
            if (is.null(samplers)) {
                samplers <- matrix(domain_range, 1, 2)
            }
            else {
                if (is.null(dim(samplers))) {
                    samplers <- matrix(samplers, nrow = 1)
                }
                if (ncol(samplers) != 2) {
                    stop("Interval description matrix must have 2 elements or be a 2-column matrix.")
                }
                if (is.null(domain)) {
                    domain <- INLA::inla.mesh.1d(sort(unique(as.vector(samplers))))
                }
            }
            ips <- list()
            if (domain$degree >= 2) {
                warning("Integration points projected onto knots may lead to instability for degree >= 2 basis functions.")
            }
            nsub <- int.args[["nsub1"]]
            u <- rep((seq_len(nsub) - 0.5)/nsub, domain$n - 1)
            int_loc <- domain$loc[rep(seq_len(domain$n - 1), each = nsub)] * 
                (1 - u) + domain$loc[rep(seq_len(domain$n - 1) + 
                                         1, each = nsub)] * u
            int_w <- (domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] - 
                      domain$loc[rep(seq_len(domain$n - 1), each = nsub)])/nsub
            for (j in 1:nrow(samplers)) {
                subsamplers <- samplers[j, ]
                if (identical(int.args[["method"]], "stable")) {
                    A_w <- INLA::inla.spde.make.A(domain, int_loc, 
                                                  weights = int_w * (int_loc >= min(subsamplers)) * 
                                                      (int_loc <= max(subsamplers)))
                    ips[[j]] <- data.frame(loc = domain$loc, weight = Matrix::colSums(A_w))
                }
                else {
                    inside <- (int_loc >= min(subsamplers)) & (int_loc <= 
                                                               max(subsamplers))
                    ips[[j]] <- data.frame(loc = int_loc[inside], 
                                           weight = int_w[inside])
                }
                colnames(ips[[j]]) <- c(name, "weight")
            }
            ips <- do.call(rbind, ips)
        }
        else if (inherits(domain, "inla.mesh") && is.null(samplers) && 
                 identical(int.args[["method"]], "stable")) {
            coord_names <- c("x", "y", "coordinateZ")
            if (!is.null(name)) {
                coord_names[seq_along(name)] <- name
            }
            if (!fm_crs_is_null(domain$crs)) {
                crs <- domain$crs
                samplers <- stransform(domain, crs = CRS("+proj=cea +units=km"))
            }
            ips <- vertices(domain)
            ips$weight <- INLA::inla.mesh.fem(domain, order = 1)$va
            if (!fm_crs_is_null(domain$crs)) {
                ips <- stransform(ips, crs = crs)
            }
            coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
        }
        else if (class(samplers) == "SpatialPoints") {
            ips <- samplers
            ips$weight <- 1
        }
        else if (class(samplers) == "SpatialPointsDataFrame") {
            if (!("weight" %in% names(samplers))) {
                warning("The integration points provided have no weight column. Setting weights to 1.")
                samplers$weight <- 1
            }
            ips <- samplers
        }
        else if (inherits(samplers, "SpatialLines") || inherits(samplers, 
                                                                "SpatialLinesDataFrame")) {
            if (inherits(samplers, "SpatialLines") && !inherits(samplers, 
                                                                "SpatialLinesDataFrame")) {
                samplers <- SpatialLinesDataFrame(samplers, data = data.frame(weight = rep(1, 
                                                                                           length(samplers))))
            }
            if (!("weight" %in% names(samplers))) {
                samplers$weight <- 1
            }
            ips <- int.slines(samplers, domain, group = group, project = identical(int.args[["method"]], 
                                                                                   "stable"))
            coord_names <- c("x", "y", "coordinateZ")
            if (!is.null(coordnames(samplers))) {
                coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
            }
            else if (!is.null(name)) {
                coord_names[seq_along(name)] <- name
            }
            coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
        }
        else if (is_2d && (inherits(samplers, c("SpatialPolygons", 
                                                "SpatialPolygonsDataFrame")) || is.null(samplers))) {
            if (is.null(samplers)) {
                stop("Direct integration scheme for mesh domain with no samplers is not yet implemented.")
            }
            if (class(samplers)[1] == "SpatialPolygons") {
                samplers <- SpatialPolygonsDataFrame(samplers, data = data.frame(weight = rep(1, 
                                                                                              length(samplers))), match.ID = FALSE)
            }
            else if (is.null(samplers@data[["weight"]])) {
                samplers@data[["weight"]] <- 1
            }
            cnames <- coordnames(samplers)
            samplers_crs <- fm_sp_get_crs(samplers)
            if (!fm_crs_is_null(domain$crs)) {
                samplers <- stransform(samplers, crs = sp::CRS("+proj=cea +units=km"))
            }
            polyloc <- do.call(rbind, lapply(1:length(samplers), 
                                             function(k) {
                                                 cbind(x = rev(coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 
                                                                                                                 1]), y = rev(coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 
                                                                                                                                                                                2]), group = k)
                                             }))
            poly_segm <- INLA::inla.sp2segment(samplers, join = FALSE)
            poly_segm <- lapply(seq_along(poly_segm), function(k) {
                segm <- poly_segm[[k]]
                segm[["grp"]] <- rep(k, NROW(segm[["idx"]]))
                segm[["is.bnd"]] <- TRUE
                segm
            })
            if (is.null(domain)) {
                warning("Computing integration points from polygon; specify a mesh for better numerical control.")
                max.edge <- max(diff(range(polyloc[, 1])), diff(range(polyloc[, 
                                                                              2])))/20
                domain <- INLA::inla.mesh.2d(boundary = samplers, 
                                             max.edge = max.edge)
                domain$crs <- fm_sp_get_crs(samplers)
            }
            else {
                if (!fm_crs_is_null(domain$crs)) {
                    domain <- stransform(domain, crs = CRS("+proj=cea +units=km"))
                }
            }
            domain_crs <- fm_ensure_crs(domain$crs)
            if (identical(int.args[["poly_method"]], "legacy")) {
                ips <- int.polygon(domain, loc = polyloc[, 1:2], 
                                   group = polyloc[, 3], method = int.args$method, 
                                   nsub = int.args$nsub2)
            }
            else {
                ips <- bru_int_polygon(domain, poly_segm, method = int.args$method, 
                                       nsub = int.args$nsub2)
            }
            df <- data.frame(samplers@data[ips$group, pregroup, drop = FALSE], 
                             weight = ips[, "weight"] * samplers@data[ips$group, 
                                                                      "weight"])
            ips <- SpatialPointsDataFrame(ips[, c("x", "y")], data = df, 
                                          match.ID = FALSE, proj4string = domain_crs)
            if (!fm_crs_is_null(domain_crs) && !fm_crs_is_null(samplers_crs)) {
                ips <- stransform(ips, crs = samplers_crs)
            }
            coord_names <- c("x", "y", "coordinateZ")
            if (!is.null(coordnames(samplers))) {
                coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
            }
            else if (!is.null(name)) {
                coord_names[seq_along(name)] <- name
            }
            coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
        }
        else {
            stop("No integration handling code reached; please notify the package developer.")
        }
        ips
    }





    spde.posterior <- function (result, name, what = "range") 
    {
        spdespec <- result$bru_info$model$effects[[name]]$main$model
        spderesult <- INLA::inla.spde.result(result, name, spdespec)
        if (what == "matern.correlation" || what == "matern.covariance") {
            xmax <- exp(spderesult$summary.log.range.nominal[["0.975quant"]]) * 
                1.2
            x <- seq(0, xmax, length = 200)
            log.range <- list(mean = spderesult$summary.log.range.nominal[["mean"]], 
                              sd = spderesult$summary.log.range.nominal[["sd"]])
            log.variance <- list(mean = spderesult$summary.log.variance.nominal[["mean"]], 
                                 sd = spderesult$summary.log.variance.nominal[["sd"]])
            if (what == "matern.correlation") {
                corr <- TRUE
                ylab <- "Matern Correlation"
                out <- materncov.bands(result$bru_info$model$effects[[name]]$main$mapper$mesh, 
                                       dist = x, log.range = log.range, log.variance = NULL, 
                                       alpha = 2, quantile = 0.95)
            }
            else {
                corr <- FALSE
                ylab <- "Matern Covariance"
                out <- materncov.bands(result$bru_info$model$effects[[name]]$main$mapper$mesh, 
                                       dist = x, log.range = log.range, log.variance = log.variance, 
                                       alpha = 2, quantile = 0.95)
            }
            df <- data.frame(x = x, median = out$median, q0.025 = out$lower, 
                             q0.975 = out$upper)
            attr(df, "type") <- "1d"
            attr(df, "misc") <- list(dims = "x", predictor = c("distance", 
                                                               ylab))
            class(df) <- list("prediction", "data.frame")
            df
        }
        else {
            marg <- switch(what, range = spderesult$marginals.range.nominal[[1]], 
                           log.range = spderesult$marginals.log.range.nominal[[1]], 
                           variance = spderesult$marginals.variance.nominal[[1]], 
                           log.variance = spderesult$marginals.log.variance.nominal[[1]])
            if (is.null(marg)) 
                stop("Invalid varname: ", what, ". must be one of 'range', \n                           'log.range',  'variance',  'log.variance', \n                           'matern.correlation', matern.covariance")
            med <- INLA::inla.qmarginal(0.5, marg)
            uq <- INLA::inla.qmarginal(0.975, marg)
            lq <- INLA::inla.qmarginal(0.025, marg)
            inner.x <- seq(lq, uq, length.out = 100)
            inner.marg <- data.frame(x = inner.x, y = INLA::inla.dmarginal(inner.x, 
                                                                           marg))
            colnames(inner.marg) <- c(what, "pdf")
            df <- data.frame(marg)
            colnames(df) <- c(what, "pdf")
            attr(df, "type") <- "0d"
            attr(df, "summary") <- list(uq = uq, lq = lq, median = med, 
                                        inner.marg = inner.marg)
            class(df) <- list("prediction", "data.frame")
            df
        }
    }

}

