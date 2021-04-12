## 
## seed for reproducibility
##
## set.seed(111086) 
## !! quilt.plot
## 
## load packages
##
counter <- 0
## sim     <- FALSE
## library(Rcpp)
library(tidyverse)
library(purrr)
library(INLA)   #install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
inla.setOption(pardiso.license = "/Users/ipapasta/pardiso.lic")
library(inlabru)
library(sp)
library(fields)
library(nloptr)
## require(rgdal, quietly=TRUE)
## scp -r /home/ipapasta/Dropbox/org/Research/Software/R/grid_fields/R/* ipapasta@xserver2.maths.ed.ac.uk:/home/ipapasta/Software/R/grid_fields/R/
## 
## source R and cpp functions.
## 
source("load_data.R")
source("Functions.R")
source("osc_precision.R")
source("hd_precision.R")
source("objective.R")
source("temp_precision.R")
source("priorbetaXZ_osc_temp.R")
source("priortheta_osc_temp.R")              
source("gradient_osc_temp.R")
source("hessian_osc_temp.R")
source("llik.R")
source("marginalposterior.R")

k    <- 3
mesh      <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
## 
fem.mesh  <- inla.mesh.fem(mesh, order = 2)
Bspatial.phi0 = matrix(c(0,1,0,0), nrow=1)
Bspatial.phi1 = matrix(c(0,0,1,0), nrow=1)
Bspatial.phi2 = matrix(c(0,0,0,1), nrow=1)
M0.spatial = fem.mesh$c0 # C
M1.spatial = fem.mesh$g1
M2.spatial = fem.mesh$g2
## 
p         <- mesh$n                         
theta.nodes <- seq(0, 2*pi, len=20)
mesh.hd     <- inla.mesh.1d(theta.nodes, boundary="cyclic", degree=1)
fem.mesh.hd <- inla.mesh.fem(mesh.hd, order = 2)
##
Bhd.phi0 = matrix(c(0,1,0), nrow=1)
Bhd.phi1 = matrix(c(0,0,1), nrow=1)
M0.hd = fem.mesh.hd$c0
M1.hd = fem.mesh.hd$g1
M2.hd = fem.mesh.hd$g2
## 
nodes       <- c(mesh.hd$loc, 2*pi)
intervals   <- head(cbind(nodes, lead(nodes)), -1)

Ypos.tmp <- data.frame(
    hd=X$hd, time=X$synced_time,
    coords=I(lapply(as.list(apply(cbind(X$position_x, X$position_y),1, as.list)), unlist))) %>%
    mutate(coords.lead = lead(coords)) %>%
    mutate(time.lead = lead(X$synced_time)) %>%
    mutate(hd.lead = lead(X$hd)) %>%
    head(-1)

## options(warn=2)


## o <- lapply(1:nrow(Ypos.tmp), function(k) Ypos.tmp$hd == Ypos.tmp$hd.lead)


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

Ypos.tmp <- Ypos.tmp %>% dplyr::select(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead)%>%
    unnest(cols=c(new.time, new.time.lead, new.hd, new.hd.lead, new.coords, new.coords.lead))

names(Ypos.tmp) <- c("time", "time.lead", "hd", "hd.lead", "coords", "coords.lead")
## Ypos.tmp$Time.split[3][[1]]
    

line.segments <- split.lines(mesh, sp=Ypos.tmp$coords,
                             filter.zero.length=FALSE,
                             ep=Ypos.tmp$coords.lead, tol=.0)

## line.segments.filtered <- split.lines(mesh, sp=do.call("rbind", Ypos.tmp$coords),
##                              filter.zero.length=TRUE,
##                              ep=do.call("rbind",Ypos.tmp$coords.lead), tol=1e-4)
## dim(line.segments.filtered$sp)
## filter.index contains the # of lines used in the integration

df <- data.frame(origin=line.segments$split.origin,
                 filter.index=line.segments$filter.index,
                 sp=I(lapply(as.list(apply(line.segments$sp, 1, as.list)), unlist)),
                 ep=I(lapply(as.list(apply(line.segments$ep, 1, as.list)), unlist))) %>%
    group_by(origin) %>%
    summarize(sp=list(sp), ep=list(ep), filter.index=list(filter.index)) %>%
    mutate(sp = lapply(sp, function(x) do.call("rbind", x))) %>%
    mutate(ep = lapply(ep, function(x) do.call("rbind", x))) 




## length of line segments
## duration of time intervals + data as attribute
## change in head direction + data as attribute


## oo <- inner_join(Ypos.tmp %>%
##                    mutate(origin=1:nrow(Ypos.tmp)), df)


Ypos <- inner_join(Ypos.tmp %>%
                   mutate(origin=1:nrow(Ypos.tmp)), df) %>%    
    mutate(Li = map2(sp, ep,
                     function(x, y) apply(y-x, 1, function(z) sqrt(sum(z^2))))) %>%  
    mutate(Ti = pmap(list(time, time.lead, Li), interpolate)) %>%
    mutate(HDi = pmap(list(hd, hd.lead, Li), interpolate )) 



filter.index  <- do.call("c", Ypos$filter.index)
## no filter
coords.trap  <- rbind(do.call("rbind",Ypos$sp)[filter.index,], tail(do.call("rbind",Ypos$ep),1))
HD.data      <- c(do.call("c", (Ypos %>% mutate(HD=lapply(HDi, function(x) attr(x, "data"))))$HD), tail(Ypos$hd.lead, 1))
T.data       <- c(do.call("c", (Ypos %>% mutate(T=lapply(Ti, function(x) attr(x, "data"))))$T), tail(Ypos$time.lead, 1))



## ---------------------------------------
## SpatialPointsDataFrame and SpatialLines
## ---------------------------------------
Y.spdf    <- SpatialPointsDataFrame(SpatialPoints(cbind(Y$position_x, Y$position_y)), as.data.frame(Y%>%select(-c(position_x, position_y))))
Ypos.sldf <- SpatialLinesDataFrame(SpatialLines(lapply(as.list(1:nrow(Ypos)),function(k) Lines(list(Line(cbind(c(Ypos$coords[k,1],
                                                                                    Ypos$coords.lead[k,1]),
                                                                                  c(Ypos$coords[k,2],
                                                                                    Ypos$coords.lead[k,2])))),
                                                                                  ID=k))),
                                   Ypos %>% select(-c(coords, coords.lead)))

data <- list(Ypos=Ypos, Y=Y, Yspdf=Y.spdf, Ypos.sldf = Ypos.sldf)




## SpatialPointsDataFrame and SpatialLines



## dim(coords.trap)
## length(functions.multiplicity)
      



## circular domain and temporal domain meshes
mesh1d  <- inla.mesh.1d(loc=T.data[seq(1, length(T.data), by = 300)], order=2)
print(paste("mesh1d:", head(diff(mesh1d$loc))))

## POSITIONAL - Used for INTEGRAL
Atilde <- inla.mesh.projector(mesh1d, loc=T.data)$proj$A
Aosc   <- inla.mesh.projector(mesh, loc=coords.trap)$proj$A
Ahd    <- inla.mesh.projector(mesh.hd, loc=HD.data)$proj$A
A      <- inla.row.kron(Ahd, Aosc)


## OBSERVED
Aosc.obs  <- inla.spde.make.A(mesh=mesh,
                              loc=as.matrix(data$Y %>%
                                            dplyr:: select(position_x, position_y)))
Ahd.obs   <- inla.spde.make.A(mesh=mesh.hd, data$Y$hd)

Aobs      <- inla.row.kron(Ahd.obs, Aosc.obs)
Atildeobs    <- inla.spde.make.A(mesh=mesh1d, data$Y$firing_times)

## c(3.1305657, -0.8904750, -3.3000715,  1.2177724,  1.1621642,  3.0267431, -0.2246174)




## nrow(A)==nrow(At); 92264 each row of above matrices contains
## non-zero values at knots wrapping a distinct line segment.
## Ck <- sapply(dGamma, function(x) rep(x, 6))
dGamma <- c(do.call("c", Ypos$Li))
dT  <- diff(T.data)
Atmp <- as(A, "dgTMatrix")
A.indices <- cbind(cbind(Atmp@i+1, Atmp@j+1), Atmp@x)
A.indices <- A.indices[order(A.indices[,1]),] %>% as.data.frame
Attmp <- as(Atilde, "dgTMatrix")
At.indices <- cbind(cbind(Attmp@i+1, Attmp@j+1), Attmp@x)
At.indices <- At.indices[order(At.indices[,1]),] %>% as.data.frame
names(A.indices) <- c("tk", "i", "psi.ot") #ot: omega x theta
names(At.indices) <- c("tk", "l", "psi.t")

df <- full_join(At.indices %>% group_by(tk) %>% nest(),
                A.indices %>% group_by(tk) %>% nest(), by="tk") %>%
    ungroup %>% 
    mutate(dGamma=c(dGamma,0),
           dGamma.lead = lead(dGamma),
           dGamma.lag = lag(dGamma),
           val = pmap(list(data.x, data.y, dGamma), function(x, y, z){
        oo <- expand.grid(1:nrow(x), 1:nrow(y))
        ooo <- unlist(lapply(1:(nrow(x) * nrow(y)), function(k) {
            x$psi.t[oo[k,1]] * y$psi.ot[oo[k,2]]}))
        oooo <- data.frame(l=x$l[oo[,1]], i=y$i[oo[,2]],val=ooo)
        oooo
    })) %>% dplyr::select(-c("data.x", "data.y")) %>%
    unnest(val)

df.W <- rbind(df %>% mutate(group=tk,
                            dGamma.lag=0),
              df %>% filter(tk!=1) %>%
              mutate(group=tk-1,
                     dGamma=0)) %>%
    arrange(group) %>%
    mutate(dGamma.trap = dGamma + dGamma.lag)

tol <- 0
df.dGamma.sum.k.kplus1 <- df.W %>% group_by(group, l, i) %>%
    summarize(val = sum(max(dGamma.trap*val, tol)))
## df.dGamma.sum.k.kplus1$val[df.dGamma.sum.k.kplus1$val==0] <- 1e-16

## the integration scheme is stable, which can be achieved by ensuring
## positive weights on all basis functions that interact with the
## line/curve of integration.

W <- sparseMatrix(i=df.dGamma.sum.k.kplus1$l,
                  j=df.dGamma.sum.k.kplus1$i,
                  x=df.dGamma.sum.k.kplus1$val/2)

W         <- W %>% cbind(Matrix(0, nrow=nrow(W), ncol=ncol(A)-ncol(W)))


if(FALSE){
    Wtmp <- as(W, "dgTMatrix")
    Wtmp.triplets <- cbind(cbind(Wtmp@i+1, Wtmp@j+1), Wtmp@x)
    ## 
    quilt.plot(Wtmp.triplets[,1], Wtmp.triplets[,2], Wtmp.triplets[,3],
               col=ocean.balance(100), asp=1/100)
}



load("modelfit.RData")
## par.theta <- c(log(27),
##                log(0.2),
##                -log((1-(-0.95))/(1+(-0.95))),
##                log(4.3),
##                log(1),
##                log(5),
##                log(9,5))

## par.theta <- opt.theta.no.prior$par

Xest       <-   Xinit    <- Matrix(rep(0, ncol(A)), ncol=1)
Zest       <-   Zinit    <- Matrix(rep(0, mesh1d$n), ncol=1)

betaest       <-   betainit <- 0
## nrow(Y)/sum((Ypos %>% mutate(speed.lead = lead(speed), dt=c(diff(time)[1], diff(time) ))  %>%
##            head(-1) %>%
##            mutate(val=dt*((speed + speed.lead)/2)))$val)
                                        #(number of firing events)/
                                        #(\int_\Gamma ds)
gradest       <- NULL
Hessianest    <- NULL


## debugonce(pthetapc.prop.marg.post_osc_temp)
## debugonce(grad.objective_osc_temp)
## debugonce(hessian.objective_osc_temp)
## ldmvnorm(hessian.objective_osc_temp)

set.seed(111086)




## opt.theta.one <- optim(par=opt.theta$par, fn = pthetapc.prop.marg.post_osc_temp, data=data, 
##                        X=Xinit, Z=Zinit, beta=betaest,
##                        mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
##                        A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                        W=W, acc=1e-3,
##                        print.verbose=FALSE,
##                        control=list(maxit=1))

par.theta <- c(log(21.77),
               log(0.48),
               -log((1-(-0.91))/(1+(-0.91))),
               log(3.23),
               log(3.01),
               log(23.36),
               log(0.79))



opt.theta <- optim(par=par.theta, fn = pthetapc.prop.marg.post_osc_temp, data=data, 
                   X=Xinit, Z=Zinit, beta=betainit,
                   mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
                   A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
                   W=W, acc=1e-2,
                   print.verbose=TRUE,
                   control=list(maxit=5000))


## opt.theta <- optim(par=opt.theta$par, fn = pthetapc.prop.marg.post_osc_temp, data=data, 
##                    X=Xinit, Z=Zinit, beta=betainit,
##                    mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
##                    A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                    W=W, acc=1e-3,
##                    print.verbose=TRUE,method="BFGS")

save(opt.theta, Xest, Zest, betaest, Hessianest, file="fitted_model.RData")



## 
## inlabru implementation
##

if(FALSE){
    library(inlabru)
    B.phi0.matern = matrix(c(0,1,0), nrow=1)
    B.phi1.matern = matrix(c(0,0,1), nrow=1)
    B.phi0.oscillating = matrix(c(0,1,0,0), nrow=1)
    B.phi1.oscillating = matrix(c(0,0,1,0), nrow=1)
    B.phi2.oscillating = matrix(c(0,0,0,1), nrow=1)
    fem.mesh <- inla.mesh.fem(mesh, order = 2)
    fem.mesh.hd <- inla.mesh.fem(mesh.hd, order = 2)
    M0 = fem.mesh$c0
    M1 = fem.mesh$g1
    M2 = fem.mesh$g2
    M0.hd = fem.mesh.hd$c0
    M1.hd = fem.mesh.hd$g1
    M2.hd = fem.mesh.hd$g2
    
    source("rgeneric_models.R")
    
    ## 
    ## rgeneric
    oscillating.rgeneric <- inla.rgeneric.define(oscillating.model, M = list(M0=M0, M1=M1, M2=M2))
    circular.rgeneric <- inla.rgeneric.define(circular1D.model, M=list(M0=M0.hd, M1=M1.hd, M2=M2.hd))
    ## cmp.oscillating.rgeneric <- hd2 + coordinates ~ circular(hd2, model=circular.rgeneric) +
    ##     mySPDE(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
    cmp.oscillating.rgeneric <- coordinates ~ 
        spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
        + Intercept
    
    ## +        direction(Y.spdf$hd)
    ## hd 
    ## form <- coordinates  ~ spde2 + hd + Intercept 
    fit.oscillating.rgeneric <- lgcp(cmp.oscillating.rgeneric, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh),
                                     ## formula=form,
                                     options=list(verbose = TRUE))

    ## -------------------------------------
    ## temporal
    ## -------------------------------------
    cmp.temporal <- firing_times ~ 
        temporal(firing_times, model = "ar1", mapper=bru_mapper(mesh1d, indexed=TRUE)) +
        + Intercept
    fit.temporal <- lgcp(cmp.temporal, data = Y.spdf, samplers = Ypos.sldf, domain = list(firing_times = mesh1d),
                                     ## formula=form,
                         options=list(verbose = TRUE))

    ## -------------------------------------
    ## directional
    ## -------------------------------------
    circular.rgeneric <- inla.rgeneric.define(circular1D.model, M=list(M0=M0.hd, M1=M1.hd, M2=M2.hd))
    cmp.circular <- hd ~ 
        directional(hd, model = circular.rgeneric, mapper=bru_mapper(mesh.hd, indexed=TRUE)) +
        + Intercept

    fit.circular <- lgcp(cmp.circular, data = Y.spdf, samplers = Ypos.sldf, domain = list(hd = mesh.hd),
                         options=list(verbose = TRUE))

    ## example
    df.ex <- tibble(x=rnorm(100), y = 0.8 * x + sqrt(1-.8^2)*rnorm(100), th = atan2(y,x) + pi)
    A.hd <- inla.spde.make.A(mesh = mesh.hd, loc = df.ex$th)
    indexs <- inla.spde.make.index("dir", mesh.hd$n)
    stk.e <- inla.stack(
        tag = "est",
        data = list(y = df.ex$y),
        A = list(1, A.hd),
        effects = list(data.frame(b0 = rep(1, nrow(df.ex))), s = indexs)
    )
    formula <- y ~ 0 + b0 + f(dir, model = circular.rgeneric)
    inla(formula, family="normal",   data = inla.stack.data(stk.e),
         control.predictor = list(
             compute = TRUE,
             A = inla.stack.A(stk.e)
         ))
    
    ## -------------------------------------
    ## oscillating and temporal - incorrect
    ## -------------------------------------
    cmp.oscillating.temporal <- coordinates + firing_times ~
        spde2(coordinates, model = oscillating.rgeneric, mapper=bru_mapper(mesh, indexed=TRUE)) +
        f(firing_times, model = "ar1", mapper=bru_mapper(mesh1d, indexed=TRUE)) +
        + Intercept
    fit.oscillating.temporal <- lgcp(cmp.oscillating.temporal, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh, firing_times = mesh1d),
                                     ## formula=form,
                                     options=list(verbose = TRUE))    ## 
    ## 
    
    plot(fit.oscillating.rgeneric$marginals.hyperpar[[1]], type="l", xlim=c(-7,10))
    lines(fit.oscillating.rgeneric$marginals.hyperpar[[2]], lty=2)
    lines(fit.oscillating.rgeneric$marginals.hyperpar[[3]], lty=2)
    
    pxl <- pixels(mesh, nx=1000, ny=1000)
    pr.int <- predict(fit.oscillating.rgeneric, pxl, ~ mySPDE)

    library(pals)
    ggplot() +
        gg(pr.int) +
        gg(mycoords, color="red", size=0.2)+
        scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar")+
        xlim(0,100)+
        ylim(0,100)+
        ## gg(mexdolphin$ppoly) +
        ## gg(mexdolphin$samplers, color = "grey") +
        ## gg(mexdolphin$points, size = 0.2, alpha = 1) +
        ## noyticks +
        ## noxticks +
        ## theme(legend.key.width = unit(x = 0.2, "cm"), legend.key.height = unit(x = 0.3, "cm")) +
        ## theme(legend.text = element_text(size = 6)) +
        ## guides(fill = FALSE) +
        coord_equal() + theme_classic()
    ## 
    
    matern.spde2 <- inla.spde2.pcmatern(mesh,
                                        prior.sigma = c(2, 0.01),
                                        prior.range = c(20, 0.01))

    ## fit.oscillating.rgeneric.no.samplers <- lgcp(cmp.oscillating.rgeneric,
    ## data = Y.spdf, samplers = NULL, domain = list(coordinates = mesh), options=list(verbose = TRUE))    
    ## pr.int <- predict(fit.oscillating.rgeneric.no.samplers, pxl, ~ mySPDE + Intercept)
    
    ## 
    ## cmp.matern <- coordinates ~ mySPDE(coordinates, model = matern.spde2, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
    ## fit.matern <- lgcp(cmp.matern, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh), options=list(verbose = TRUE))    

    ## 
    ## spde2
    ##
    oscillating.spde2 = inla.spde2.generic(M0 = M0, M1 = M1, M2 = M2, 
                                           B0 = B.phi0.oscillating, B1 = B.phi1.oscillating, B2 = B.phi2.oscillating, theta.mu = c(0, log(20), -5), 
                                           theta.Q = diag(c(10, 10, 10)), transform = "log")
    cmp.oscillating.spde2 <- coordinates ~ mySPDE(coordinates, model = oscillating.spde2, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
    fit.oscillating.spde2 <- lgcp(cmp.oscillating.spde2, data = Y.spdf, samplers = Ypos.sldf, domain = list(coordinates = mesh), options=list(verbose = TRUE))

    plot(fit.oscillating.spde2$marginals.hyperpar[[1]], type="l", xlim=c(-5,5))
    lines(fit.oscillating.spde2$marginals.hyperpar[[2]], lty=2)
    lines(fit.oscillating.spde2$marginals.hyperpar[[3]], lty=2)
}





## opt.theta.copyla <- cobyla(x0=opt.theta$par, fn = pthetapc.prop.marg.post_osc_temp, data=data, 
##                            X=Xest, Z=Zest, beta=betaest,
##                            mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
##                            A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                            W=W, acc=1e-3,
##                            print.verbose=FALSE,
##                            control = list(xtol_rel = 1e-4, maxeval = 2000))

## opt.theta.full.mle <- opt.theta

 ## opt.theta
## $par
## [1]  3.1495294 -0.7172584 -4.5455201  2.9199646 -0.1594968  1.7720353 -0.2435957

## $value
## [1] 8045.094

## $counts
## function gradient 
##      328       NA 

## $convergence
## [1] 0

## $message
## NULL

## opt.theta.no.prior 
## opt.theta.no.prior
## $par
## [1]  3.1336597 -0.6915424 -4.0642065  2.4643263  0.1568164  1.4272924 -0.7964659

## $value
## [1] 3690.094

## $counts
## function gradient 
##      570       NA 

## $convergence
## [1] 0

## $message
## NULL


## opt.theta.strong.prior
## $par
## [1]  3.136220401 -0.389107783 -3.818003347  1.805182245  0.009856284
## [6]  1.340381863 -0.824562382

## $value
## [1] 3707.501

## $counts
## function gradient 
##      714       NA 

## $convergence
## [1] 0

## $message
## NULL






## theta.hd <- seq(2.5,7.5,len=5)
## theta.temp <- seq(1,10,len=10)
## theta.var <- seq(0.1,4, len=4)
theta.rho <- seq(15,25,len=10)
val.prof <- NULL
for(i in 1:length(theta.rho)){    
    par.theta.profile <- c(log(22.9),
                           log(0.41),
                           -log((1-(-0.92))/(1+(-0.92))),
                           log(3.37),
                           log(3.19),
                           log(theta.rho[i]),
                           log(0.79))
    val.prof[i] <- pthetapc.prop.marg.post_osc_temp(par.theta.profile, data=data, 
                                                    X=Xinit, Z=Zinit, beta=betaest,
                                                    mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
                                                    A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
                                                    W=W, 
                                                    acc=1e-3)    
}


## theta.rho.hd <- seq(5, 8, len=5)
## theta.rho.temp <- seq(3, 6, len=5)
## par.theta.profile <- c(log(22.898),
##                        log(0.9002),
##                        -log((1-(-0.954))/(1+(-0.954))),
##                        log(6.203),
##                        log(0.757),
##                        log(3.852),
##                        log(0.438))

## mat.prof  <- matrix(0, nrow=nrow(theta.rho.hd), ncol=nrow(theta.rho.temp))
## for(i in 1:length(theta.rho)){
##     par.theta.profile <- c(log(21),
##                            log(1.2),
##                            -log((1-(-0.99))/(1+(-0.99))),
##                            log(5),
##                            log(0.9),
##                            log(3.5),
##                            log(0.5))
##     mat.prof[i,j] <- pthetapc.prop.marg.post_osc_temp(par.theta.profile, data=data, 
##                                                     X=Xinit, Z=Zinit, beta=betaest,
##                                                     mesh=mesh, mesh.theta=mesh.hd, mesh1d=mesh1d,
##                                                     A=A, Atilde=Atilde, Aobs=Aobs, Atildeobs=Atildeobs,
##                                                     W=W, 
##                                                     acc=1e-3)    
## }

## ## rho is:   sigma is:  phi is:   rho_hd is:   sigma_hd is:  rho_temporal is:   sigma_temporal is:  

## plot(theta.var, val.prof)

## plot(theta.var, val.prof)

## logrhol    <- log(10)
## logsigmaU  <- log(1) #                                                                                                                              
## alpha1     <- .001                                                                                                                                  
## alpha2     <- 1e-6                                                                                                                                  
## logrho.tempL    <- log(10)                                                                                                                          
## logsigma.tempU  <- log(1) #                                                                                                                         
## alphaT1     <- .001                                                                                                                                 
## alphaT2     <- 1e-6
## hyperpar   <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU),
##                    alphaT1=alphaT1, alphaT2=alphaT2, rho.tempL=exp(logrho.tempL), sigma.tempU=exp(logsigma.tempU))
## ## par=par.theta

## fn = pthetapc.prop.marg.post_osc_temp; hyperpar=hyperpar; data=data
## X=Xinit; Z=Zinit; beta=betaest; mesh=mesh; mesh.theta=mesh.hd; mesh1d=mesh1d; A=A; Atilde=Atilde; Aobs=Aobs; Atildeobs=Atildeobs;
## W=W;  acc=1e-3; print.verbose=TRUE; control=list(maxit=1000)

## save(opt.theta, file="fitted_model.RData")

## if(FALSE){
##         ## multiplicities for trapezoidal rule
##         intermediate.knots <- cbind(do.call("rbind",Ypos$sp)[filter.index,1][-1] -
##                                     do.call("rbind",Ypos$ep)[filter.index,1][-sum(filter.index)])
##         functions.multiplicity <- intermediate.knots

##         for(i in 1:length(functions.multiplicity)) {
##             if(i==1){
##                 if(functions.multiplicity[i]!=0) functions.multiplicity[i] <- 1
##                 else functions.multiplicity[i] <- 2
##             }else{
##                 if(functions.multiplicity[i]!=0) {
##                     functions.multiplicity[i] <- 1
##                     functions.multiplicity[i-1] <- 1} else{
##                                                         functions.multiplicity[i] <- 2
##                                                     }
##             }
##         }

##         functions.multiplicity <- c(1, functions.multiplicity, 1)
## }




