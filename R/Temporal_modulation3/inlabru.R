options("rgdal_show_exportToProj4_warnings"="none")
## -----------------
## GORILLAS examples
## -----------------
## reprex({
library(INLA)
library(inlabru)
data("gorillas")
mesh   <- gorillas$mesh
nests  <- gorillas$nests
boundary <- gorillas$boundary
fem.mesh.spatial <- inla.mesh.fem(mesh, order = 2)
B.phi0.spatial = matrix(c(0,1,0,0), nrow=1)
B.phi1.spatial = matrix(c(0,0,1,0), nrow=1)
B.phi2.spatial = matrix(c(0,0,0,1), nrow=1)
M0 = fem.mesh.spatial$c0
M1 = fem.mesh.spatial$g1
M2 = fem.mesh.spatial$g2
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


oscillating.rgeneric <- inla.rgeneric.define(oscillating.model, M = list(M0=M0, M1=M1, M2=M2)) #arguments that need to be passed are M0, M1 and M2
cmp.rgeneric <- coordinates ~ mySmooth(coordinates, model = oscillating.rgeneric,
                                       mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept

fit <- lgcp(cmp.rgeneric, nests, samplers = boundary, domain = list(coordinates = mesh), options=list(verbose = TRUE))


lambda <- predict(fit,
                  pixels(gorillas$mesh, mask = gorillas$boundary),
                  ~ exp(mySmooth + Intercept))


library(pals)
ggplot() + 
    gg(lambda) +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar")+
    gg(gorillas$nests, color = "red", size = 0.2) +
    coord_equal() +
    ggtitle("Nest intensity per km squared")

    ## traceback()
## }



## ,style=TRUE,session_info=TRUE, advertise = FALSE, html_preview    = TRUE,
## comment         = "#;-)", tidyverse_quiet = FALSE, std_out_err     = TRUE
## )



## -----------------
## MRSEA examples
## -----------------

str(mrsea$points) # this contains the data
dim(mrsea$points@data)
head(mrsea$points@data)
str(mrsea$samplers) #this contains line transects
dim(mrsea$samplers@data)
head(mrsea$samplers@data)

mrsea$points@coords %>% dim
mrsea$samplers@lines

ips <- ipoints(mrsea$samplers, mrsea$mesh, group = "season")


ggplot() +
  ## gg(mrsea$mesh) +
  ## gg(mrsea$boundary) +
  ## gg(mrsea$samplers) +
  gg(mrsea$points, size = 0.5) +
  coord_fixed() +
  facet_wrap(~season) +
  ggtitle("MRSea observation seasons")





oscillating.spde2 = inla.spde2.generic(M0 = M0.spatial, M1 = M1.spatial, M2 = M2.spatial, 
                                 B0 = B.phi0.spatial, B1 = B.phi1.spatial, B2 = B.phi2.spatial, theta.mu = c(0, 0, 0), 
                                 theta.Q = diag(c(10, 10, 10)), transform = "log")


B.matern.phi0.spatial = matrix(c(0,1,0), nrow=1)
B.matern.phi1.spatial = matrix(c(0,0,1), nrow=1)
matern.spde2 = inla.spde2.generic(M0 = M0.spatial, M1 = M1.spatial, M2 = M2.spatial, 
                                 B0 = B.matern.phi0.spatial, B1 = B.matern.phi1.spatial, B2 = 1, theta.mu = c(0, 0), 
                                 theta.Q = diag(c(10, 10)))

cmp.matern <- coordinates ~ mySmooth(coordinates, model = matern.spde2, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
fit.matern <- lgcp(cmp.matern, nests, samplers = boundary, domain = list(coordinates = mesh), options=list(verbose = TRUE))








cmp.matern <- coordinates ~ mySmooth(coordinates, model = matern.spde2, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
fit.matern <- lgcp(cmp.matern, nests, samplers = boundary, domain = list(coordinates = mesh), options=list(verbose = TRUE))



oscillating.rgeneric <- inla.rgeneric.define(oscillating.model, M = list(M0=M0.spatial, M1=M1.spatial, M2=M2.spatial)) 
cmp.rgeneric <- coordinates ~ mySmooth(coordinates, model = f(coordinates, model = oscillating.rgeneric),
                                       mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept
fit <- lgcp(cmp.rgeneric, nests, samplers = boundary, domain = list(coordinates = mesh), options=list(verbose = TRUE))



colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE)
  )
}
## ---------------------------
## construction of spde models
## ---------------------------

## Oscillating precision
## -------------------------
d      <- 2
alpha  <- 2
nu     <- alpha-d/2

par.theta <- c(0,0,-1)

theta1 <- par.theta[1]
theta2 <- par.theta[2]
theta3 <- par.theta[3]
rho     <- exp(theta1)
sigma   <- exp(theta2)
kappa   <- sqrt(8)/rho
phi     <- (1-exp(-theta3))/(1+exp(-theta3))
sincpth <- sqrt(1-phi^2)/acos(phi)
tausq   <- 1/(4*pi*(sigma^2)*(kappa^2)*sincpth)
## ---------------------------------------
mat.gen <- tausq*((kappa^4)*(fem$c0) + (2*(kappa^2)*phi*(fem$g1)) + fem$g2)


fem       <- inla.mesh.fem(mesh, order = 2)
##
param <- list(B.tau=tausq, B.kappa=kappa^2, B.phi = phi, theta.prior.prec=0.1, theta.prior.mean=0)
M0 = fem$c0
M1 = fem$g1
M2 = fem$g2


## example oscillating
B.phi0.spatial = matrix(c(0,1,0,0), nrow=1)
B.phi1.spatial = matrix(c(0,0,1,0), nrow=1)
B.phi2.spatial = matrix(c(0,0,0,1), nrow=1)
M0.spatial = fem$c0
M1.spatial = fem$g1
M2.spatial = fem$g2

## head directional
theta.nodes.irregular <- c(0,sort(runif(18, 0, 2*pi)), 2*pi)
theta.nodes <- seq(0, 2*pi, len=1000)
h <- diff(theta.nodes)[1]
mesh.hd     <- inla.mesh.1d(theta.nodes, boundary="cyclic", degree=1)

fem.mesh.hd <- inla.mesh.fem(mesh.hd, order = 2)
##
B.phi0.hd = matrix(c(0,1,0), nrow=1)
B.phi1.hd = matrix(c(0,0,1), nrow=1)
M0.hd = fem.mesh.hd$c0
M1.hd = fem.mesh.hd$g1
M2.hd = fem.mesh.hd$g2

M0.hd[1,1] == h
M0.hd[1,1] - h
M1.hd[1,1] == 2/h
M1.hd[1,1] - 2/h
M1.hd[1,2] == -1/h
M2.hd[1,1] == (6/h^3) # numerical?
M2.hd[1,1] - (6/h^3)
M2.hd[1,2] == -(4/h^3)
M2.hd[1,19] == -(4/h^3) # numerical?
M2.hd[1,19] + (4/h^3)   #wrong calculations here, this has been fixed
oscillating = inla.spde2.generic(M0 = M0.spatial, M1 = M1.spatial, M2 = M2.spatial, 
                                 B0 = B.phi0.spatial, B1 = B.phi1.spatial, B2 = B.phi2.spatial, theta.mu = c(3.5,-3.5,-10), 
                                 theta.Q = diag(c(0,0,0)), transform = "log")

directional = inla.spde2.generic(M0 = M0.hd, M1 = M1.hd, M2 = M2.hd, 
                                 B0 = B.phi0.hd, B1 = B.phi1.hd, B2 = 1, theta.mu = c(0,0), 
                                 theta.Q = c(.1,.1))

theta1.spatial <- 3.5
theta2.spatial <- -3.5
theta3.spatial <- 3
true_Q.spatial <- inla.spde.precision(spde=oscillating, theta = c(theta1.spatial, theta2.spatial, theta3.spatial))

theta1.directional <- 1
theta2.directional <- 1.10
kappa.dir <- sqrt(exp(theta2.directional))
tau.dir   <- exp(theta1.directional)

true_Q.directional <- (tau.dir)^2*((kappa.dir^4)*(fem.mesh.hd$c0) + (2*(kappa.dir^2)*(fem.mesh.hd$g1)) + fem.mesh.hd$g2)
true_S.directional <- solve(true_Q.directional)
true_field.directional <- inla.qsample(1000, true_Q.directional)

theta <- seq(0, 2*pi, len=100)
r1 <- NULL
for(i in 1:length(theta)){
    A0 <- inla.mesh.projector(mesh.hd, loc=0)$proj$A
    Atheta <- inla.mesh.projector(mesh.hd, loc=theta[i])$proj$A
    r1[i] <- as.numeric(A0 %*% (true_S.directional %*% t(Atheta)))
}

true_cov <- function(theta, kappa, tausq){
    (1/(tausq * sinh(pi*kappa)^2*(2*kappa)^2)) *
        ((theta/2)*cosh((2*pi-theta)*kappa) +
         (pi - (theta/2))*cosh(theta*kappa) + (cosh((pi-theta)*kappa)*sinh(pi*kappa))/kappa)
}

pdf("/Users/ipapasta/Desktop/alpha2_covariance_check.pdf", height=7, width=11)
plot(theta,r1, type="l", ylab="Covariance", xlab="angle")
points(theta, true_cov(theta, kappa.dir, tausq = exp(1)^2), col=2)
points(0, (pi + sinh(2 * pi * kappa.dir)/(2*kappa.dir))/((sinh(pi*kappa.dir)^2)*(exp(1)^2)*(2*kappa.dir)^2), pch=16, col=1)
dev.off()

true_Q.directional <- inla.spde.precision(spde=directional, theta = c(theta1.directional, theta2.directional))
true_S.directional <- solve(true_Q.directional)

## Sigma_11 : 6.877131e-04
a <-  (kappa.dir^4 * h) + (4 * kappa.dir^2/h) + 6/h^3
b <- -(2*kappa.dir^2/h)-4/h^3
d <- b
## 
x   <- c(a, b, rep(0, mesh.hd$n-3), d)
mat <- (tau.dir^2)*circulant(x)
## a^2 > 4*b*d
n  <- length(x)
th <- -atan(sqrt(4*b*d-a^2)/(-a))
t  <- sqrt(b/d)
j  <- 0:(n-1)


true_Q.directional[1,1]
(tau.dir^2)*(h*kappa.dir^4 + 4*(kappa.dir^2)/h + 6/h^3)

true_Q.directional[1,2]
(tau.dir^2)*(-2*((kappa.dir^2)/h) - 4/h^3)


## with t = 1
oo <- ((1/tau.dir)^2) * circulant(((t^(j+1))*(sin(j * th) + (t^n)*sin((n-j)*th))/(b*sin(th)*(1-2*(t^n)*cos(n*th) + t^(2*n))))) - true_S.directional

S.searle <- ((1/tau.dir)^2) * circulant(((t^(j+1))*(sin(j * th) + (t^n)*sin((n-j)*th))/(b*sin(th)*(1-2*(t^n)*cos(n*th) + t^(2*n)))))






## generate a sample from the model

true_field.spatial <- inla.qsample(1, true_Q.spatial)[, 1]
true_field.directional <- inla.qsample(1, true_Q.directional)[, 1]

truth.spatial <- expand.grid(
  lon = seq(0, 100, length = 100),
  lat = seq(0, 100, length = 100)
)

truth.directional <- expand.grid(
    hd = seq(0, 2*pi, length = 100)
)

truth.spatial$field <- inla.mesh.project(mesh,
                                         loc = as.matrix(truth.spatial),
                                         field = true_field.spatial
                                         )

truth.directional$field <- inla.mesh.project(mesh.hd,
                                         loc = as.matrix(truth.directional),
                                         field = true_field.directional
                                         )


coordinates(truth.spatial) <- c("lon", "lat")
truth.spatial <- as(truth.spatial, "SpatialPixelsDataFrame")
pl_truth.spatial <- ggplot() +
  gg(truth.spatial, mapping = aes_string("lon", "lat", fill = "field")) +
    coord_equal() +
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(truth.spatial$field),max(truth.spatial$field)))+
    ggtitle("True field") + theme_classic()
## 
pl_truth.directional <- ggplot(data=truth.directional) +
    geom_line(aes(x=hd, y=field))+
    ## coord_equal() +
    ggtitle("True field") + theme_classic()

pl_truth.spatial
pl_truth.directional

csc <- colsc(truth$field)

multiplot(pl_truth, pl_truth + csc, cols = 2)







## GORILLAS
cmp <- coordinates ~ mySmooth(coordinates, model = oscillating, mapper=bru_mapper(mesh, indexed=TRUE)) + Intercept(1)

## Fit LGCP model
fit <- lgcp(cmp,
            data = SpatialPoints(cbind(data$Y$position_x, data$Y$position_y)),
            ## samplers = gorillas$boundary,
            domain = list(coordinates = mesh),
            options = list(control.inla = list(int.strategy = "eb")))



## VARIANCE of directional field
library(pracma)
theta1.directional <- 0
theta2.directional <- -0.10
true_Q.directional <- inla.spde.precision(spde=directional, theta = c(theta1.directional, theta2.directional))
rho     <- exp(theta2.directional)
kappa   <- sqrt(8*(3/2))/rho
tausq   <- exp(theta1.directional)
## tausq   <- 1/sigma

(1/(2*tausq*(kappa^4)*pi)) + (1/(tausq*pi))*(-2 + pi * kappa * (coth(pi * kappa) + pi * kappa * csch(pi * kappa)^2))/(4 * kappa^4)

    



((pi* coth(pi*kappa) + (pi^2)*kappa*csch((pi*kappa)))^2/kappa^3)/(tausq * sqrt(2*pi))
true_field.directional <- inla.qsample(1, true_Q.directional)[, 1]
directional.MC <- inla.qsample(1000, true_Q.directional)
var(directional.MC[10,])








## example matern
## B.phi0 = param$B.tau
## B.phi1 = 2 * param$B.kappa
## B.phi2 = param$B.phi




## ------------------------------------
## inlabru
## ------------------------------------

data(gorillas, package = "inlabru")
nests <- gorillas$nests
mesh <- gorillas$mesh
boundary <- gorillas$boundary


matern <- inla.spde2.pcmatern(mesh,
  prior.sigma = c(0.1, 0.01),
  prior.range = c(5, 0.01)
)




traceback()



##
data(mrsea, package = "inlabru")
ips <- ipoints(mrsea$samplers, mrsea$mesh, group = "season")



kappa <- 1
phi <- .1
h <- 0.1

solve(circulant(c((kappa^2)*h + (2*phi/h), -phi/h,0,0,0,0,0,-phi/h)))
