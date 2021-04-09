## Simulation from head directional model

source("hd_precision.R")
## variance of oscillating process - tests

N <- 1000
mesh.test     <- inla.mesh.1d(seq(0, 2*pi, len=100), boundary="cyclic", degree=1)

theta1 <- log(8)
theta2 <- log(1)
rho     <- exp(theta1)
kappa   <- sqrt(8*(3/2))/rho
sigma   <- exp(theta2)
tausq   <- 1/(4*(sigma^2)*(kappa^3))
Q1     <- hd.bsp.precision(c(theta1, theta2), mesh=mesh.test)
x1 <- matrix(nrow=N, ncol=mesh.test$n)
for(i in 1:N){
    x1[i,]  <- inla.qsample(n=1, Q=Q1)
}

plot(apply(x1, 1, var), type="l")
matplot(mesh.test$loc, t(x1), type="l", lty=1, col=1)

coth <- function(x) (exp(2*x) + 1)/(exp(2*x) - 1)
csch <- function(x) (2*exp(2*x))/(exp(2*x) - 1)

## variance
(1/((tausq)*sqrt(2*pi)))*(pi*coth(pi*kappa)+(pi^2)*kappa*csch((pi*kappa^2)))/kappa^3




theta1 <- c(log(15.14), log(2.22))
theta2 <- c(2, 2)
theta3 <- c(3, 2)
theta4 <- c(.04, 2)
theta.nodes <- seq(0, 2*pi, len=25)
mesh.hd <- inla.mesh.1d(theta.nodes, boundary="cyclic", degree=1)
Q1   <- hd.bsp.precision(theta1, mesh=mesh.hd);  Q2   <- hd.bsp.precision(theta2, mesh=mesh.hd)
Q3   <- hd.bsp.precision(theta3, mesh=mesh.hd);  Q4   <- hd.bsp.precision(theta4, mesh=mesh.hd)
x1  <- inla.qsample(n=1, Q=Q1);  x2  <- inla.qsample(n=1, Q=Q2)
x3  <- inla.qsample(n=1, Q=Q3);  x4  <- inla.qsample(n=1, Q=Q4)
## 
coords.theta    <- seq(0, 2*pi, length=100)
Ahd.pred  <- inla.mesh.projector(mesh.hd, loc=coords.theta)$proj$A
field <- cbind(as.numeric(Ahd.pred %*% x1), as.numeric(Ahd.pred %*% x2),  as.numeric(Ahd.pred %*% x3), as.numeric(Ahd.pred %*% x4))
df.plot    <- data.frame(time=coords.theta, field1=field[,1], field2=field[,2], field3=field[,3], field4=field[,4])
lims  <- c(min(c(df.plot$field1,df.plot$field2,df.plot$field3,df.plot$field4)), max(c(df.plot$field1,df.plot$field2,df.plot$field3,df.plot$field4)))
p1    <- ggplot(data = df.plot) + geom_line(aes(time,field1))  + ylim(lims) +  theme_classic()
p2    <- ggplot(data = df.plot) + geom_line(aes(time,field2))  + ylim(lims) +  theme_classic()
p3    <- ggplot(data = df.plot) + geom_line(aes(time,field3))  + ylim(lims) +  theme_classic()
p4    <- ggplot(data = df.plot) + geom_line(aes(time,field4))  + ylim(lims) +  theme_classic()
a <- grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
plot(a)

p1    <- ggplot(data = df.plot) + geom_line(aes(time,field1))  + theme_classic()

ggsave(filename="hd_fields.svg", plot=a )
dev.off()



## Oscillating model

source("osc_precision.R")
library(ggplot2)
library(gridExtra)
library(pals)
theta1 <- c(3.2, -1.51, -5.5)
## print(paste("rho is:", exp(par.theta[1]), " sigma is:", exp(par.theta[2]),
##             "phi is: ", (1-exp(-par.theta[3]))/(1+exp(-par.theta[3])),
##             "rho_hd is:", exp(par.theta[4]), " sigma_hd is:", exp(par.theta[5]),
##             "rho_temporal is: ", exp(par.theta[6]), "sigma_temporal is: ", exp(par.theta[7])))   
## rho is: 24.6315067497884  sigma is: 0.223412648495406 phi is:  -0.984156846971105 rho_hd is: 2.34785300331212  sigma_hd i\
## s: 1.85599826799832 rho_temporal is:  1.09712442633471 sigma_temporal is:  0.537161073174284"
Q1   <- osc.precision(theta1, mesh)
x1  <- inla.qsample(n=1, Q=Q1)
## 
coords  <- expand.grid(seq(-20, 120, seq=100), seq(-20, 120, seq=100))
proj.s.pred  <- inla.mesh.projector(mesh, loc=as.matrix(coords))
Aosc.pred      <- proj.s.pred$proj$A
field <- cbind(as.numeric(Aosc.pred %*% x1))
df.plot       <- data.frame(x=coords[,1], y=coords[,2], field1=field[,1])
lims <- c(min(c(df.plot$field1, df.plot$field2, df.plot$field3, df.plot$field4)),max(c(df.plot$field1, df.plot$field2, df.plot$field3, df.plot$field4)))
p1  <- ggplot(data = df.plot) + geom_raster(aes(x,y,fill=exp(field1)), interpolate=TRUE) + 
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar") + coord_fixed() + theme_classic()
plot(p1)



## kronecker covariance
theta.theta    <- c(0.84, 1)
theta.omega <- c(3.2, -1, -9)
Q.theta  <- hd.bsp.precision(theta.theta, mesh=mesh.hd)
Q.omega  <- osc.precision(theta.omega, mesh)
Q        <- kronecker(Q.theta, Q.omega)

sum(diag(expand(Cholesky(Q))$L)>0) 

