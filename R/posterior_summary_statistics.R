source("Functions.R")
## 
## 
## 
posterior.gridness.score.M0         <- posterior.spatial.gridness.score(fit.space, theta.mapping=function(x) theta.2.phi(x, l=l, u=u))
posterior.gridness.score.M1         <- posterior.spatial.gridness.score(fit.space.direction, theta.mapping=function(x) theta.2.phi(x, l=l, u=u))
posterior.spatial.stdev.M0          <- posterior.spatial.standard.deviation(fit.space)
posterior.spatial.stdev.M1          <- posterior.spatial.standard.deviation(fit.space.direction)
posterior.spatial.range.M0          <- posterior.spatial.range(fit.space)
posterior.spatial.range.M1          <- posterior.spatial.range(fit.space.direction)
## 
posterior.range.M0     <- posterior.spatial.range(fit.space)
posterior.range.M1     <- posterior.spatial.range(fit.space.direction)
## 
posterior.directional.stdev.M1          <- posterior.directional.standard.deviation(fit.space.direction)
posterior.directional.range.M1          <- posterior.directional.range(fit.space.direction)
posterior.directional.kappa.M1          <- posterior.directional.kappa(fit.space.direction)
## 
## 
##


par(mfrow=c(1,3))
plot(posterior.gridness.score.M0, type="l", xlim=c(-1,1), xlab="phi", ylab="density", main="M0",
     ylim=c(min(posterior.gridness.score.M0[,2]), max(posterior.gridness.score.M0[,2])))
## abline(v=c(theta.2.phi(initial.space$theta3, l=l, u=u),l))
lines(seq(-.99, 1, len=1000), prior.phi_osc(seq(-.99, 1, len=1000), a=a.par.phi.prior.spatial.oscillating, b=b.par.phi.prior.spatial.oscillating, l=l,u=u,lg=FALSE), lty=2)
## 
plot(posterior.spatial.stdev.M0, type="l", xlim=c(0,10), xlab="sigma.spatial", ylab="density", main="M0",
     ylim=c(min(posterior.spatial.stdev.M0[,2]), max(posterior.spatial.stdev.M0[,2])))
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), 1/2), lty=2)
## 
plot(posterior.spatial.range.M0, type="l", xlim=c(0,60), xlab="rho.spatial", ylab="density", main="M0",
     ylim=c(min(posterior.spatial.range.M0[,2]), max(posterior.spatial.range.M0[,2])))
lines(seq(5, 60, len=1000), dlnorm(seq(5, 60, len=1000)-5, log(mu.range.spatial.oscillating), sigma.range.spatial.oscillating, log=FALSE), lty=2)



pdf("posterior_hyperparameters2.pdf")

par(mfrow=c(4,2))
plot(posterior.gridness.score.M0, type="l", xlim=c(-1,1), ylim=c(0, 10), xlab="phi", ylab="density", main="M0")
## abline(v=c(theta.2.phi(initial.space$theta3, l=l, u=u),l))
lines(seq(-.99, 1, len=1000), prior.phi_osc(seq(-.99, 1, len=1000), a=a.par.phi.prior.spatial.oscillating, b=b.par.phi.prior.spatial.oscillating, l=l,u=u,lg=FALSE), lty=2)
abline(v=l)
plot(posterior.gridness.score.M1, type="l", xlim=c(-1,1), ylim=c(0, 10), xlab="phi", ylab="density", main="M1")
lines(seq(-.99, 1, len=1000), prior.phi_osc(seq(-.99, 1, len=1000), a=a.par.phi.prior.spatial.oscillating, b=b.par.phi.prior.spatial.oscillating, l=l,u=u,lg=FALSE), lty=2)
## 
plot(posterior.spatial.stdev.M0, type="l", xlim=c(0,10), ylim=c(0, 2), xlab="sigma.spatial", ylab="density", main="M0")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), 1/2), lty=2)
plot(posterior.spatial.stdev.M1, type="l", xlim=c(0,10), ylim=c(0, 2), xlab="sigma.spatial", ylab="density", main="M1")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), 1/2), lty=2)
## 
plot(posterior.spatial.range.M0, type="l", xlim=c(0,60), ylim=c(0, .61), xlab="rho.spatial", ylab="density", main="M0")
lines(seq(5, 60, len=1000), dlnorm(seq(5, 60, len=1000)-5, log(mu.range.spatial.oscillating), sigma.range.spatial.oscillating, log=FALSE), lty=2)
## lines(seq(5, 100, len=1000), , lty=2)
plot(posterior.spatial.range.M1, type="l", xlim=c(0,60), ylim=c(0, .61), xlab="rho.spatial", ylab="density", main="M1")
lines(seq(0, 100, len=1000), dlnorm(seq(0, 100, len=1000), log(mu.range.spatial.oscillating), sigma.range.spatial.oscillating), lty=2)
## par(mfrow=c(1,3))
plot(posterior.directional.range.M1, type="l", xlim=c(0,10), ylim=c(0, 4), xlab="rho.directional", ylab="density", main="M1")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), rho.directional), lty=2)
plot(posterior.directional.stdev.M1, type="l", xlim=c(0,10), ylim=c(0, 2), xlab="sigma.directional", ylab="density", main="M1")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), sigma.directional), lty=2)

dev.off()


## 
plot(posterior.directional.kappa.M1, type="l", xlim=c(0,10), ylim=c(0, 6), xlab="kappa.directional", ylab="density", main="M1")
lines(seq(0, 10, len=1000), sqrt(8*(3/2))*dexp(sqrt(8*(3/2))/seq(0, 10, len=1000), rho.directional)/(seq(0, 10, len=1000))^2, lty=2)


## lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), rho.directional), lty=2)

## 
## par(mfrow=c(1,2))

