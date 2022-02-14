source("Functions.R")

## 
## 
## 
posterior.gridness.score.M0 <- posterior.spatial.gridness.score(fit.space)
posterior.gridness.score.M1 <- posterior.spatial.gridness.score(fit.space.direction)
posterior.spatial.stdev.M0          <- posterior.spatial.standard.deviation(fit.space)
posterior.spatial.stdev.M1          <- posterior.spatial.standard.deviation(fit.space.direction)
posterior.spatial.range.M0          <- posterior.spatial.range(fit.space)
posterior.spatial.range.M1          <- posterior.spatial.range(fit.space.direction)
## 
posterior.range.M0          <- posterior.spatial.range(fit.space)
posterior.range.M1          <- posterior.spatial.range(fit.space.direction)
## 
posterior.directional.stdev.M1          <- posterior.directional.standard.deviation(fit.space.direction)
posterior.directional.range.M1          <- posterior.directional.range(fit.space.direction)
posterior.directional.kappa.M1          <- posterior.directional.kappa(fit.space.direction)


## 
## 
## 
par(mfrow=c(3,2))
plot(posterior.gridness.score.M0, type="l", xlim=c(-1,1), ylim=c(0, 10), xlab="phi", ylab="density", main="M0")
lines(seq(-.99, 1, len=1000), prior.phi_osc(seq(-.99, 1, len=1000), a=2, b=10, lg=FALSE), lty=2)
plot(posterior.gridness.score.M1, type="l", xlim=c(-1,1), ylim=c(0, 10), xlab="phi", ylab="density", main="M1")
lines(seq(-.99, 1, len=1000), prior.phi_osc(seq(-.99, 1, len=1000), a=2, b=10, lg=FALSE), lty=2)
## 
## par(mfrow=c(1,2))
plot(posterior.spatial.stdev.M0, type="l", xlim=c(0,10), ylim=c(0, 2), xlab="sigma.spatial", ylab="density", main="M0")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), 1/2), lty=2)
plot(posterior.spatial.stdev.M1, type="l", xlim=c(0,10), ylim=c(0, 2), xlab="sigma.spatial", ylab="density", main="M1")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), 1/2), lty=2)
## 
## par(mfrow=c(1,2))
plot(posterior.spatial.range.M0, type="l", xlim=c(0,60), ylim=c(0, .61), xlab="rho.spatial", ylab="density", main="M0")
lines(seq(0, 100, len=1000), dlnorm(seq(0, 100, len=1000), log(25), 3), lty=2)
plot(posterior.spatial.range.M1, type="l", xlim=c(0,60), ylim=c(0, .61), xlab="rho.spatial", ylab="density", main="M1")
lines(seq(0, 100, len=1000), dlnorm(seq(0, 100, len=1000), log(25), 3), lty=2)


par(mfrow=c(1,3))
plot(posterior.directional.range.M1, type="l", xlim=c(0,10), ylim=c(0, 4), xlab="rho.directional", ylab="density", main="M1")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), rho.directional), lty=2)
plot(posterior.directional.stdev.M1, type="l", xlim=c(0,10), ylim=c(0, 2), xlab="sigma.directional", ylab="density", main="M1")
lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), sigma.directional), lty=2)
## TODO: add prior of kappa
plot(posterior.directional.kappa.M1, type="l", xlim=c(0,10), ylim=c(0, 6), xlab="kappa.directional", ylab="density", main="M1")
lines(seq(0, 10, len=1000), sqrt(8*(3/2))*dexp(sqrt(8*(3/2))/seq(0, 10, len=1000), rho.directional)/(seq(0, 10, len=1000))^2, lty=2)


## lines(seq(0, 10, len=1000), dexp(seq(0, 10, len=1000), rho.directional), lty=2)

## 
## par(mfrow=c(1,2))

