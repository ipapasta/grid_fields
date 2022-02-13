## posterior distribution of gridness score
## input is a fitted model with inlabru (e.g. fit.space and fit.space.direction)
posterior.gridness.score <- function(inlabru.fitted.object){
    marg         <- inla.tmarginal(theta.2.phi, inlabru.fitted.object$marginals.hyperpar[["Theta3 for spde2"]])
    summaries    <- inla.zmarginal(marg, silent=TRUE)
    hpd.interval <- inla.hpdmarginal(0.95, marg)
    attr(marg, "summary") <- list(interval.estimate.hpd = hpd.interval, point.estimates = summaries)
    return(marg)
}

## posterior.summaries.gridness.score <- function(inlabru.fitted.object){
##     marg         <- posterior.gridness.score(inlabru.fitted.object = inlabru.fitted.object)

## }

posterior.gridness.score.M0 <- posterior.gridness.score(fit.space)
posterior.gridness.score.M1 <- posterior.gridness.score(fit.space.direction)


par(mfrow=c(1,2))
plot(posterior.gridness.score.M0, type="l", xlim=c(-1,1), ylim=c(0, 10))
lines(seq(-.99, 1, len=100), prior.phi_osc(seq(-.99, 1, len=100), a=2, b=10, lg=FALSE), lty=2)
## abline(v=attr(posterior.gridness.score.M0, "summary")$interval.estimate.hpd)
plot(posterior.gridness.score.M1, type="l", xlim=c(-1,1), ylim=c(0, 10))
lines(seq(-.99, 1, len=100), prior.phi_osc(seq(-.99, 1, len=100), a=2, b=10, lg=FALSE), lty=2)
## abline(v=attr(posterior.gridness.score.M1, "summary")$interval.estimate.hpd)
