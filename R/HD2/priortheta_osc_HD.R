prior.phi_osc <- function(phi, a, b, l=(-0.99), u=1, lg=TRUE){
    if(lg) {
        o  <- -log(u-l)+dbeta((phi-l)/(u-l), shape1=a, shape2=b, log=TRUE)
        return(o)
    }
    if(!lg) {
        o  <- (1/(u-l))*dbeta((phi-l)/(u-l), shape1=a, shape2=b)
        return(o)
    }
}

prior.theta_osc_HD <- function(rho, sigma, phi, rhoL, sigmaU, alpha1, alpha2, # hyperparameters of spatial component
                               rho.HD, sigma.HD, rho.HDL, sigma.HDU, alphaHD1, alphaHD2, 
                               lg=TRUE)
{
    ## input:
    ##       range: range parameter of spatial field
    ##       sigma: sigma parameter of spatial field
    ##       rhoL: lower quantile of prior distribution of rho, P(rho < rhoL) = alpha1
    ##       sigmaU: upper quantile of prior distribution of sigma, P(sigma < sigmaL) = alpha2
    ## penalised complexity for Gaussian Markov random field where the
    ## base model has infinite range and zero variance
    sigmaLN <- 0.3
    murho   <- 20
    ## sigmaTLN <- 0.2
    ## murhoT   <- 10
    l1tilde   <- -log(alpha1)*rhoL
    l2tilde   <- -log(alpha2)/sigmaU
    l1HDtilde <- -log(alphaHD1)*rho.HDL
    l2HDtilde <- -log(alphaHD2)/sigma.HDU
    lpphi         <- prior.phi_osc(phi, a=2, b=20, lg=TRUE)
    ## lpsigma    <- dlnorm(sigma, meanlog=log(1), sdlog = .5, log=TRUE) #check
    ## lpsigma       <- dexp(sigma, rate=0.01, log=TRUE) #check
    lpsigmaLN     <- dlnorm(sigma, meanlog=log(.5), sdlog=.3, log=TRUE)
    ## lprho.HD   <- dlnorm(rho.HD, meanlog = log(.4), sdlog=1., log=TRUE)
    const         <- -(1/(sqrt(2*pi)))*log(0.9)
    ## lprho.HD      <- -log(2) - (3/2)*log(rho.HD) + log(1/const) - (1/const)*sqrt(rho.HD)
    lprho.HD      <- dlnorm(rho.HD, meanlog=log(pi/6), sdlog=.5, log=TRUE)
    ## lpsigma.HD <- dlnorm(sigma.HD, meanlog = log(.5), sdlog=.5, log=TRUE)
    ## lpsigma.HD <- dexp(sigma.HD, rate=1, log=TRUE) #check
    lpsigma.HD    <- dlnorm(sigma.HD, meanlog=log(.5), sdlog=.3, log=TRUE)
    ## lphdrho    <- prior.phi_osc(rho.HD, a=10, b=2, l=0, u=2*pi, lg=TRUE)
    ## outl       <- log(l1tilde) + log(l2tilde) - l1tilde/rho - l2tilde*sigma ## + lpphi 
    lprhoLN       <- dlnorm(rho, meanlog=log(22), sdlog=.005, log=TRUE)
    const2        <- -30 * log(1-1e-8)
    ## lprho      <- log(const2)-2*log(rho) - const2/rho
    ## lprho        <- dexp(rho, rate=2, log=TRUE)
    outl        <- lprhoLN + lpsigmaLN + lpphi + lpsigma.HD + lprho.HD 
    ## out      <- l1tilde*l2tilde*(rho^(-2))*exp(-l1tilde/rho-l2tilde*sigma)
    ## out      <- l2tilde*exp(-l2tilde*sigma) * dlnorm(rho, log(murho), sigmaLN, log=FALSE)
    if(lg) return(outl) else return(out)
    ## if(lg) return(0) else return(out)
}

## 
## construct HD prior based on number of stationary points.
##

## THE FOLLOWING FUNCTION IS USED BELOW IN VISUALIZATIONS
prior.theta_osc <- function(rho, sigma, phi, rhoL, sigmaU, alpha1, alpha2, lg=TRUE)
{
    ## input:
    ##       range: range parameter of spatial field
    ##       sigma: sigma parameter of spatial field
    ##       rhoL: lower quantile of prior distribution of rho, P(rho < rhoL) = alpha1
    ##       sigmaU: upper quantile of prior distribution of sigma, P(sigma < sigmaL) = alpha2
    ## penalised complexity for Gaussian Markov random field where the
    ## base model has infinite range and constant variance
    sigmaLN <- 0.1
    murho   <- 20
    l1tilde <- -log(alpha1)*rhoL
    l2tilde <- -log(alpha2)/sigmaU
    lpphi   <- prior.phi_osc(phi, a=2, b=20, lg=TRUE)
    ## outl    <- log(l1tilde) + log(l2tilde) - l1tilde/rho - l2tilde*sigma ## + lpphi 
    lrhoLN  <- dlnorm(rho, log(murho), sigmaLN, log=TRUE)
    ## outl    <- lrhoLN + log(l2tilde) - l2tilde*sigma + lpphi
    outl    <- lrhoLN + dlnorm(sigma, meanlog=log(2), sdlog = .4, log=TRUE) + lpphi
    ## out     <- l1tilde*l2tilde*(rho^(-2))*exp(-l1tilde/rho-l2tilde*sigma)
    ## out     <- l2tilde*exp(-l2tilde*sigma) * dlnorm(rho, log(murho), sigmaLN, log=FALSE)* dlnorm(sigma, log(2), sdlog=.4, log=FALSE)
    ## out     <- exp(outl)
    if(lg) return(outl) else return(out)
}

## 
## Experimentation
##

if(FALSE)
{
    
    ## plot((1-exp(-v))/(1+exp(-v)))   
    ## grid of rho and sigma values to evaluate prior
    rhoo   <- seq(1, 80, len= 200)       
    sigmaa <- seq(1e-6, 10, len= 200)       
    rhoolg   <- seq(1e-1, 1000, len= 200)       
    sigmaalg <- seq(0.001, 5, len= 200)       
    mat    <- matrix(nrow=200, ncol=200)
    matlg    <- matrix(nrow=200, ncol=200)
    logrhoL    <- log(10)
    logsigmaU  <- log(1) #
    alpha1     <- 5## .001
    alpha2     <- 0.4
    ## logrhoL    <- log(10)
    ## logsigmaU  <- log(.25)
    ## alpha1     <- .001
    ## alpha2     <- 1e-80
    hyperpar   <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU))
    ## 
    for(i in 1:length(rhoo))
        for(j in 1:length(sigmaa))
            mat[i,j] <- prior.theta_osc(rho=rhoo[i], sigma=sigmaa[j], phi=-.99,
                                    rhoL=hyperpar$rhoL, sigmaU=hyperpar$sigmaU,
                                    alpha1=hyperpar$alpha1, alpha2=hyperpar$alpha2, lg=FALSE)
    ## 
    for(i in 1:length(rhoo))
        for(j in 1:length(sigmaa))
            matlg[i,j] <- prior.theta_osc(rho=rhoolg[i], sigma=sigmaalg[j], phi=0,
                                        rhoL=hyperpar$rhoL, sigmaU=hyperpar$sigmaU,
                                        alpha1=hyperpar$alpha1, alpha2=hyperpar$alpha2, lg=TRUE)    
    ##
    palette   <- colorRampPalette(gray(seq(1, .1, len=100)))
    palettelg <- colorRampPalette(gray(seq(1, 0.01, len=100)))
    ##
    ## postscript("prior.ps", heigh=600, width=900)
    par(mfrow=c(1,2))
    filled.contour(rhoo, sigmaa, mat, xlab="rho", ylab="st deviation", main="joint prior",
                   ylim=c(0, 10), nlevels=50, color.palette = palette,
                   plot.axes = { contour(rhoo, sigmaa, mat,  nlevels = 50, lwd=.5,
                                         drawlabels = TRUE, axes = FALSE, 
                                         frame.plot = FALSE, add = TRUE);
                                         axis(1); axis(2) } )
    
    ##    
    filled.contour(log(rhoolg), log(sigmaalg), matlg, xlab="rho", ylab="st deviation", main="joint (log) prior",
                   nlevels=30, color.palette = palettelg, 
                   plot.axes = { contour(rhoo, sigmaa, matlg,  nlevels = 30, lwd=.5,
                                         drawlabels = TRUE, axes = FALSE, 
                                         frame.plot = FALSE, add = TRUE);
                                         axis(1); axis(2) } )
    
    ## ## log = TRUE
    ## filled.contour(rhoo, sigmaa, mat, xlab="range", ylab="st deviation", main="penalised complexity prior",
    ##                ylim=c(0, 5),
    ##                plot.axes = { contour(rhoo, sigmaa, mat,  nlevels = 20, lwd=.5,
    ##                                      drawlabels = TRUE, axes = FALSE, 
    ##                                      frame.plot = FALSE, add = TRUE);
                                         ## axis(1); axis(2) } )    
}








##############################################################################
## reference: 1) https://www.tandfonline.com/doi/pdf/10.1080/01621459.2017.1415907?needAccess=true
## 2) TBP: computing with oscillating covariance functions
##############################################################################

## -----------------------------------
## some guidelines for prior selection
## -----------------------------------
## https://groups.google.com/forum/#!topic/r-inla-discussion-group/dunoXK_yAco

## The parameters "prior.range" and "prior.sd" control the joint prior
## on range and standard deviation of the spatial field.

## "prior.range" contains two values: range0 and Prange, such that the
## resulting prior satisfies P( range < range0) = Prange. The logic is
## that prior is constructed to "shrink" towards high values of range
## and that the user must select which ranges are so short that they
## are unfeasible for the problem at hand. For example, with range0 =
## 0.01 and Prange = 0.05 there is an a priori probability of only 5%
## that the range will be less than 0.01. From a theoretical point of
## view, it is challenging to come up with clear guidelines for how
## range0 and Prange should be selected, but, from a practical point
## of view, it is advisable not to favour ranges that are smaller than
## the resolution of the mesh. In the case that the range is similar
## to the edge lengths or smaller there are large discretisation
## errors in the SPDE model. However, the specification does not need
## to be based on the mesh size. If you have reason to believe that it
## is unreasonable that the range is less than, say 0.05, you would
## base the prior on that and not on a, potentially, much smaller
## triangle size in the mesh.

## It is possible to force stronger smoothing (longer ranges) by
## making the prior on range stricter. For example, (range0 = 0.05,
## Prange = 0.05) will result in higher ranges than (range0 = 0.01,
## Prange = 0.05).

## "prior.sigma" contains two values: sigma0 and Psigma, such that the
## resulting prior satisfies P(sigma > sigma_0) = Psigma. The logic is
## that the prior is constructed to shrink towards small values of the
## standard deviation for the spatial field and that the user must
## select which standard deviations for the spatial field that are so
## high that they are unfeasible for the problem at hand. This is
## dependent on the scaling of the data and the specific application
## so I do not know any clear guidelines to follow. If you have other
## similar data available, it could be used to come up with an idea
## about how large the spatial effect can be and then encode this into
## the prior by making standard deviation corresponding to
## unreasonably large spatial effects unlikely.


## let theta_3 = \sinc(\pi \tilde{\theta}) \in (0, 1).
##
## For optimization purposes, the following parameterisation is used
## phi = {1-exp(-theta_3)}/{1+exp(-theta_3)} \in (-\infty, infty) : 
##
## if \phi = \cos(\pi \tilde{\theta}) then the following relation can be used
##
## tau^2   =  {4* \pi* sigma^2* kappa* sinc(pi\tilde{\theta})}^{-1}
##
## with sinc( \pi \tilde{\theta}) = \sqrt{(1-\phi^2)}/{\arccos \phi}
##
## 

## 
## "prior.phi" contains two values: phi0 and Pphi, such that the
## resulting prior satisfies P(phi > sigma_0) = Psigma. The logic is
## that the prior is constructed to shrink towards small values of the
## standard deviation for the spatial field and that the user must
## select which standard deviations for the spatial field that are so
## high that they are unfeasible for the problem at hand. This is
## dependent on the scaling of the data and the specific application
## so I do not know any clear guidelines to follow. If you have other
## similar data available, it could be used to come up with an idea
## about how large the spatial effect can be and then encode this into
## the prior by making standard deviation corresponding to
## unreasonably large spatial effects unlikely.
