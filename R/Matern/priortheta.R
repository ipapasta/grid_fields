##############################################################################
## reference: https://www.tandfonline.com/doi/pdf/10.1080/01621459.2017.1415907?needAccess=true
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

## At this point I think one has to do some experimenting with the
## priors. As far as I know, we still do not have enough experience
## with the priors to come up with clearer guidelines.


prior.theta <- function(rho, sigma, rhoL, sigmaU, alpha1, alpha2, lg=TRUE)
{
    ## input:
    ##       range: range parameter of spatial field
    ##       sigma: sigma parameter of spatial field
    ##       rhoL: lower quantile of prior distribution of rho, P(rho < rhoL) = alpha1
    ##       sigmaU: upper quantile of prior distribution of sigma, P(sigma < sigmaL) = alpha2
    ## penalised complexity for Gaussian Markov random field where the
    ## base model has infinite range and constant variance
    l1tilde <- -log(alpha1)*rhoL
    l2tilde <- -log(alpha2)/sigmaU
    outl    <- log(l1tilde) + log(l2tilde) - l1tilde/rho - l2tilde*sigma
    out     <- l1tilde*l2tilde*(rho^(-2))*exp(-l1tilde/rho-l2tilde*sigma)
    if(lg) return(outl) else return(out)
}


## 
## Experimentation
##

if(FALSE)
{
    
        ## grid of rho and sigma values to evaluate prior
    rhoo   <- seq(0.01, 50, len= 100)       
    sigmaa <- seq(0.01, 20, len= 100)       
    mat    <- matrix(nrow=100, ncol=100)

    ## 3
    logrhoL    <- log(20)
    logsigmaU  <- log(5)
    alpha1     <- .4
    alpha2     <- .01
    hyperpar   <- list(alpha1=alpha1, alpha2=alpha2, rhoL=exp(logrhoL), sigmaU=exp(logsigmaU))        

    for(i in 1:length(rhoo))
        for(j in 1:length(sigmaa))
            mat[i,j] <- prior.theta(rho=rhoo[i], sigma=sigmaa[j],
                                    rhoL=hyperpar$rhoL, sigmaU=hyperpar$sigmaU,
                                    alpha1=hyperpar$alpha1, alpha2=hyperpar$alpha2, lg=FALSE)


    palette <- colorRampPalette(gray(seq(1, .1, len=100)))
    ##
    postscript("prior.ps", heigh=600, width=900)
    filled.contour(rhoo, sigmaa, mat, xlab="range", ylab="st deviation", main="penalised complexity prior",
                   ylim=c(0, 5), nlevels=20, color.palette = palette,
                   plot.axes = { contour(rhoo, sigmaa, mat,  nlevels = 20, lwd=.5,
                                         drawlabels = TRUE, axes = FALSE, 
                                         frame.plot = FALSE, add = TRUE);
                                         axis(1); axis(2) } )
    dev.off()


    filled.contour(rhoo, sigmaa, mat, xlab="range", ylab="variance", main="penalised complexity prior",
                   nlevels=30, ylim=c(0,10), color.palette = palette,
                   plot.axes = { contour(rhoo, sigmaa, mat,  nlevels = 30, lwd=.5,
                                         drawlabels = TRUE, axes = FALSE, 
                                         frame.plot = FALSE, add = TRUE);
                                         axis(1); axis(2) } )
    

    
    prior.theta(10, .1,
                rhoL=hyperpar$rhoL,
                sigmaU=hyperpar$sigmaU,
                alpha1=hyperpar$alpha1,
                alpha2=hyperpar$alpha2, lg=TRUE)

    prior.theta(10, .01,
                rhoL=hyperpar$rhoL,
                sigmaU=hyperpar$sigmaU,
                alpha1=hyperpar$alpha1,
                alpha2=hyperpar$alpha2, lg=TRUE)
    
    }






