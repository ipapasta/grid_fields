fake.locations = matrix(c(0,0,10,10, 0, 10, 10, 0), nrow = 4, byrow = T)
mesh.sim = inla.mesh.2d(loc = fake.locations, max.edge=c(0.4, 1))

spde = inla.spde2.pcmatern(mesh.sim, prior.range = c(.5, .5), prior.sigma = c(.5, .5))

inla.seed = sample.int(n=1E6, size=1)
sigma.u = 1.5
# - the marginal standard deviation of the spatial field
range = 2

Qu = inla.spde.precision(spde, theta=c(log(range), log(sigma.u)))
u = inla.qsample(n=1, Q=Qu, seed = inla.seed)
u = u[,1]

local.plot.field(u, mesh.sim)
len = range
# - the true range
arrows(5-0.5*len, 6, 5+0.5*len, 6, length=0.05, angle=90, code=3, lwd=3)


## -----------------------

RFgetModelNames(type="positive definite", domain="single variable",
                iso="isotropic") 

## our choice is the exponential model;
## the model includes nugget effect and the mean:
model <- RMgauss(var=5, scale=1) + # with variance 4 and scale 10
    RMnugget(var=1) + # nugget
    RMtrend(mean=.0) # and mean

## define the locations:
from <- 0
to <- 20
N <- 30
x.seq <- seq(from, to, length=N) 
y.seq <- seq(from, to, length=N)

simu1 <- RFsimulate(model, x=x.seq, y=y.seq)
df1 <- data.frame(cbind(expand.grid(1:N, 1:N), simu1$variable1))
names(df1) <- c("x", "y", "log_lambda")
simu2 <- RFsimulate(model, x=x.seq, y=y.seq)
df2 <- data.frame(cbind(expand.grid(1:N, 1:N), simu2$variable1))
names(df2) <- c("x", "y", "log_lambda")

p1 <- ggplot(data=df) + geom_raster(aes(x=x, y=y, fill=log_lambda), interpolate=TRUE) +
    scale_fill_gradient2(low = "darkblue", mid = "green", high = "red", midpoint = 0) +
    coord_fixed()+ ## geom_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
    ## aes(x=position_x, y=position_y),color="red") +
    theme_classic() + theme(axis.text=element_text(size=20),
                            axis.title=element_text(size=20),
                            text = element_text(size=24)) +
    xlab("x") + ylab("y")
p2 <- ggplot(data=df2) + geom_raster(aes(x=x, y=y, fill=log_lambda), interpolate=TRUE) +
    scale_fill_gradient2(low = "darkblue", mid = "green", high = "red", midpoint = 0) +
    coord_fixed()+ ## geonm_point(data=data$Y, size=(rank(lambdafit)/(1+n)),
    ## aes(x=position_x, y=position_y),color="red") + 
    theme_classic() + theme(axis.text=element_text(size=20),
                            axis.title=element_text(size=20),
                            text = element_text(size=24)) +
    xlab("x") + ylab("y")


a <- grid.arrange(p1, p2, nrow=1)

ggsave(filename="matern_field.ps", a, width = 55, height = 25, units = "cm")

x11()
image(simu)

image(gridfield, mat, asp=1)
