library(GeneralizedHyperbolic)
library(rmutil)
delta <- 1
alpha <- 1
beta  <- alpha/3
rho   <- 0.7
n     <- 10000
eps   <- rnig(n, delta=1, alpha=alpha, beta=beta)
X1   <- rnig(n, delta=delta/(1-rho), alpha=alpha, beta=beta)
X2   <- rho*X1 + eps

Y1 <- qlaplace(pnig(X1, delta=delta/(1-rho), alpha=alpha, beta=beta))
Y2 <- qlaplace(pnig(X2, delta=delta/(1-rho), alpha=alpha, beta=beta))


par(mfrow=c(1,2))
plot(X1, X2, main="NIG")
abline(a=0, b=1)
abline(a=0, b=rho)
abline(a=0, b=1/rho)
plot(Y1, Y2, main="Laplace")
abline(a=0, b=1)
abline(a=0, b=rho)
x <- seq(0, 100, len=1000)
lines(x, rho*x, col=3, lty=2)


## ------------------
## tail approximation
## ------------------
delta <- 1
alpha <- 1
beta  <- alpha/3

gamma <- sqrt(alpha^2-beta^2)
x <- seq(-20, 20, len=100)
uterm <- alpha*sqrt(delta^2 + x^2)
dens <- (alpha*delta/pi)*(alpha*besselK(uterm, 1)/(uterm))*exp(delta*gamma + beta*x)
plot(x, dnig(x, mu=0,delta=1,alpha=alpha,beta=beta), type="l")
lines(x, dens, col=2)

plot(x, exp(delta*gamma)*(2*alpha*delta/pi)*(2*beta*uterm - alpha*x*(besselK(uterm,0)+besselK(uterm,2))/besselK(uterm,1))^(-1), ylim=c(0, 0.2))
lines(x, 1-pnig(x, delta=delta, alpha=alpha, beta=beta), col=2)

