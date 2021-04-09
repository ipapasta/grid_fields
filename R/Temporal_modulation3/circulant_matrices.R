circulant <- function(x) {
    n <- length(x)
    mat <- matrix(NA, n, n)
    for (i in 1:n) {
        mat[i, ] <- c(x[-(1:(n + 1 - i))], x[1:(n + 1 - i)])
    }
    return(mat)
}


a <-  2
b <- -8
d <- -2
## 
x   <- c(a, b, 0, 0, d)
mat <- circulant(x)
## a^2 > 4*b*d  FALSE
n  <- length(x)
th <- acos(-a/(2*sqrt(b*d)))
th2 <- asin(sqrt(1-a^2/(4*b*d)))
th <- -atan(sqrt(4*b*d-a^2)/(-a))
t  <- sqrt(b/d)
j <- 0:(n-1)

## with t = 1
circulant(((sin(j * th) + sin((n-j)*th))/(b*sin(th)*(1-2*cos(n*th) + 1)))) - solve(mat)
## with t generic
circulant(((t^(j+1))*(sin(j * th) + (t^n)*sin((n-j)*th))/(b*sin(th)*(1-2*(t^n)*cos(n*th) + t^(2*n))))) - solve(mat)

(sin(j * th) + sin((n-j)*th))/(b*sin(th)*(-2*cos(n*th)))

circulant((sin(j * th) + sin((n-j)*th))/(b*sin(th)*(-2*cos(n*th))))

circulant(((t^(j+1))*(sin(j * th) + (t^n) * sin((n-j)*th))/(b*sin(th)*(1-2*cos(n*th) + t^(2*n))))) - solve(mat)
circulant(((sin(j * th) + sin((n-j)*th))/(-b*sin(th)*(2*cos(n*th)))))

solve(mat)


solve(mat)


acos(-(a)/(2*(b*d)^(1/2)))
sin((1-((a^2)/(4*b*d)))^(1/2))
