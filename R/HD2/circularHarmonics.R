C.circular.harmonics <- function(order, inverse=FALSE){
    if(!inverse) return(Diagonal(n = 1+2*order, x=c(2*pi, rep(pi, 2*order))))
    if(inverse) return(Diagonal(n = 1+2*order, x=c(1/(2*pi), rep(1/pi, 2*order))))
}

G.circular.harmonics <- function(order){
    Diagonal(n = 1+2*order, x=c(0, rep(pi, 2*order))*c(1, sort(rep((1:order)^2, 2))))
}

circular.harmonics <- function(hd, order){
    as.numeric(c(1, sapply((1:order)*hd, function(x) c(cos(x), sin(x)))))
}

circular.make.A <- function(hd.vec, order){
    hd.mat <- matrix(hd.vec, ncol=1)
    t(apply(hd.mat, 1, circular.harmonics, order=order))
}


## C         <- C.circular.harmonics(order, inverse=FALSE)
## C.inverse <- C.circular.harmonics(order, inverse=TRUE)
## G         <- G.circular.harmonics(order)
## (G %*% C.inverse %*% G)/pi: diagonal with squares of squares
