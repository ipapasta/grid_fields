if(FALSE){
    vec.temp <- data.frame(X=rank(Xbeta[1]+as.numeric(A[1,]%*%Xbeta[-1]))/(1+length(Xbeta[-1])))

    df.temp <- data.frame(do.call("rbind", Ypos$s.midpoints), X=rank(Xbeta[1]+as.numeric(A%*%Xbeta[-1]))/(1+length(Xbeta[-1])))
    df.temp <- data.frame(do.call("rbind", Ypos$s.midpoints), X=Xbeta[1]+as.numeric(A%*%Xbeta[-1]))
    names(df.temp) <- c("x", "y", "X")
    ggplot(data=df.temp[20:100,])  + geom_raster(aes(x=x,y=y,fill=X)) + coord_fixed()+theme_classic() 


    plot(df.temp[1:1000,1:3], add=TRUE)

    par(mfrow=c(1,2))
    plot(df.temp[1:100,1:2], add=TRUE)
    plot(df.temp[1:100,1:2], add=TRUE)
    for(i in 1:500){
        plot(df.temp[,1:3]) 
    }
    k <- 1.
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
    ##         k <- 1.5
    ## mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03,120), cutoff=k/2)
    ## 
    k <- 1.
    mesh <- inla.mesh.2d(mycoords, max.edge=c(k, 25*k), offset=c(0.03, 120), cutoff=k/2)
    ## 
    par(mfrow=c(1,2))
    plot(mesh, asp=1)
    plot(mycoords, add=TRUE, col=2, pch=16)
    ## 
    plot(trajectory)
    points(mycoords, pch=16, cex=0.7, col=2)
    ## plot(trajectory, add=TRUE)
    save(mesh, file="mesh.RData")

    df.W <- cbind(data.frame(line.id=1:nrow(s.psi)),
                  data.frame(origin=line.segments$split.origin),
                  data.frame(Li=do.call("c", Ypos$Li)),
                  s.tv,
                  data.frame(knot1=1:length(do.call("c", Ypos$Li)),
                             knot2=2:(1+length(do.call("c", Ypos$Li)))),
                  s.psi,
                  data.frame(psi.tilde.1 = rep(1/2, length(do.call("c", Ypos$Li)))),
                  data.frame(psi.tilde.2 = rep(1/2, length(do.call("c", Ypos$Li))))) %>%
        mutate(jk = pmap(list(knot1, knot2, v1, v2, v3), function(x, y, z, w, k) expand.grid(c(x, y), c(z, w, k)))) %>%
        mutate(x  = pmap(list(psi.tilde.1, psi.tilde.2, psi.1, psi.2, psi.3, Li), function(x, y, z, w, k, L){
            matrix(apply(expand.grid(c(x, y), c(z, w, k)), 1, prod)*L, ncol=1) 
        })) %>%
        mutate(jkx = map2(jk, x, function(x1, x2) {
            df <- cbind(x1, x2)
            names(df) <- c("knot", "v", "value")
            df
        })) %>%
        as_tibble

    zero.entries <- data.frame(expand.grid(c(1,2), setdiff(1:p, unique(do.call("rbind", df.W$jkx)$v))),
                               rep(0, dim(expand.grid(c(1,2), setdiff(1:p, unique(do.call("rbind", df.W$jkx)$v))))[1]))
    names(zero.entries) <- c("knot", "v", "value")
    W.input.to.sparse.matrix <- rbind(do.call("rbind", df.W$jkx), zero.entries) %>% arrange(v)
    W  <- sparseMatrix(i=W.input.to.sparse.matrix$knot, j=W.input.to.sparse.matrix$v, x=W.input.to.sparse.matrix$value)
    object_size(W)
    dim(W)
}
