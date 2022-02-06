res.x = excursions(alpha=.95, u=0, mu=Apred%*%Xest,
                   Q=Apred %*% (Hessianest[(2:(length(Xest)+1)), (2:(length(Xest)+1))]%*% t(Apred)), 
                   type='!=', verbose=1, max.threads=0)


res.x = excursions(alpha=.95, u=0, mu=Xest[2:(mesh$n+1)],
                   Q=Hessianest[(2:(mesh$n+1)), (2:(mesh$n+1))], 
                   type='!=', verbose=1, max.threads=0)


map <- contourmap(mu=Xest[2:(mesh$n+1)],
                  Q=Hessianest[(2:(mesh$n+1)), (2:(mesh$n+1))], n.levels=10,
                  max.threads=0)

plot(map$map)

map<-continuous(res.x, mesh)


reo = mesh$idx$loc
cols = contourmap.colors(map, col=heat.colors(100, 1),
                         credible.col = grey(0.5, 1))
names(cols) = as.character(-1:2)

res.z = excursions(alpha=.95, u=0, mu=Zest,
                   Q=Hessianest[(length(Xest)+2):ncol(Hessianest), (length(Xest)+2):ncol(Hessianest)], 
                        type='!=', verbose=1, max.threads=0)


df2      <- data.frame(x=mesh$loc[,1], y=mesh$loc[,2],
                       rho=res.x$rho, F=res.x$rho)


p1  <- ggplot()  +
    geom_point(data=df2, size=4.5, colour=gray(1-df2$scaledX), aes(x,y)) +
    geom_point(data=data$Y, aes(x=position_x, y=position_y),color="red") + coord_fixed() +
    ylim(0,100) + xlim(0,100) + theme_classic()

plot(Zest, type="l")
plot(res.z$F)
lines(res.z$rho, col=2)
