## --------------------------------------------
dx.dt <- diff(X$position_x)/diff(X$synced_time)
dy.dt <- diff(X$position_y)/diff(X$synced_time)

speed <- sqrt(dx.dt^2 + dy.dt^2)

tmp <- cbind(speed, X$speed[-1])        # row match
plot(tmp)
abline(a=0,b=1)

plot(dx.dt, dy.dt)


## Y matrix

dx.dt <- diff(Y$position_x)/diff(Y$firing_times)
dy.dt <- diff(Y$position_y)/diff(Y$firing_times)

speed <- sqrt(dx.dt^2 + dy.dt^2)

tmp <- cbind(speed, X$speed[-1])        # row match
abline(a=0,b=1)

plot(dx.dt, dy.dt)




Ypos$speed
Ypos <- Ypos %>% mutate(speed.lead = lead(speed), dt=c(diff(synced_time)[1], diff(synced_time) ))  %>% head(nrow(Ypos)-1) %>% 
    mutate(val=dt*((speed + speed.lead)/2))


sum(X$val)
