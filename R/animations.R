N <- 100
coords.dir   <- expand.grid(seq(0, 100, len=N), seq(0, 100, len=100), seq(0,2*pi,len=50)) %>% unname
predict.data.space.direction.full  <- data.frame(hd=coords.dir[,3], coords.x1=coords.dir[,1], coords.x2=coords.dir[,2])
lambda.space.direction.full        <- predict(fit.space.direction, predict.data.space.direction.full, ~ Intercept + spde2)

df.animation.varying.direction <- lambda.space.direction.full %>% mutate(hd=as.factor(hd))
df.animation.varying.coordinates <- lambda.space.direction.full %>%
    mutate(coord.diag = (coords.x1==coords.x2)) %>% dplyr::filter(coord.diag==TRUE) %>%
    mutate(coords.x1  = as.factor(coords.x1))



## animations side by side
## animation of a path


p.space.direction.varying.direction  <- ggplot(df.animation.varying.direction %>% group_by(coords.x1, coords.x2) %>%
                                               summarize(mean=mean(mean)), aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    guides(fill = guide_colourbar(title="log firing rate"))+
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.direction.full$mean),max(lambda.space.direction.full$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11))


p.space.direction.varying.direction + geom_path(data = df.animation.varying.coordinates, aes(x=as.numeric(coords.x1), y=coords.x2))

transition_states(hd, transition_length = 1, state_length = 1)


## 
## code for animation
## 


p.space.direction.varying.direction  <- ggplot(df.animation.varying.direction, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    guides(fill = guide_colourbar(title="log firing rate"))+
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.direction.full$mean),max(lambda.space.direction.full$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11)) + transition_states(hd, transition_length = 1, state_length = 1)

p.space.direction.varying.coord <- ggplot(df.animation.varying.coordinates) + 
    geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    xlab("head direction") +
    ylab("log firing rate")+
    geom_line(aes(x=hd, y=mean)) +
    geom_hline(yintercept=0, colour="grey")+
    coord_polar(start = pi, direction=1) + theme_minimal()  +
    transition_states(coords.x1, transition_length = 1, state_length = 1)

## ggplot(df.animation, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
##     scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
##                          limits=c(min(lambda.space.direction.full$mean),max(lambda.space.direction.full$mean)))+
##     coord_fixed()+ 
##     theme_classic() + theme(legend.text=element_text(size=11)) + transition_states(hd, transition_length = 2, state_length = 1)

## p.space.direction  <- ggplot(df.animation, aes(coords.x1,coords.x2)) + geom_raster(aes(fill=mean), interpolate=TRUE) +
##     scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
##                          limits=c(min(lambda.space.direction.full$mean),max(lambda.space.direction.full$mean)))+
##     coord_fixed()+ 
##     theme_classic() + theme(legend.text=element_text(size=11)) + transition_states(hd, transition_length = 2, state_length = 1)
p.space.direction.varying.direction
p.space.direction.varying.coord

anim_save("anim_space_direction.varying.direction.gif", p.space.direction.varying.direction)
anim_save("anim_space_direction.varying.coord.gif", p.space.direction.varying.coord)



## 
## animations side by side for coordinates varying on the diagonal of Omega
##
p.space.direction.varying.direction  <- ggplot(df.animation.varying.direction, aes(coords.x1,coords.x2)) +
    geom_raster(aes(fill=mean), interpolate=TRUE) +
    guides(fill = guide_colourbar(title="log firing rate"))+
    scale_fill_gradientn(colours=ocean.balance(100), guide = "colourbar",
                         limits=c(min(lambda.space.direction.full$mean),max(lambda.space.direction.full$mean)))+
    coord_fixed()+ 
    theme_classic() + theme(legend.text=element_text(size=11)) + transition_states(hd, transition_length = 1, state_length = 1)

p.space.direction.varying.coord <- ggplot(df.animation.varying.coordinates) + 
    geom_ribbon(aes(x= hd, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70")+
    scale_x_continuous(breaks=seq(0,2*pi - pi/3,pi/3),
                       labels=paste0(0:5,expression("pi"),"/",3)) +
    xlab("head direction") +
    ylab("log firing rate")+
    ## scale_y_continuous(breaks=seq(min(lambda.space.direction.average.coord$q0.025) %>% floor,
    ##                               max(lambda.space.direction.average.coord$q0.975) %>% ceiling, by=0.5),
    ##                    limits=c(min(lambda.space.direction.average.coord$q0.025) %>% floor,
    ##                             max(lambda.space.direction.average.coord$q0.975) %>% ceiling))+
    geom_line(aes(x=hd, y=mean)) +
    geom_hline(yintercept=0, colour="grey")+
    coord_polar(start = pi, direction=1) + theme_minimal()  +
    transition_states(coords.x1, transition_length = 1, state_length = 1)



## animation of temporal intensity (mouse walks on the path)
predict.intensity.on.path  <- data.frame(coords.x1=Ypos$coords[,1], coords.x2=Ypos$coords[,2], hd=Ypos$hd, firing_times=Ypos$time)
## lambda.M0.path
## lambda.M1.path

window_width <- nrow(lambda.M0.path)/1500  # How much of the whole data to show at once
frames       <- 5000   # Increase to make smoother animation & bigger file
shift_per_frame = (nrow(lambda.M0.path) - window_width) / frames

# This bit of purrr copies the whole data frame [frames] times, identifying each with "id"
lambda.M0.path_copied <- map_df(seq_len(frames), ~lambda.M0.path, .id = "id") %>%
    mutate(id = as.integer(id)) %>%
    dplyr::filter(time >= id * shift_per_frame,
                  time <= id * shift_per_frame + window_width)

lambda.M1.path_copied <- map_df(seq_len(frames), ~lambda.M1.path, .id = "id") %>%
    mutate(id = as.integer(id)) %>%
    dplyr::filter(time >= id * shift_per_frame,
                  time <= id * shift_per_frame + window_width)

p.M0.path.time.anim  <- ggplot(lambda.M0.path_copied) +
    geom_ribbon(aes(x = time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") +
    geom_line(aes(x=time,y=mean)) +
    ## geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025, length(firing_times)))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11)) +
    transition_manual(id) +
    view_follow(fixed_y=TRUE)

p.M1.path.time.anim  <- ggplot(lambda.M1.path_copied) +
    geom_ribbon(aes(x = time, ymin=q0.025, ymax=q0.975), alpha=0.4, colour="grey70") +
    geom_line(aes(x=time,y=mean)) +
    ## geom_point(data=Y, aes(x=firing_times, rep(min(lambda.M0.path$q0.025, length(firing_times)))), colour="red") + 
    theme_classic() + theme(legend.text=element_text(size=11)) +
    transition_manual(id) +
    view_follow(fixed_y=TRUE)

  ##   ggplot(lambda.M0.path_copied, aes(x = firing_times, y = )) + 
  ## geom_line() +


animate(p.M0.path.time.anim, nframes = frames)
anim_save("anim_temporal_intensity.M0.gif", p.M0.path.time.anim)
