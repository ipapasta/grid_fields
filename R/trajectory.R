ggplot() + geom_path(data=data.frame(X), aes(position_x, position_y)) + coord_fixed() +
    geom_point(data=Y, aes(position_x, position_y), colour=2)+theme_classic()
