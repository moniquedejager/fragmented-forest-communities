
ix   <- seq(-5, 5, 1)
dx   <- rep(ix, length(ix))
dy   <- sort(dx)
dist <- sqrt(dx^2 + dy^2)
Pm   <- c(0.1, 0.2, 0.3, 0.4)

df <- data.frame(x = vector(length=0),
                 y = vector(length=0),
                 dist = vector(length=0),
                 Pm   = vector(length=0),
                 Pm2  = vector(length=0))

for (iPm in unique(Pm[Pm > 0])) {
  Pm2 <- (2*pi*iPm^2)^-1 * exp(-1*dist / iPm) 
  Pm2 <- Pm2 / sum(Pm2)
  
  df2 <- data.frame(x = dx,
                   y = dy,
                   dist = dist,
                   Pm   = iPm,
                   Pm2  = Pm2)
  df <- rbind(df, df2)
}

ggplot(df, aes(x=dist, y=Pm2, color=as.factor(Pm))) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(trans='log10') + 
  xlab('Distance from parental subcommunity') + 
  ylab('Dispersal probability') + 
  scale_color_discrete(name = expression(lambda)) + 
  theme_bw() + 
  theme(legend.position = 'top')

