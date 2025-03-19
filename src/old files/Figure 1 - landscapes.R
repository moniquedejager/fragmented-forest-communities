# show the patterning of the landscapes
library(ggplot2)

df <- data.frame(mu   = vector(length = 0),
                 loss = vector(length = 0),
                 x    = vector(length = 0),
                 y    = vector(length = 0),
                 hab  = vector(length = 0))

x <- rep(1:20, 20)
y <- sort(x)

for (mu in c(1, 3, 5)){
  filename <- paste('results/simulation results different dispersal 20240611/landscapes/landscapes_4_', mu, '.txt', sep='')
  landscape<- as.matrix(read.table(filename))  
  
  for (i in 1:20){
    loss <- seq(0, 95, 5)[i]
    df2  <- data.frame(mu   = mu,
                       loss = loss,
                       x    = x,
                       y    = y,
                       hab  = landscape[,i])
    df   <- rbind(df, df2)
  }
}

df$loss2 <- paste(df$loss, '% habitat loss', sep='')
df$mu2             <- "μ = 1"
df$mu2[df$mu == 3] <- "μ = 3"
df$mu2[df$mu == 5] <- "μ = 5"


sel <- (df$loss %in% c(10, 30, 50, 70, 90)) & (df$hab == 1)
ggplot(df[sel,], aes(x=x, y=y)) +
  geom_raster(fill='green4') +
  facet_grid(cols=vars(loss2), rows=vars(mu2), switch = 'y') + 
  theme_bw() + 
  xlab('') + 
  ylab('') + 
  theme(strip.placement = "outside", 
        strip.background = element_blank())

df2 <- df
df2$clustering <- 'Random'
df2$clustering[df2$mu == 3] <- 'Fractal'
df2$clustering[df2$mu == 5] <- 'Clustered'

sel <- (df2$loss %in% c(10, 30, 50, 70, 90)) & (df2$hab == 1)
p1 <- ggplot(df2[sel,], aes(x=x, y=y)) +
  geom_raster(fill='green4') +
  facet_grid(rows=vars(loss2), cols=vars(clustering), switch = 'y') + 
  theme_bw() + 
  xlab('') + 
  ylab('') + 
  theme(strip.placement = "outside", 
        strip.background = element_blank())

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 1.tif', 
     width = 4, height = 6, units = 'in', res = 600)
p1
dev.off()

