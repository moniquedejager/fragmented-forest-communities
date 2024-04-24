# test how similar the subcommunities are; are subcommunities that are 
# spatially closer to one another also more similar? 

# hypothesis: Euclidian distance between subpopulations is negatively related 
# to similarity in community composition

filename <- 'results/simulation results 04-2024/dissimilarity/dissimilarity_data_1.txt'
df <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('results/simulation results 04-2024/dissimilarity/dissimilarity_data_', i, '.txt', sep='')
  df2 <- read.table(filename, header = TRUE)
  df  <- rbind(df, df2)
}

df$mu2 <- paste("μ = ", df$mu, sep='')
df$f2  <- paste(df$f_hab_loss, "% habitat loss", sep='')
df$f2  <- factor(df$f2, levels=unique(df$f2))

df$dist2 <- round(df$distance + 5, -1)
df$dist3 <- paste(df$dist2 - 10, '-', df$dist2, ' cells apart', sep='')

# plot using geom_raster:
df$diss2 <- round(df$dissimilarity, 2)
group <- paste(df$dist2, df$mu, df$f_hab_loss, df$diss2, sep='-')
df2   <- data.frame(mu            = tapply(df$mu, group, mean),
                    f_hab_loss    = tapply(df$f_hab_loss, group, mean),
                    dissimilarity = tapply(df$diss2, group, mean),
                    dist          = tapply(df$dist2, group, mean),
                    z             = tapply(df$diss2, group, length))

df2$mu2   <- paste("μ = ", df2$mu, sep='')
df2$dist3 <- paste(df2$dist - 10, '-', df2$dist, ' cells apart', sep='')
ggplot(df2, aes(x=f_hab_loss, y=dissimilarity, fill=z)) + 
  geom_raster() + 
  facet_grid(rows=vars(mu2),cols=vars(dist3)) + 
  scale_fill_gradientn(colors = c('yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       na.value = 'white',
                       trans = 'log10',
                       name="# Occasions") +
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  theme_classic() + 
  guides(fill = guide_colorbar(barwidth = 30)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

sel <- df2$dist %in% c(10, 20, 30)
ggplot(df2[sel,], aes(x=f_hab_loss, y=dissimilarity, fill=z)) + 
  geom_raster() + 
  facet_grid(cols=vars(mu2),rows=vars(dist3), switch = 'y') + 
  scale_fill_gradientn(colors = c('lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       na.value = 'white',
                       trans = 'log10',
                       name="# Occurrences") +
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 20)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

df3 <- df2[sel,]
df3 <- df3[df3$mu != 2,]
df3$clustering <- 'Random'
df3$clustering[df3$mu == 3] <- 'Fractal'
df3$clustering[df3$mu == 5] <- 'Clustered'

ggplot(df3, aes(x=f_hab_loss, y=dissimilarity, fill=z)) + 
  geom_raster() + 
  facet_grid(cols=vars(clustering),rows=vars(dist3), switch = 'y') + 
  scale_fill_gradientn(colors = c('lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       na.value = 'white',
                       trans = 'log10',
                       name="# Occurrences") +
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 12)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

# is there a difference in dissimilarity between distances? 
sel <- df2$f_hab_loss %in% c(0, 10, 30, 50, 70, 90)
ggplot(df2[sel,], aes(x=dist, y=dissimilarity, fill=z)) + 
  geom_raster() + 
  facet_grid(cols=vars(mu),rows=vars(f_hab_loss), switch = 'y') + 
  scale_fill_gradientn(colors = c('lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       na.value = 'white',
                       trans = 'log10',
                       name="# Occurrences") +
  xlab('Distance between subcommunities') + 
  ylab('Bray-Curtis dissimilarity') + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 12)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
