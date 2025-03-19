# test how similar the subcommunities are; are subcommunities that are 
# spatially closer to one another also more similar? 

# hypothesis: Euclidian distance between subpopulations is negatively related 
# to similarity in community composition

filename <- 'results/simulation results different dispersal 20240611/dissimilarity/dissimilarity_data_1.txt'
df <- read.table(filename, header = TRUE)
df$dispersal_type <- 2
df <- df[(df$f_hab_loss > 0)&(df$mutation_rate == 0.0003),]

filename <- 'results/simulation results similar dispersal 20240612/dissimilarity/dissimilarity_data_1.txt'
df2 <- read.table(filename, header = TRUE)
df2$dispersal_type <- 1
df  <- rbind(df, df2)

# plot using geom_raster:
df$diss2 <- round(df$dissimilarity, 2)
group <- paste(df$mu, df$f_hab_loss, df$diss2, df$dispersal_type, sep='-')
df2   <- data.frame(mu            = tapply(df$mu, group, mean),
                    f_hab_loss    = tapply(df$f_hab_loss, group, mean),
                    dissimilarity = tapply(df$diss2, group, mean),
                    z             = tapply(df$diss2, group, length),
                    disp_type     = tapply(df$dispersal_type, group, mean))

df3 <- df2
df3$clustering <- 'Random'
df3$clustering[df3$mu == 3] <- 'Fractal'
df3$clustering[df3$mu == 5] <- 'Clustered'

df3$dispersal_type <- 'Same dispersal'
df3$dispersal_type[df3$disp_type == 2] <- 'Different dispersal'
df3$dispersal_type <- factor(df3$dispersal_type, 
                             levels=c('Same dispersal', 
                                      'Different dispersal'))

#tapply(df$dissimilarity, paste(df$mu, df$f_hab_loss, df$dispersal_type), length)
ggplot(df3, aes(x=f_hab_loss, y=dissimilarity, fill=z/1900)) + 
  geom_raster() + 
  facet_grid(cols=vars(clustering),rows=vars(dispersal_type)) + 
  scale_fill_gradientn(colors = c('lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       labels = scales::percent_format(),
                       na.value = 'white',
                       trans = 'log10',
                       name="% Subcommunity pairs",
                       limits=c(0.001, 1)) +
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 12)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())


