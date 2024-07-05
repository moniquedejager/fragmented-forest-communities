# test how similar the subcommunities are; are subcommunities that are 
# spatially closer to one another also more similar? 

# hypothesis: Euclidian distance between subpopulations is negatively related 
# to similarity in community composition

filename <- 'results/dissimilarity/dissimilarity_data_1.txt'
df <- read.table(filename, header = TRUE)

for (i in 1:10){
  filename <- paste('results/dissimilarity/dissimilarity_data_', i, '.txt', sep='')
  df2 <- read.table(filename, header = TRUE)
  #df <- df2[df2$dispersal == 'different',]
  #write.table(df, filename, append = FALSE, row.names = FALSE, col.names = TRUE)
  df <- rbind(df, df2)
}

df$diss2 <- round(df$dissimilarity, 2)
df$disp <- df$dispersal == 'different'
group <- paste(df$mu, df$diss2,df$f_hab_loss, df$disp,sep='-')
df2   <- data.frame(mu            = tapply(df$mu, group, mean),
                    f_hab_loss    = tapply(df$f_hab_loss, group, mean),
                    dissimilarity = tapply(df$diss2, group, mean),
                    z             = tapply(df$diss2, group, length),
                    disp          = tapply(df$disp, group, mean))
df2$clustering <- 'Random'
df2$clustering[df2$mu == 3] <- 'Fractal'
df2$clustering[df2$mu == 5] <- 'Clustered'

df2$disp2 <- 'Same dispersal'
df2$disp2[df2$disp == 1] <- 'Different dispersal'
df2$disp2 <- factor(df2$disp2, levels=unique(df2$disp2)[2:1])
ggplot(df2, aes(x=f_hab_loss*100, y=dissimilarity, fill=z/1900)) + 
  geom_raster() + 
  facet_grid(cols=vars(clustering), rows=vars(disp2)) + 
  scale_fill_gradientn(colors = c('white','lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       labels = scales::percent_format(),
                       limits=c(1/1900, 1),
                       na.value = 'white',
                       trans = 'log10',
                       name="% Subcommunity pairs") +
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 12)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())


