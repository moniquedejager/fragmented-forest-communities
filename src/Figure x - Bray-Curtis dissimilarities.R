# test how similar the subcommunities are; are subcommunities that are 
# spatially closer to one another also more similar? 

# hypothesis: Euclidian distance between subpopulations is negatively related 
# to similarity in community composition

filename <- 'results/dissimilarity/dissimilarity_data_1.txt'
df <- read.table(filename, header = TRUE)

df3 <- df[df$mu == 6,]
for (mu in c(1, 3, 5)){
  for (f_hab_loss in ((1:19)*0.05)){
    df2 <- df[(df$mu == mu)&(abs(df$f_hab_loss - f_hab_loss) < 0.01),]
    n <- length(df2$mu)
    df2 <- df2[(n-190+1):n,]
    df3 <- rbind(df3, df2)
  }
}

df <- df3

df$diss2 <- round(df$dissimilarity, 2)
group <- paste(df$mu, df$diss2,df$f_hab_loss,sep='-')
df2   <- data.frame(mu            = tapply(df$mu, group, mean),
                    f_hab_loss    = tapply(df$f_hab_loss, group, mean),
                    dissimilarity = tapply(df$diss2, group, mean),
                    z             = tapply(df$diss2, group, length))
df2$clustering <- 'Random'
df2$clustering[df2$mu == 3] <- 'Fractal'
df2$clustering[df2$mu == 5] <- 'Clustered'

ggplot(df2, aes(x=f_hab_loss, y=dissimilarity, fill=z/190)) + 
  geom_tile() + 
  facet_grid(cols=vars(clustering)) + 
  scale_fill_gradientn(colors = c('lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       labels = scales::percent_format(),
                       limits=c(1/190, 0.5),
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



