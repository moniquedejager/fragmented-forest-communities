# How does the dispersal capacity change with increasing habitat fragmentation?
library(ggplot2)
source('./src/summarySE.R')

filename <- './results/dispersal_capacity/fragmented_dispersal capacity_data_0.05.txt'
df       <- read.table(filename, header = TRUE)
ix       <- 1:length(df$mu)
group    <- paste(df$mu, df$disp_cap, sep='-')
ix       <- tapply(ix, group, max)
df2      <- df[ix,]


for (i in seq(0.1, 0.95, 0.05)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, '.txt', sep='')
  m      <- read.table(filename, header = TRUE)
  
  ix <- 1:length(m$mu)
  group <- paste(m$mu, m$disp_cap, sep='-')
  ix <- tapply(ix, group, max)
  df <- rbind(df, m[ix,])

}

df <- df[df$disp_cap > 0,]

sdf <- summarySE(df, measurevar="perc_indiv", 
                 groupvars=c("mu","f_loss", "disp_cap"))

windows(height=4, width=7)

sdf$clustering <- 'Random'
sdf$clustering[sdf$mu == 3] <- 'Fractal'
sdf$clustering[sdf$mu == 5] <- 'Clustered'

ggplot(sdf, aes(x=f_loss, y=disp_cap, fill=perc_indiv/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering)) + 
  xlab('% Habitat loss') + 
  ylim(c(0, 0.6)) + 
  ylab('Dispersal strategy') + 
  scale_fill_gradientn(colors = c('white','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       labels = scales::percent_format(),
                       #trans = 'log10',
                       name="% Individuals",
                       na.value = 'white',
                       limits=c(0, 0.5)) + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 10)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

max(sdf$perc_indiv)

