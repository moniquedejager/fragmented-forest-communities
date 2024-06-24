# How does the dispersal capacity change with increasing habitat fragmentation?
library(ggplot2)
source('./src/summarySE.R')

filename <- './results/dispersal_capacity/fragmented_dispersal capacity_data_0.05.txt'
df       <- read.table(filename, header = TRUE)

for (i in seq(0.1, 0.95, 0.05)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, '.txt', sep='')
  m        <- read.table(filename, header = TRUE)
  m        <- m[m$mutation_rate == 0.0003,]
  
  df    <- rbind(df, m)
}
df <- df[df$dispersal == 'different',]

sdf <- summarySE(df, measurevar="perc_indiv", 
                 groupvars=c("mu","f_loss", "disp_cap"))

windows(height=4, width=7)

sdf$clustering <- 'Random'
sdf$clustering[sdf$mu == 3] <- 'Fractal'
sdf$clustering[sdf$mu == 5] <- 'Clustered'

sdf$perc_indiv[sdf$perc_indiv == 0] <- NA
ggplot(sdf, aes(x=f_loss*100, y=disp_cap, fill=perc_indiv/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering)) + 
  xlab('% Habitat loss') + 
  ylim(c(0, 1)) + 
  ylab('Dispersal strategy') + 
  scale_fill_gradientn(colors = c('lightyellow','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue'),
                       labels = scales::percent_format(),
                       #trans = 'log10',
                       name="% Individuals",
                       na.value = 'transparent',
                       limits=c(0.1, 1)) + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 10)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

max(sdf$perc_indiv)

