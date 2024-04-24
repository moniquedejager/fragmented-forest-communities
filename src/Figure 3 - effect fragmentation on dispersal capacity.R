# How does the dispersal capacity change with increasing habitat fragmentation?
library(ggplot2)
source('./src 03-2024/summarySE.R')

filename <- './results/simulation results 04-2024/dispersal_capacity/dispersal capacity_data_1.txt'
df       <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('./results/simulation results 04-2024/dispersal_capacity/dispersal capacity_data_', i, '.txt', sep='')
  df2      <- read.table(filename, header = TRUE)
  df       <- rbind(df, df2)
}

df <- df[df$disp_cap > 0,]

sdf <- summarySE(df, measurevar="perc_indiv", 
                 groupvars=c("mu","f_hab_loss", "disp_cap"))

windows(height=4, width=7)

sdf$clustering <- 'Random'
sdf$clustering[sdf$mu == 3] <- 'Fractal'
sdf$clustering[sdf$mu == 5] <- 'Clustered'

ggplot(sdf, aes(x=f_hab_loss, y=disp_cap, fill=perc_indiv/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering)) + 
  xlab('% Habitat loss') + 
  ylab('Dispersal strategy') + 
  scale_fill_gradientn(colors = c('white','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue', 'black'),
                       labels = scales::percent_format(),
                       #trans = 'log10',
                       name="% Individuals",
                       na.value = 'white',
                       limits=c(0, 1)) + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 10)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())


