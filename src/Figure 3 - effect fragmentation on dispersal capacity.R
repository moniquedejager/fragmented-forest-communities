# How does the dispersal capacity change with increasing habitat fragmentation?
library(ggplot2)
source('./src/summarySE.R')

filename <- './results/dispersal_capacity/fragmented_dispersal capacity_data_0.txt'
df       <- read.table(filename, header = TRUE)
df       <- df[df$mutation_rate == 0.0003,]
df       <- df[df$n_iterations == 10050,]
df       <- df[40:49,]
df2      <- df
df2$mu   <- 3
df3      <- df
df3$mu   <- 5
df       <- rbind(df, df2, df3)

for (i in seq(0.05, 0.95, 0.05)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, '.txt', sep='')
  m        <- read.table(filename, header = TRUE)
  m        <- m[m$mutation_rate == 0.0003,]
  
  df    <- rbind(df, m)
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
  ylim(c(0, 1)) + 
  ylab('Dispersal strategy') + 
  scale_fill_gradientn(colors = c('white','yellow', 'seagreen', 
                                  'darkcyan', 'midnightblue'),
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

max(sdf$perc_indiv)

