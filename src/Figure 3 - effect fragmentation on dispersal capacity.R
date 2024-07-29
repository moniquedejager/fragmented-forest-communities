# How does the dispersal capacity change with increasing habitat fragmentation?
library(ggplot2)
source('./src/summarySE.R')

start <- TRUE
for (i in seq(0, 0.95, 0.05)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m        <- read.table(filename, header = TRUE)
    m        <- m[m$mutation_rate == 0.0003,]
    
    if (start == TRUE){
      df <- m
      start <- FALSE
    } else {
      df    <- rbind(df, m)
    }
  }
}
df <- df[df$dispersal == 'different',]

sdf <- summarySE(df, measurevar="perc_indiv", 
                 groupvars=c("mu","f_loss", "disp_cap"))

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

table(df$mu, df$f_loss, df$sim_nr)

# plot the average dispersal strategy:
df$y <- df$perc_indiv/10 * df$disp_cap

sdf <- summarySE(df, measurevar="y", 
                 groupvars=c("mu","f_loss", "sim_nr"))

sdf <- summarySE(sdf, measurevar="y", 
                 groupvars=c("mu","f_loss"))
sdf$clustering <- 'Random'
sdf$clustering[sdf$mu == 3] <- 'Fractal'
sdf$clustering[sdf$mu == 5] <- 'Clustered'

pd <- position_dodge(0.01) 
p1 <- ggplot(sdf, aes(x=f_loss*100, y=y, color=clustering)) + 
  geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('% Habitat loss') + 
  ylab(expression(paste('Average dispersal capacity ( ', lambda, ')'))) +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 3.tif', 
     width = 4, height = 3, units = 'in', res = 600)
p1
dev.off()
  
