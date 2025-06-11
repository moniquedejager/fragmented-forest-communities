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

# plot the average dispersal strategy:
df$y <- df$perc_indiv/10 * df$disp_cap

sdf <- summarySE(df, measurevar="y", 
                 groupvars=c("mu","f_loss"))
sdf$clustering <- 'Random'
sdf$clustering[sdf$mu == 3] <- 'Fractal'
sdf$clustering[sdf$mu == 5] <- 'Clustered'

pd <- position_dodge(0.0) 

p1 <- ggplot(sdf, aes(x=f_loss*100, y=y)) + 
  geom_bar(stat='identity', color='black', fill="#E69F00", alpha=0.7) +
  geom_errorbar(aes(ymin=y-se, ymax=y+se), width=0, position=pd) + 
  facet_wrap(vars(clustering)) + 
  xlab('% Habitat loss') + 
  ylab(expression(paste('Average dispersal capacity ( ', lambda, ')'))) +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())
p1  

# 8x3
# tiff file 600 dpi:
tiff(filename = 'figures/Figure 2 dd20250603.tif', 
     width = 8, height = 3, units = 'in', res = 600)
p1
dev.off()
  
