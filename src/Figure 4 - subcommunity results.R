library(ggplot2)
library(ggpubr)

filename <- 'results/subcommunity_data/fragmented_subcommunity_data_0.05.txt'
df <- read.table(filename, header = TRUE)

for (i in seq(0.1, 0.95, 0.05)){
  filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    df <- rbind(df, m)
    
    #df <- m[m$dispersal == 'different',]
    #write.table(df, filename, append = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

df$METE_fit[df$METE_fit < 0] <- 0
#df$Shannon <- df$Shannon / log(df$n_species)

##### use averages instead:
source('./src/summarySE.R')
#df$METE_fit[df$METE_fit < 0] <- 0

# average number of species per subcommunity:
sdf <- summarySE(df, measurevar='n_species',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

m1 <- sdf

# average % ancestors from elsewhere:
sdf <- summarySE(df, measurevar='Pm',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

m2 <- sdf

# average METE fit to SAD per subcommunity:
sdf <- summarySE(df, measurevar='METE_fit',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

m3 <- sdf

# plot them together
names(m1)[names(m1) == 'n_species'] <- 'y'
names(m2)[names(m2) == 'Pm'] <- 'y'
m2$y <- 1 - m2$y
names(m3)[names(m3) == 'METE_fit'] <- 'y'

m1$ylab <- '# Species per subcommunity'
m2$ylab <- '% Ancestors from elsewhere'
m3$ylab <- 'METE fit to SAD'

sdf <- rbind(m1, m2, m3)

sdf$ylab <- factor(sdf$ylab, levels=unique(sdf$ylab))
p1 <- ggplot(sdf, aes(x=f_loss * 100, y=y, color=mu2)) + 
  geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  facet_grid(cols=vars(disp_type), rows=vars(ylab), scales='free_y', switch='y') + 
  xlab('% Habitat loss') + 
  ylab('') +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank(), 
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines"))

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 4.tif', 
     width = 5, height = 7, units = 'in', res = 600)
p1
dev.off()

