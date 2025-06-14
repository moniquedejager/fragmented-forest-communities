library(ggplot2)
library(ggpubr)

filename <- 'results/subcommunity_data/fragmented_subcommunity_data_0.05immediate_after_hab_loss.txt'
df <- read.table(filename, header = TRUE)

for (i in seq(0.1, 0.95, 0.05)){
  filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, 'immediate_after_hab_loss.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    df <- rbind(df, m)
    
    #df <- m[m$dispersal == 'different',]
    #write.table(df, filename, append = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

##### use averages instead:
source('./src/summarySE.R')
#df$METE_fit[df$METE_fit < 0] <- 0

# average number of species per subcommunity:
sdf <- summarySE(df, measurevar='n_species',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same Dispersal (SD)'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different Dispersal (DD)'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same Dispersal (SD)', 'Different Dispersal (DD)'))

m1 <- sdf

# average % ancestors from elsewhere:
sdf <- summarySE(df, measurevar='Pm',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same Dispersal (SD)'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different Dispersal (DD)'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same Dispersal (SD)', 'Different Dispersal (DD)'))

m2 <- sdf

# average dissimilarity:
filename <- 'results/dissimilarity/dissimilarity_data_1.txt'
df3 <- read.table(filename, header = TRUE)

for (i in 1:10){
  filename <- paste('results/dissimilarity/dissimilarity_data_', i, '.txt', sep='')
  df2 <- read.table(filename, header = TRUE)
  df3 <- rbind(df3, df2)
}
df4 <- data.frame(mu = df3$mu,
                  f_loss = df3$f_hab_loss,
                  dispersal = df3$dispersal,
                  y = df3$dissimilarity)
sdf <- summarySE(df4, measurevar='y',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same Dispersal (SD)'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different Dispersal (DD)'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same Dispersal (SD)', 'Different Dispersal (DD)'))

m3 <- sdf

# plot them together
names(m1)[names(m1) == 'n_species'] <- 'y'
names(m2)[names(m2) == 'Pm'] <- 'y'
m2$y <- 1 - m2$y

m1$ylab <- '# Species per subcommunity'
m2$ylab <- '% Ancestors from elsewhere'
m3$ylab <- 'Bray-Curtis dissimilarity'

sdf <- rbind(m1, m2, m3)

pd <- position_dodge(1) 
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
p1

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 4.tif', 
     width = 5, height = 7, units = 'in', res = 600)
p1
dev.off()

# For presentation: 
sel <- sdf$ylab == '# Species per subcommunity'
p1 <- ggplot(sdf[sel,], aes(x=f_loss * 100, y=y, color=mu2)) + 
  geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  facet_grid(cols=vars(disp_type), rows=vars(ylab), scales='free_y', switch='y') + 
  xlab('% Habitat loss') + 
  #ylim(c(0, 400)) + 
  ylab('') +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank(), 
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines"))
p1

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 4 no species only immediate.tif', 
     width = 5, height = 3, units = 'in', res = 600)
p1
dev.off()

sel <- sdf$ylab == '% Ancestors from elsewhere'
p1 <- ggplot(sdf[sel,], aes(x=f_loss * 100, y=y, color=mu2)) + 
  geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  facet_grid(cols=vars(disp_type), rows=vars(ylab), scales='free_y', switch='y') + 
  xlab('% Habitat loss') + 
  ylab('') +
  #ylim(c(0,1)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank(), 
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines"))
p1

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 4 connectivity only immediate.tif', 
     width = 5, height = 3, units = 'in', res = 600)
p1
dev.off()
