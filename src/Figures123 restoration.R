library(ggplot2)

df <- read.table('results/data per iteration/sim_nr=2clustering=1f_loss=0.95mutation_rate=3e-04restoration1.txt',
                 header=TRUE)
df <- df[df$sim_nr == 1,]
for (i in round(c(seq(0.05, 0.95, 0.15), 1), 2)){
  for (j in 1:10){
    for (k in c(1, 5)){
      filename <- paste('results/data per iteration/sim_nr=',
                        j, 'clustering=',
                        k, 'f_loss=0.95mutation_rate=3e-04restoration',
                        i, '.txt', sep='')
      df2      <- read.table(filename, header=TRUE)
      df       <- rbind(df, df2) 
    }
  }
}

df$clustering2 <- 'Random restoration'
df$clustering2[df$clustering_restored == 3] <- 'Fractal restoration'
df$clustering2[df$clustering_restored == 5] <- 'Clustered restoration'

sel <- df$clustering == 5
p1 <- ggplot(df[sel,], aes(x=iteration_nr, y=m_species, color=as.factor(hab_cover))) + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2)) + 
  scale_color_viridis_d(name='Restored habitat cover (%)') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('Average number of species per subcommunity') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p1  

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 1 restoration.tif', 
     width = 10, height = 5, units = 'in', res = 600)
p1
dev.off()


p2 <- ggplot(df[sel,], aes(x=iteration_nr, y=total_species, color=as.factor(hab_cover))) + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2)) + 
  scale_color_viridis_d(name='Restored habitat cover') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('Total number of species in metacommunity') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p2

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 2 restoration.tif', 
     width = 10, height = 5, units = 'in', res = 600)
p2
dev.off()

p3 <- ggplot(df[sel,], aes(x=iteration_nr, y=m_Pm, color=as.factor(hab_cover))) + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2)) + 
  scale_color_viridis_d(name='Restored habitat cover') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('Average dispersal capacity') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p3

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 3 restoration.tif', 
     width = 10, height = 5, units = 'in', res = 600)
p3
dev.off()


# figures 1-3 together:

df2 <- data.frame(iteration_nr = rep(df$iteration_nr, 3),
                  y = c(df$m_species, df$total_species, df$m_Pm),
                  hab_cover = rep(df$hab_cover, 3),
                  clustering = rep(df$clustering, 3),
                  clustering2 = rep(df$clustering2, 3),
                  y_labs = c(rep('Average number of species per subcommunity', length(df$sim_nr)),
                             rep('Total number of species in metacommunity', length(df$sim_nr)),
                             rep('Average dispersal capacity', length(df$sim_nr))))
df2$y_labs <- factor(df2$y_labs, levels=unique(df2$y_labs))

sel <- df2$clustering == 5
p4 <- ggplot(df2, aes(x=iteration_nr, y=y, color=as.factor(hab_cover))) + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2), rows=vars(y_labs), scales='free_y', switch='y') + 
  scale_color_viridis_d(name='Restored habitat cover') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p4

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 1-3 restoration_clustered destruction.tif', 
     width = 10, height = 9, units = 'in', res = 600)
p4
dev.off()
