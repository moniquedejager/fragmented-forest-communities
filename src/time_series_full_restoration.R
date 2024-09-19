library(ggplot2)
source('src/summarySE.R')

dat <- expand.grid(n_ind = 1000, 
                   Pm_range = 1, 
                   clustering = c(1,3,5), 
                   mutation_rate = 0.0003, 
                   max_mutation = 0,  
                   sim_nr = 1:10, #1:10,
                   f_loss = 0.95, #round(seq(0.05, 0.95, 0.05), 2),
                   dispersal = 'different', #c('similar', 'different'),
                   hab_cover = round(c(seq(0.05, 0.95, 0.15), 1), 2),#round(seq(0.05, 0.95, 0.05), 2),
                   clustering_restored = c(1, 3, 5),
                   n_iterations = 50)

df <- read.table('results/data per iteration/sim_nr=2clustering=1f_loss=0.95mutation_rate=3e-04restoration1b.txt',
                 header=TRUE)
df <- df[df$sim_nr == 1,]
for (i in round(c(seq(0.05, 0.95, 0.15), 1), 2)){ # habitat cover
  for (j in 1:10){ # simulation nr
    for (k in c(1, 3, 5)){ # clustering of destruction
      filename <- paste('results/data per iteration/sim_nr=',
                        j, 'clustering=',
                        k, 'f_loss=0.95mutation_rate=3e-04restoration',
                        i, 'b.txt', sep='')
      df2      <- read.table(filename, header=TRUE)
      df       <- rbind(df, df2) 
    }
  }
}


sdf <- summarySE(df, measurevar=c('m_METE_fit'),
                 groupvars=c('clustering', 'hab_cover', 'clustering2', 'iteration_nr'))

p1 <- ggplot(sdf, aes(x=iteration_nr, y=m_METE_fit, color=as.factor(hab_cover))) + 
  geom_line() + 
  #geom_point(alpha=0.5, size=0.6) + 
  facet_grid(cols=vars(clustering2), rows=vars(clustering)) + 
  #scale_x_continuous(trans='log10') + 
  scale_color_viridis_d(name='Restored habitat cover (%)') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('Average METE fit to SAD') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p1  


df$clustering2 <- 'Random restoration'
df$clustering2[df$clustering_restored == 3] <- 'Fractal restoration'
df$clustering2[df$clustering_restored == 5] <- 'Clustered restoration'

sel <- df$clustering > 0
p1 <- ggplot(df[sel,], aes(x=iteration_nr, y=m_species, color=as.factor(hab_cover))) + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2), rows=vars(clustering)) + 
  #scale_x_continuous(trans='log10') + 
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


# how does this relate to the average number of species per subcommunity in undisturbed systems?
df2 <- read.table('results/simulation_data/fragmented_simulation_data_0.txt', 
                  header=TRUE)
df2 <- df2[df2$dispersal == 'different',]

p1 <- ggplot(df[sel,], aes(x=iteration_nr, y=m_species, color=as.factor(hab_cover))) + 
  geom_hline(yintercept=min(df2$m_nspecies), linetype='dashed', color='grey') + 
  geom_hline(yintercept=max(df2$m_nspecies), linetype='dashed', color='grey') + 
  geom_hline(yintercept=mean(df2$m_nspecies), color='grey30') + 
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

# how long does it take to get back to the minimum number of species per subcommunity?
df3 <- df[df$m_species >= min(df2$m_nspecies),]
tapply(df3$iteration_nr, df3$clustering2, min)

df3 <- df[df$m_species >= mean(df2$m_nspecies),]
tapply(df3$iteration_nr, df3$clustering2, min)

p2 <- ggplot(df[sel,], aes(x=iteration_nr, y=total_species, color=as.factor(hab_cover))) + 
  geom_hline(yintercept=min(df2$n_species), linetype='dashed', color='grey') + 
  geom_hline(yintercept=max(df2$n_species), linetype='dashed', color='grey') + 
  geom_hline(yintercept=mean(df2$n_species), color='grey30') + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2), rows=vars(clustering)) + 
  scale_color_viridis_d(name='Restored habitat cover (%)') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('Total number of species') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p2 

# how long does it take to get back to the minimum total number of species?
df3 <- df[df$total_species >= min(df2$n_species),]
tapply(df3$iteration_nr, df3$clustering2, min)

df3 <- df[df$total_species >= mean(df2$n_species),]
tapply(df3$iteration_nr, df3$clustering2, min)

p1 <- ggplot(df, aes(x=iteration_nr, y=m_METE_fit, color=as.factor(hab_cover))) + 
  geom_point(alpha=0.15, size=0.6) + 
  facet_grid(cols=vars(clustering2), rows=vars(clustering)) + 
  #scale_x_continuous(trans='log10') + 
  scale_color_viridis_d(name='Restored habitat cover (%)') + 
  xlab('Time since restoration (# model iterations)') + 
  ylab('Average METE fit to SAD') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "inches"),
        panel.spacing.x = unit(1, "lines")) + 
  guides(color = guide_legend(nrow=1, override.aes = list(size = 6, alpha=1))) 
p1  
