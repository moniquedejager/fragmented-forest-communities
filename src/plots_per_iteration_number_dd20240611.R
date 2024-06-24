sim_nr <- 1
mutation_rate <- 0.0003
filename <- paste('results/data per iteration/sim_nr=', sim_nr, 
                  'clustering=1f_loss=0.1mutation_rate=', mutation_rate,'.txt', sep='')
df <- read.table(filename, header=TRUE)   
df <- df[df$sim_nr == 2,]

f_loss <- seq(0.05, 0.95, 0.05)
clustering <- c(1, 3, 5)
sim_nr <- 1:10

for (i in f_loss){
  for (j in clustering){
    for (k in sim_nr){
      filename <- paste('results/data per iteration/sim_nr=', k, 
                        'clustering=',j,
                        'f_loss=', i,
                        'mutation_rate=', mutation_rate,'.txt', sep='')
      df2 <- read.table(filename, header=TRUE)  
      df2$dispersal <- 'different'
      filename2 <- paste('results/data per iteration/sim_nr=', k, 
                        'clustering=',j,
                        'f_loss=', i,
                        'mutation_rate=', mutation_rate,'different.txt', sep='')
      #write.table(df2, filename2, append=FALSE, 
      #            row.names = FALSE, col.names = TRUE)
    }
  }
}




library(ggplot2)

df$f_loss2 <- paste(df$f_loss*100, '% Habitat loss', sep='')
df$f_loss2 <- factor(df$f_loss2, levels = unique(df$f_loss2))
df$clustering2 <- 'Random'
df$clustering2[df$clustering == 3] <- 'Fractal'
df$clustering2[df$clustering == 5] <- 'Clustered'

new_df <- data.frame(iteration_nr = rep(df$iteration_nr, 4),
                     f_loss = rep(df$f_loss, 4),
                     clustering = rep(df$clustering2, 4),
                     y = c(rep('Mean # species per subcommunity', length(df$sim_nr)),
                           rep('Total # species', length(df$sim_nr)),
                           rep('Mean dispersal strategy', length(df$sim_nr)),
                           rep('% Ancestors from same subcommunity', length(df$sim_nr))),
                     z = c(df$m_species, df$total_species, df$m_Pm,df$f_t50_same_subcom*100),
                     disp = rep(df$dispersal, 4))
new_df$y <- factor(new_df$y, levels=unique(new_df$y))

sel = new_df$disp == 'similar'
ggplot(new_df[sel,], aes(x=iteration_nr, y=z, color=f_loss)) + 
  geom_point(alpha = 0.3)  +
  scale_color_viridis_c(name='% Habitat loss') + 
  facet_grid(cols=vars(clustering), 
             rows=vars(y), 
             scales='free_y',
             switch = 'y') +  
  xlab('Iteration number') + 
  ylab('') + 
  theme_bw() +
  guides(color = guide_colorbar(barwidth = 25)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
