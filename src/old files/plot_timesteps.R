# check the data per iteration:
library(ggplot2)

sim_nr <- 1
clustering <- 1
f_loss <- 0.1

filename <- paste('results/data per iteration/sim_nr=', sim_nr,
                  'clustering=', clustering, 'f_loss=', f_loss, '.txt', sep='')

df <- read.table(filename, header=TRUE)
df$timestep <- 1:length(df$sim_nr)

ggplot(df, aes(x=timestep, y=m_species)) + geom_point()

ggplot(df, aes(x=timestep, y=total_species)) + geom_point()

ggplot(df, aes(x=timestep, y=m_Pm)) + geom_point()


# now, lets make the same figures, but for all three spatial configurations 
# of habitat destruction and for four levels of f_loss (0.1, 0.4, 0.7, 0.9):

df <- df[df$sim_nr == 10,]
for (clustering in c(1, 3, 5)){
  for (f_loss in seq(0.05, 0.95, 0.05)){#c(0.1, 0.4, 0.7, 0.9)){
    filename <- paste('results/data per iteration/sim_nr=', sim_nr,
                      'clustering=', clustering, 'f_loss=', f_loss, '.txt', sep='')
    
    df2 <- read.table(filename, header=TRUE)
    df2$timestep <- 1:length(df2$sim_nr)
    df <- rbind(df, df2)
  }
}

t_max <- min(tapply(df$timestep, paste(df$f_loss, df$clustering, sep='-'), max))

ggplot(df[df$timestep <= t_max,], aes(x=timestep, y=f_loss, fill=m_species)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering)) + 
  scale_fill_viridis_c()

ggplot(df[df$timestep <= t_max,], aes(x=timestep, y=f_loss, fill=total_species)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering)) + 
  scale_fill_viridis_c()

ggplot(df[df$timestep <= t_max,], aes(x=timestep, y=f_loss, fill=m_Pm)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering)) + 
  scale_fill_viridis_c()



p1 <- ggplot(df) + 
  geom_line(aes(x=timestep, y=m_species, color=as.factor(clustering))) + 
  facet_grid(rows = vars(f_loss)) + 
  xlab('Iteration number') + 
  ylab('Mean number of species per subcommunity') + 
  theme_bw()

p2 <- ggplot(df) + 
  geom_line(aes(x=timestep, y=total_species, color=as.factor(clustering))) + 
  facet_grid(rows = vars(f_loss)) + 
  xlab('Iteration number') + 
  ylab('Total number of species') + 
  theme_bw()

ggplot(df, aes(x=timestep, y=total_species)) + 
  geom_point() + 
  facet_grid(cols = vars(f_loss), rows=vars(clustering), scales = 'free_x') + 
  xlab('Iteration number') + 
  ylab('Total number of species') + 
  theme_bw()

ggplot(df, aes(x=timestep, y=m_Pm)) + 
  geom_point() + 
  facet_grid(cols = vars(f_loss), rows=vars(clustering), scales = 'free_x') + 
  xlab('Iteration number') + 
  ylab('Mean dispersal capacity') + 
  theme_bw()


