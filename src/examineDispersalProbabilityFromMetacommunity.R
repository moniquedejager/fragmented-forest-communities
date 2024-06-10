# examine what rate of metacommunity dispersal into the subcommunities
# works best to stabilize the subcommunities:

sim_nr <- 1
mutation_rate <- c(0, 0.0001, 0.0003, 0.0005, 0.001, 0.003, 
                  0.005, 0.01, 0.03, 0.05)
filename <- paste('results/data per iteration/sim_nr=', sim_nr, 
                  'clustering=1f_loss=0mutation_rate=', mutation_rate[1],'.txt', sep='')
df <- read.table(filename, header=TRUE)                  

for (i in 2:10){
  filename <- paste('results/data per iteration/sim_nr=', sim_nr, 
                    'clustering=1f_loss=0mutation_rate=', mutation_rate[i],'.txt', sep='')
  df2 <- read.table(filename, header=TRUE)  
  df <- rbind(df, df2)
}

for (i in 1:10){
  filename <- paste('results/data per iteration/sim_nr=', sim_nr, 
                    'clustering=1f_loss=0.95mutation_rate=', mutation_rate[i],'.txt', sep='')
  df2 <- read.table(filename, header=TRUE)  
  df <- rbind(df, df2)
}



library(ggplot2)
library(ggpubr)

df$f_loss2 <- paste(df$f_loss*100, '% Habitat loss', sep='')
df$f_loss2 <- factor(df$f_loss2, levels = unique(df$f_loss2))

p1 <- ggplot(df[], aes(x=iteration_nr, y=m_species, color=as.factor(mutation_rate))) + 
  geom_point() +
  scale_color_viridis_d(name='Dispersal probability from metacommunity') + 
  facet_grid(cols=vars(f_loss2)) + 
  xlab('Iteration number') + 
  ylab('Mean # species per subcommunity') + 
  theme_bw()
p1

p2 <- ggplot(df, aes(x=iteration_nr, y=total_species, color=as.factor(mutation_rate))) + 
  geom_point()  +
  scale_color_viridis_d(name='Dispersal probability from metacommunity') + 
  facet_grid(cols=vars(f_loss2)) + 
  xlab('Iteration number') + 
  ylab('Total # species') + 
  theme_bw()
p2

ggplot(df, aes(x=iteration_nr, y=f_t50_same_subcom, color=as.factor(mutation_rate))) + 
  geom_point()  +
  scale_color_viridis_d(name='Dispersal probability from metacommunity') + 
  facet_grid(cols=vars(f_loss2)) + 
  xlab('Iteration number') + 
  ylab('Total # species') + 
  theme_bw()

df2 <- read.table('results/subcommunity_data/fragmented_subcommunity_data_0.txt',
                  header=TRUE)
df3 <- read.table('results/subcommunity_data/fragmented_subcommunity_data_0.95.txt',
                  header=TRUE)
df2 <- rbind(df2, df3)

df2$f_loss2 <- paste(df2$f_loss*100, '% Habitat loss', sep='')
df2$f_loss2 <- factor(df2$f_loss2, levels = unique(df2$f_loss2))
p3 <- ggplot(df2, aes(x=as.factor(mutation_rate), y=METE_fit, color=as.factor(mutation_rate))) + 
  geom_boxplot() + 
  scale_color_viridis_d(name='Dispersal probability from metacommunity') + 
  facet_grid(cols=vars(f_loss2)) + 
  xlab('Dispersal probability from metacommunity') + 
  ylab('METE fit to SAD') + 
  theme_bw()
p3

ggarrange(p1, p2, p3, labels=c('A', 'B', 'C'), common.legend = TRUE, ncol=1, nrow=3)

tapply(df$iteration_nr, paste(df$mutation_rate, df$f_loss, sep='-'), max)

