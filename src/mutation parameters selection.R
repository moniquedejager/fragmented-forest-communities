library(ggplot2)
library(ggpubr)

filename <- 'results/simulation results 03-2024/simulation_data/simulation_data_1.txt'
df       <- read.table(filename, header=T)

p1 <- ggplot(df, aes(x=as.factor(mutation_rate), y=max_mutation, fill=m_METE_fit)) + 
  geom_raster() + 
  scale_fill_gradient2(low='red', mid='white', high='green', 
                       name='METE fit', midpoint = mean(df$m_METE_fit)) + 
  facet_wrap(vars(mu)) + 
  xlab('Mutation rate') + 
  ylab('Maximum mutation size')
p1

p2 <- ggplot(df, aes(x=as.factor(mutation_rate), y=max_mutation, fill=n_iterations)) + 
  geom_raster() + 
  scale_fill_gradient2(low='green', mid='white', high='red', 
                       name='# Iterations', midpoint = mean(df$n_iterations)) + 
  facet_wrap(vars(mu)) + 
  xlab('Mutation rate') + 
  ylab('Maximum mutation size')
p2

p3 <- ggplot(df, aes(x=as.factor(mutation_rate), y=max_mutation, fill=m_nspecies)) + 
  geom_raster() + 
  scale_fill_gradient2(low='red', mid='white', high='green', 
                       name='Mean # species',  midpoint=mean(df$m_nspecies)) + 
  facet_wrap(vars(mu)) + 
  xlab('Mutation rate') + 
  ylab('Maximum mutation size')
p3

p4 <- ggplot(df, aes(x=as.factor(mutation_rate), y=max_mutation, fill=n_species)) + 
  geom_raster() + 
  scale_fill_continuous(trans='log10', name='Total # species') + 
  facet_wrap(vars(mu)) + 
  xlab('Mutation rate') + 
  ylab('Maximum mutation size')
p4

range(df$n_species)
ggplot(df, aes(x=as.factor(mutation_rate), y=n_species)) + 
  geom_violin() + 
  scale_y_continuous(trans='log10')

ggarrange(p1, p2, p3, p4, nrow=4, ncol=1, align='hv')


# is there an effect of distance to center on n_species and METE fit?
filename <- 'results/simulation results 03-2024/subcommunity_data/subcommunity_data_1.txt'
df2 <- read.table(filename, header=T)
df2 <- df2[(df2$mutation_rate == 0.001)&(df2$max_mutation == 0.05),]

df2$dist_to_center <- sqrt((df2$x - 10)^2 + (df2$y - 10)^2)
ggplot(df2, aes(x=x, y=y, fill=dist_to_center)) + 
  geom_raster()

ggplot(df2, aes(x=dist_to_center, y=n_species, color=x)) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mu))

ggplot(df2, aes(x=dist_to_center, y=METE_fit)) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mu))

ggplot(df2, aes(x=dist_to_center, y=m_dist)) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mu))

ggplot(df2, aes(x=n_species, y=Shannon)) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mu))

ggplot(df2, aes(x=x, y=y, fill=METE_fit)) + 
  geom_raster() + 
  facet_wrap(vars(mu))

ggplot(df2, aes(x=x, y=y, fill=n_species)) + 
  geom_raster() + 
  facet_wrap(vars(mu))

ggplot(df2, aes(x=x, y=y, fill=m_dist)) + 
  geom_raster() + 
  facet_wrap(vars(mu))

x = sort(rep(-4:25, 30))
y = rep(-4:25, 30)

ID = 1:900

ID2 = (x - min(x))*(30) + (y - (min(y))) + 1
plot(ID, ID2)






