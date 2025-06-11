library(ggplot2)

df        <- read.table('results/example of seed dispersal.txt')
names(df) <- c('n_cell', 'n_indiv', 'n_seeds', 
               'f_hab', 'simnr', 'strategy', 'generation')

ggplot(df, aes(x=n_cell, y=n_indiv, color=strategy)) + 
  geom_jitter() + 
  facet_grid(cols=vars(n_seeds), rows=vars(f_hab)) + 
  scale_color_viridis_c()

group <- paste(df$n_cell, df$n_indiv, df$n_seeds, df$f_hab, sep='-')
sdf <- data.frame(n_cell = tapply(df$n_cell, group, mean),
                  n_indiv = tapply(df$n_indiv, group, mean),
                  n_seeds = tapply(df$n_seeds, group, mean),
                  f_hab = tapply(df$f_hab, group, mean),
                  strategy = tapply(df$strategy, group, mean),
                  generation = tapply(df$generation, group, mean),
                  sd         = tapply(df$strategy, group, sd),
                  min_strat  = tapply(df$strategy, group, min),
                  max_strat  = tapply(df$strategy, group, max))

sdf$d_strat <- sdf$max_strat - sdf$min_strat
ggplot(sdf, aes(x=factor(n_indiv), y=factor(n_cell), fill=d_strat)) + 
  geom_raster() + 
  facet_grid(cols=vars(f_hab), rows=vars(n_seeds)) + 
  scale_fill_viridis_c()

cols <- c('grey', 'pink', 'seagreen', 'darkblue')
ggplot(sdf, aes(x=factor(f_hab), y=factor(n_indiv), color=factor(round(strategy)))) + 
  geom_point(aes(size=(1-sd))) + 
  facet_grid(cols=vars(n_cell), rows=vars(n_seeds)) + 
  scale_color_manual(values = cols) + 
  theme_bw()

ggplot(sdf, aes(x=factor(f_hab), y=factor(n_indiv), fill=factor(round(strategy)))) + 
  geom_raster() + 
  facet_grid(cols=vars(n_cell), rows=vars(n_seeds)) + 
  scale_fill_manual(values = cols) + 
  theme_bw()

cbp1 <- c("#009E73", "#E69F00", "#999999")
df$n_cell2 <- paste(df$n_cell, ' cells')
df$n_cell2 <- factor(df$n_cell2, levels=unique(df$n_cell2))
df$n_indiv2 <- paste(df$n_indiv, ' indiv. per cell')
df$n_indiv2 <- factor(df$n_indiv2, levels=unique(df$n_indiv2))
p1 <- ggplot(df, aes(x=(1 - f_hab)*100, y=strategy, color=factor(n_seeds))) + 
  geom_line(aes(linetype=factor(simnr))) + 
  geom_point() + 
  facet_grid(cols=vars(n_cell2), rows=vars(n_indiv2)) + 
  scale_color_manual(values=cbp1, name='# Seeds per indiv.') + 
  xlab('% Habitat loss') + 
  ylab('Dispersal strategy') + 
  theme_bw() +
  guides(linetype = "none") +  
  theme(strip.background = element_blank(),
        legend.position = 'top')

# 6x5
# tiff file 600 dpi:
tiff(filename = 'figures/example of seed dispersal 20250528.tif', 
     width = 6, height = 5, units = 'in', res = 600)
p1
dev.off()

# results: 
# in a continuous environment, long-distance dispersal is always more advantageous. 
# when the number of individuals per cell is sufficient (i.e. the patch is large
# enough), local dispersal takes over when (even a small amount of) habitat loss occurs. 

# these results show that local dispersal can easily evolve when habitat loss 
# occurs, but habitat patches remain sufficiently large (the negative effects of
# kin competition are smaller than the negative effects of potentially dispersing
# to non-habitat areas. 


# also plot the dispersal strategies:
df2 <- data.frame(probs = c(1, 0, 0, 
                            0.8, 0.2, 0,
                            0.6, 0.3, 0.1,
                            0.33, 0.33, 0.33),
                  strategy = c(1, 1, 1, 
                               2, 2, 2,
                               3, 3, 3, 
                               4, 4, 4),
                  distance = c(1, 2, 3,
                               1, 2, 3,
                               1, 2, 3,
                               1, 2, 3) - 1)

df2$strategy2 <- paste('Strategy ', df2$strategy)
df2$strategy2 <- factor(df2$strategy2, levels= unique(df2$strategy2))
p2 <- ggplot(df2, aes(x=distance, y=probs)) + 
  geom_bar(stat = "identity", color='black', fill="#E69F00") + 
  facet_wrap(vars(strategy2)) + 
  xlab('Dispersal distance (cells)') + 
  ylab('Probability') + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        legend.position = 'top')
p2
# 4x3:
tiff(filename = 'figures/example of seed dispersal 20250528 b.tif', 
     width = 4, height = 3, units = 'in', res = 600)
p2
dev.off()

