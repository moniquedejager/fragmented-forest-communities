# truncated pareto distribution: 

lmin = 0.01
lmax = 5.5

df <- expand.grid(mu = 1:1000/1000,
                  x  = 1:1000/1000)


df$y = ((1 - df$x*(1 - (lmin/lmax)^df$mu))/lmin^df$mu)^(-1/df$mu)

p <- ggplot(df, aes(x=mu, y=x, fill=as.factor(round(y)))) + 
  geom_raster() + 
  scale_fill_viridis_d(name = 'Dispersal distance (d)',
                       guide=guide_legend(nrow=1)) + 
  xlab('Dispersal exponent (α)') + 
  ylab('Cumulative distribution function (Pr(D <= d))') + 
  theme_bw() + 
  theme(legend.position = 'top')
 
# 6 x 5:
tiff(filename = 'Fragmented-forest-communities/x64/Release/Figures/Suppl Figure Bounded pareto distribution.tif', 
     width = 6, height = 5, units = 'in', res = 600)
p
dev.off() 
