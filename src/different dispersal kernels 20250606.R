# show the three different dispersal kernels next to each other:
df <- data.frame(dx = rep(-10:10, length(-10:10)))
df$dy   <- sort(df$dx)
df$dist <- sqrt(df$dx^2 + df$dy^2)

df2 <- df
df2$kernel <- 0
df2$P      <- 0
df2$a      <- 0

rho  = 0 # parameter for 2-D gaussian distribution
lmin = 1 # parameter for 2-D pareto distribution
lmax = 20 # parameter for 2-D pareto distribution

for (a in c(0.1, 0.5, 1)){
  # 2-D exponential:
  P <- 1/(2*pi*a) * exp(-df$dist/a)
  P <- P/sum(P)
  
  df3        <- df
  df3$kernel <- '2-D exponential'
  df3$P      <- P
  df3$a      <- a
  
  # 2-D Gaussian:
  P = 1/(2*pi) * exp(-1/2 * ((df$dx/(a*2))^2 + (df$dy/(a*2))^2))
  P2 = 1/(2*pi) * exp(-1/2 * (df$dist^2/(a*2)^2))
  P <- P/sum(P)
  
  df4        <- df
  df4$kernel <- '2-D Gaussian'
  df4$P      <- P
  df4$a      <- a
  
  # 2-D pareto:
  mu = 1/a
  mu[mu == 2] <- 2.001
  
  P <- 1/(2*pi) * (2 - mu)/(lmax^(2-mu) - lmin^(2-mu)) * (df$dist+1)^(-mu)
  P <- P/sum(P)
  
  df5        <- df
  df5$kernel <- '2-D Pareto'
  df5$P      <- P
  df5$a      <- a
  
  df2 <- rbind(df2, df3, df4, df5)
}
df2 <- df2[df2$kernel != '0',]

sel <- df2$dy == 0
df2$P2 <- df2$P
df2$P2[df2$P2 < 0.000001] <- 0

cbp1 <- c("#009E73", "#E69F00", "#999999")
df2$a2 <- paste('λ = ', df2$a, sep='')
p1 <- ggplot(df2[sel,], aes(x=dx, y=P2)) + 
  geom_line(aes(col=a2, linetype=a2),linewidth=1.1) + 
  facet_wrap(vars(kernel)) + 
  scale_y_continuous(trans='log10') + 
  scale_color_manual(values=cbp1) + 
  xlab('Δx') + 
  ylab('Dispersal probability (at Δy = 0)') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())
p1

# 9x4
# tiff file 600 dpi:
tiff(filename = 'figures/plot different dispersal kernels 20250606.tif', 
     width = 9, height = 4, units = 'in', res = 600)
p1
dev.off()


