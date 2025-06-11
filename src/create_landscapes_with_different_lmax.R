# Patch destruction using weighted sampling. 
# here, we examine the effect of lmax on the spatial pattern. 
# in the main model, we use lmax = 50

data <- data.frame(lmax = vector(length=0),
                   mu = vector(length=0),
                   f_loss = vector(length=0), 
                   moransI = vector(length=0))

for (lmax in c(10, 30, 50, 100)){
  for (mu in seq(1, 5, 0.2)){
    for (simnr in 1:5){
      print(paste('mu = ', mu, 'and lmax = ', lmax))
      
      # start with a uniform, continuous habitat and fragment it...
      if (mu == 2){ mu <- 2.01}
      n <- 20
      x <- rep(1:n, n)
      y <- sort(x)
      h <- rep(1, n*n)
      
      mat_h <- h
      
      # first destroyed patch:
      destroyed    <- sample(1:(n*n), 1)
      h[destroyed] <- 0
      
      # weights for sampling patch destruction:
      w <- rep(0, n*n)
      
      # 2D pareto distribution of weights around the destroyed patch:
      lmin <- 1
      #lmax <- 50
      
      for (f in seq(0.95, 0.05, -0.05)){
        print(paste('f = ', f))
        
        n_patches <- f * n * n  # number of habitat patches that should remain
        while (sum(h) > n_patches){
          dist <- sqrt((x - x[destroyed])^2 + (y - y[destroyed])^2)
          w    <- w + 1/(2*pi) * ((2 - mu) / (lmax^(2-mu) - lmin^(2-mu))) * dist^(-1*mu)
          w    <- w * h
          w[is.na(w)] <- 0  # focal patch
          
          destroyed    <- sample(1:(n*n), 1, prob = w)
          h[destroyed] <- 0
        }
        
        # write the data to a dataframe and then to a file:
        mat_h <- cbind(mat_h, h)
        
        # calculate Moran's I
        # create vectors of the habitat values of each neighbor:
        zmean      <- mean(h)
        d0         <- h - zmean
        
        x2         <- x - 1
        x2[x2 < 1] <- n
        ID2        <- (y - 1)*n + x2
        n1         <- h[ID2]
        d1         <- n1 - zmean
        
        x2         <- x + 1
        x2[x2 > n] <- 1
        ID2        <- (y - 1)*n + x2
        n2         <- h[ID2]
        d2         <- n2 - zmean
        
        y2         <- y - 1
        y2[y2 < 1] <- n
        ID2        <- (y2 - 1)*n + x
        n3         <- h[ID2]
        d3         <- n3 - zmean
        
        y2         <- y + 1
        y2[y2 > n] <- 1
        ID2        <- (y2 - 1)*n + x
        n4         <- h[ID2]
        d4         <- n4 - zmean
        
        # Moran's I = 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2))  
        I <- 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2)) 
        
        data2 <- data.frame(lmax = lmax,
                            mu = mu,
                            f_loss = f, 
                            moransI = I)
        data <- rbind(data, data2)
      }
      mu[mu == 2.01] <- 2
      filename <- paste('./results/revision/landscapes_lmax', lmax, '_', mu, '.txt', sep='')
      write.table(mat_h, filename, append=FALSE, col.names = FALSE, row.names = FALSE)
    }
  }
}
filename <- './results/revision/moransI_lmax.txt'
write.table(data, filename, append=FALSE, col.names = TRUE, row.names = FALSE)

group <- paste(data$lmax, data$mu, data$f_loss, sep='-')
df <- data.frame(lmax = tapply(data$lmax, group, mean),
                 mu   = tapply(data$mu, group, mean),
                 f_cover = tapply(data$f_loss, group, mean),
                 moransI = tapply(data$moransI, group, mean))

df$lmax2 <- paste('Lmax = ', df$lmax)
df$lmax2 <- factor(df$lmax2, levels=unique(df$lmax2)[c(1, 3, 4,2)])
p1 <- ggplot(df, aes(x=round(f_cover*100), y=round(mu, 1), fill=moransI)) + 
  geom_raster() + 
  facet_wrap(vars(lmax2))+ 
  scale_fill_viridis_c(name="Moran's I", limits=c(-0.05, 1)) + 
  xlab('% Habitat cover') + 
  ylab(expression(mu)) + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 20)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

p1

sel <- (df$lmax %in% c(10, 50))
summary(lm(moransI~f_cover+mu+I(mu^2)+lmax2, data=df[sel,]))
# lmax doesn't matter (p = 0.261)

sel <- (df$lmax %in% c(30, 50))
summary(lm(moransI~f_cover+mu+I(mu^2)+lmax2, data=df[sel,]))
# lmax doesn't matter

sel <- (df$lmax %in% c(100, 50))
summary(lm(moransI~f_cover+mu+I(mu^2)+lmax2, data=df[sel,]))
# lmax doesn't matter (p = 0.261)

summary(lm(moransI~f_cover+mu+I(mu^2)+lmax2, data=df))
# lmax doesn't matter

#7x6:
# tiff file 600 dpi:
tiff(filename = 'figures/create landscapes with different lmax.tif', 
     width = 7, height = 6, units = 'in', res = 600)
p1
dev.off()
