library(ggplot2)
library(ggpubr)

# read in the community composition data and calculate part of the species 
# area curve. Fit the simulation results to the METE estimation.

SARs <- data.frame(f_loss     = vector(length=0), # fraction of habitat lost
                   clustering = vector(length=0), # clustering of habitat (mu)
                   sim_nr     = vector(length=0), # simulation number
                   area       = vector(length=0), # area size (in km2)
                   species1   = vector(length=0), # number of species in this area size, regardless of the amount of habitable patches 
                   species2   = vector(length=0), # number of species in habitable area of this size
                   type       = vector(length=0)) # Simulation output or METE prediction

GOFs <- data.frame(f_loss     = vector(length=0), # fraction of habitat lost
                   clustering = vector(length=0), # clustering of habitat patches (mu)
                   sim_nr     = vector(length=0), # simulation number
                   GOF1       = vector(length=0), # goodness-of-fit when disregarding the amount of habitable patches available
                   GOF2       = vector(length=0)) # goodness-of-fit when using habitable area size as area size

dat            <- expand.grid(x=0:44, y=0:44)
uPatchID       <- paste(dat$x, dat$y, sep='-')
areas          <- c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6)

filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition/')

if (!("GOFs.txt" %in% filenames)){
  for (j in 1:length(filenames))
  {
    filename <- filenames[j]
    
    if (grepl("composition", filename))
    {
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
      df <- read.table(filename2)
      
      names(df) <- c('x', 'y', 'species', 'n')
      df <- df[df$x != -1,]
      
      patchID <- paste(df$x, df$y, sep='-')
      uPatch2 <- unique(patchID)
      
      # For the 11 area sizes, find how many species are in the simulated
      # landscape:
      spec_observed1 <- vector(length=0)
      spec_observed2 <- vector(length=0)
     
      for (i in areas / 0.025)
      {
        nspecies       <- length(unique(df$species[patchID %in% uPatchID[1:i]]))
        spec_observed1 <- c(spec_observed1, nspecies)
        
        nspecies       <- length(unique(df$species[patchID %in% uPatch2[1:i]]))
        spec_observed2 <- c(spec_observed2, nspecies)
      }
      
      N01            <- sum(df$n[patchID %in% uPatchID[1:i]]) # total number of individuals 
      N02            <- sum(df$n[patchID %in% uPatch2[1:i]])
      
      sar1 <- meteSAR(Amin=0.025, 
                      A0=25.6, 
                      S0=spec_observed1[11], 
                      N0=N01) # N0 klopt niet, aanpassen!! 
      spec_predicted1 <- sar1$pred$S
      
      sar2 <- meteSAR(Amin=0.025, 
                      A0=25.6, 
                      S0=spec_observed2[11], 
                      N0=N02)
      spec_predicted2 <- sar2$pred$S
      
      # goodness-of-fit:
      GOF1 <- 1 - sum((spec_observed1-spec_predicted1)^2)/sum((spec_observed1-mean(spec_predicted1))^2)
      GOF2 <- 1 - sum((spec_observed2-spec_predicted2)^2)/sum((spec_observed2-mean(spec_predicted2))^2)
      
      # get f_loss, clustering, and sim_nr from the filename:
      dat        <- gsub("composition", "", filename)
      dat        <- unlist(strsplit(dat, '_'))
  
      clustering <- as.numeric(dat[1])
      sim_nr     <- as.numeric(dat[2])
      f_loss     <- as.numeric(dat[3])
      
      SARs2 <- data.frame(f_loss    = f_loss,
                         clustering = clustering,
                         sim_nr     = sim_nr,
                         area       = rep(areas, 2),
                         species1   = c(spec_observed1, spec_predicted1),
                         species2   = c(spec_observed2, spec_predicted2),
                         type       = c(rep('Simulated', length(areas)),
                                  rep('METE prediction', length(areas))))
      
      GOFs2 <- data.frame(f_loss     = f_loss,
                         clustering  = clustering,
                         sim_nr      = sim_nr,
                         GOF1        = GOF1,
                         GOF2        = GOF2)
      
      GOFs <- rbind(GOFs, GOFs2)
      SARs <- rbind(SARs, SARs2)
    }
  }
} else {
  GOFs <- read.table('Fragmented-forest-communities/x64/Release/community composition/GOFs.txt',
                     header = TRUE)
  SARs <- read.table('Fragmented-forest-communities/x64/Release/community composition/SARs.txt',
                     header = TRUE)
}

# add the data for f_loss = 0 and all clustering levels (as these are all the same):
for (i in unique(GOFs$clustering[GOFs$clustering > 1]))
{
  for (j in 1:5)
  {
    sel <- (GOFs$f_loss == 0)&(GOFs$sim_nr == j)
    gofs <- data.frame(f_loss     = 0,
                       clustering = i,
                       sim_nr     = j,
                       GOF1       = GOFs$GOF1[sel],
                       GOF2       = GOFs$GOF2[sel])
  }
  GOFs <- rbind(GOFs, gofs)
}
GOFs$clustering[GOFs$clustering == 2.01] <- 2
GOFs <- GOFs[GOFs$sim_nr == 1,]

ggplot(GOFs, aes(x=f_loss*100, y=clustering, fill=GOF2)) + 
  geom_raster() + 
  scale_fill_viridis_c(name='METE goodness-of-fit', 
                       limits=c(0,1)) + 
  xlab('% Habitat loss') + 
  ylab('Clustering of habitat patches (µ)') + 
  theme_classic() + 
  theme(legend.position = 'top',
        legend.key.width=unit(1.5,"cm"))
  
ggplot(GOFs, aes(x=f_loss*100, color=factor(clustering), y=GOF2)) + 
  geom_point() + 
  geom_line() + 
  scale_color_viridis_d(name='Clustering of habitat patches (µ)') + 
  xlab('% Habitat loss') + 
  ylab('METE goodness-of-fit') + 
  theme_classic() + 
  theme(legend.position = 'top',
        legend.key.width=unit(1,"cm"))

ggplot(GOFs, aes(x=clustering, y=GOF2, color=factor(f_loss*100))) + 
  geom_point() + 
  geom_line() + 
  scale_color_viridis_d(name='% Habitat loss') + 
  xlab('Clustering of habitat patches (µ)') + 
  ylab('METE goodness-of-fit') + 
  theme_classic() + 
  theme(legend.position = 'top',
        legend.key.width=unit(1,"cm"))

SARs$clustering[SARs$clustering == 2.01] <- 2
SARs$mu <- paste('µ = ', SARs$clustering)
sel <- (SARs$type == 'Simulated')&(SARs$f_loss > 0)
ggplot(SARs[sel,], aes(x=f_loss*100, y=area, fill=species2)) + 
  geom_raster() +
  facet_wrap(vars(mu)) + 
  scale_y_continuous(trans='log10') + 
  scale_fill_viridis_c(name='No. species') + 
  xlab('% Habitat loss') + 
  ylab(expression(paste('Area size (', km^2, ')'))) + 
  theme_classic() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.key.width=unit(1.5,"cm"))
  
ggplot(SARs[sel,], aes(color=factor(f_loss*100), x=area, y=species2)) + 
  geom_point() +
  geom_line() + 
  facet_wrap(vars(mu)) + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  scale_color_viridis_d(name='% Habitat loss') + 
  ylab('No. species') + 
  xlab(expression(paste('Area size (', km^2, ')'))) + 
  theme_classic() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.key.width=unit(1.5,"cm"))

# plot the difference between simulated and METE prediction in the number of species?
SAR2 <- SARs[SARs$type == 'Simulated',]
SAR2$dif1 <- SAR2$species1 - SARs$species1[SARs$type != 'Simulated']
SAR2$dif2 <- SAR2$species2 - SARs$species2[SARs$type != 'Simulated']

sel <- (SAR2$f_loss > 0)
ggplot(SAR2[sel,], aes(x=f_loss*100, y=area, fill=dif1)) + 
  geom_raster() +
  facet_wrap(vars(mu)) + 
  scale_y_continuous(trans='log10') + 
  scale_fill_gradient2(low='red', mid='white', high='blue', midpoint=0,
    name='Δ Species') + 
  xlab('% Habitat loss') + 
  ylab(expression(paste('Area size (', km^2, ')'))) + 
  theme_classic() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.key.width=unit(1.5,"cm"))

ggplot(SAR2[SAR2$sim_nr == 1,], aes(color=factor(f_loss*100), x=area, y=dif1)) + 
  geom_point() +
  geom_line() + 
  facet_wrap(vars(mu)) + 
  scale_color_viridis_d(name='% Habitat loss') + 
  scale_x_continuous(trans='log10') + 
  ylab('Δ Species') + 
  xlab(expression(paste('Area size (', km^2, ')'))) + 
  theme_classic() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.key.width=unit(1.5,"cm"))

SAR2$loss <- paste(SAR2$f_loss*100, '% habitat loss', sep='')
ggplot(SAR2[(SAR2$sim_nr == 1)&(SAR2$f_loss > 0),], aes(color=factor(clustering), x=area, y=dif2)) + 
  geom_point() +
  geom_line() + 
  facet_wrap(vars(loss)) + 
  scale_color_viridis_d(name='µ') + 
  scale_x_continuous(trans='log10') + 
  ylab('Δ Species') + 
  xlab(expression(paste('Area size (', km^2, ')'))) + 
  theme_classic() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.key.width=unit(0.5,"cm"))


########################################################################################

sel <- GOFs$H == 6.5
p1 <- ggplot(GOFs[sel,], aes(x=n_species, y=Pmeta, fill=GOF)) + 
  geom_raster() + 
  scale_fill_viridis_c(name = 'Goodness-of-fit', 
                       direction = -1,
                       limits=c(0,1)) + 
  xlab(expression(italic(S[max]))) + 
  ylab(expression(italic(P[meta]))) + 
  #annotate(geom='text', x=9500, y=0.0025, 
  #          label=expression(paste(italic(d[max]), "= 6.5"))) + 
  theme_classic() + 
  theme(legend.position='top',
        legend.key.width=unit(2,"cm"))
p1

sel <- GOFs$Pmeta == 0.002
p2 <- ggplot(GOFs[sel,], aes(x=n_species, y=H, fill=GOF)) + 
  geom_raster() + 
  scale_fill_viridis_c(name = 'Goodness-of-fit', 
                      direction = -1,
                      limits=c(0,1)) + 
  xlab(expression(italic(S[max]))) + 
  ylab(expression(italic(d[max]))) + 
  #annotate(geom='text', x=9500, y=7.5, 
  #         label=expression(paste(italic(P[meta]), "= 0.002"))) +
  theme_classic() + 
  theme(legend.position='top',
        legend.key.width=unit(2,"cm"))
p2

sel <- GOFs$n_species == 5500
p3 <- ggplot(GOFs[sel,], aes(x=Pmeta, y=H, fill=GOF)) + 
  geom_raster() + 
  scale_fill_viridis_c(name = 'Goodness-of-fit', 
                       direction = -1,
                       limits=c(0,1)) + 
  xlab(expression(italic(P[meta]))) + 
  ylab(expression(italic(d[max]))) + 
  #annotate(geom='text', x=0.001, y=7.5, 
  #         label=expression(paste(italic(S[max]), "= 5500"))) +
  theme_classic() + 
  theme(legend.position='top',
        legend.key.width=unit(2,"cm"))
p3

ggarrange(p1, p2, p3, labels=c('A', 'B', 'C'), nrow=1,
          common.legend = TRUE)

# 9 x 4
tiff(filename = 'Fragmented-forest-communities/x64/Release/Figures/Figure S2 Calibration.tif', 
     width = 9, height = 4, units = 'in', res = 600)
ggarrange(p1, p2, p3, labels=c('A', 'B', 'C'), nrow=1,
          common.legend = TRUE)
dev.off()

# best parameter settings:
GOFs[GOFs$GOF == max(GOFs$GOF),]
# dmax = 6.5
# Pmeta = 0.002
# Smax  = 5500
# GOF = 0.948

# range GOFs: 
range(GOFs$GOF)
# GOFs range between 0.029 and 0.948

# write the GOFs-dataframe to a file:
write.table(GOFs, 
            'Fragmented-forest-communities/x64/Release/community composition calibration/GOFs.txt',
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)
