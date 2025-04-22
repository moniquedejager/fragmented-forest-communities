library(ggplot2)

# read in the community composition data and calculate the species area curve. 
# write these data to a new file. 


filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition/')

for (j in 1:length(filenames))
{
  
  filename <- filenames[j]
  
  dat <- expand.grid(x=0:44, y=0:44)
  uPatchID <- paste(dat$x, dat$y, sep='-')
  
  if (grepl("composition", filename))
  {
    
    filename_out = gsub("composition", "SAR", filename)
    if ((filename_out %in% filenames) == FALSE)
    {
      SAR <- data.frame(Area = vector(length=0),
                        Hab_area = vector(length=0),
                        Species = vector(length=0))
      
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
      df <- read.table(filename2)
      
      names(df) <- c('x', 'y', 'species', 'n')
      df <- df[df$x != -1,]
      
      patchID <- paste(df$x, df$y, sep='-')
      #nspecies <- tapply(df$species, patchID, length)
      
      # median dispersal distance:
      #df$median_dispersal <- log(0.5)/(-1 * df$species/1000)
      
      uPatch2 <- unique(patchID)
      for (i in 1:length(uPatchID))
      {
        nspecies <- length(unique(df$species[patchID %in% uPatchID[1:i]]))
        SAR2 <- data.frame(Area = i*0.025, #patch size = 2.5 ha
                           Hab_area = length(uPatch2[uPatch2 %in% uPatchID[1:i]]) * 0.025,
                           Species = nspecies)
        SAR <- rbind(SAR, SAR2)
      }
      
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename_out, sep='')
      write.table(SAR, filename2, 
                  row.names = FALSE, 
                  col.names = TRUE, 
                  append = FALSE)
    }
  }
}


# load in data from the SAR files:
filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition/')

SAR <- data.frame(Area = vector(length=0),
                  Hab_area = vector(length=0),
                  Species = vector(length=0),
                  clustering = vector(length=0),
                  sim_nr     = vector(length=0),
                  f_loss     = vector(length=0))

for (j in 1:length(filenames))
{
  
  filename <- filenames[j]
  
  if (grepl("SAR", filename))
  {
    print(filename)
    filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
    df <- read.table(filename2, header = TRUE)
    
    dat        <- gsub("SAR", "", filename)
    dat        <- unlist(strsplit(dat, '_'))
    clustering <- as.numeric(dat[1])
    sim_nr     <- as.numeric(dat[2])
    dat        <- gsub(".txt", "", dat[3])
    f_loss     <- as.numeric(dat)
    
    SAR2 <- data.frame(Area = df$Area,
                       Hab_area = df$Hab_area,
                       Species = df$Species,
                       clustering = clustering,
                       sim_nr     = sim_nr,
                       f_loss     = f_loss)
    SAR <- rbind(SAR, SAR2)
  }
}

table(SAR$f_loss, SAR$clustering)
SAR$clustering[SAR$clustering == 2.01] = 2

SAR$loss <- paste(SAR$f_loss*100, '% habitat loss')
SAR$loss <- factor(SAR$loss, levels=unique(SAR$loss))
SAR$mu   <- paste('μ = ', SAR$clustering, sep='')
SAR$mu   <- factor(SAR$mu, levels=unique(SAR$mu))

p <- ggplot(SAR[(SAR$Species > 0)&(SAR$sim_nr < 6),], aes(x=Area, y=Species, color=factor(sim_nr))) + 
  #geom_point() + 
  geom_line() + 
  facet_grid(cols=vars(loss), rows=vars(mu)) + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  scale_color_viridis_d(direction=-1) + 
  xlab('Area size (km2)') + 
  ylab('Number of species') + 
  theme_bw() + 
  theme(legend.position = 'none', 
        strip.placement = "outside", 
        strip.background = element_blank())
p

# 12 x 8:
tiff(filename = 'Fragmented-forest-communities/x64/Release/Figures/Figure 1.tif', 
     width = 12, height = 8, units = 'in', res = 600)
p
dev.off()


p <- ggplot(SAR[(SAR$Species > 0)&(SAR$sim_nr < 6),], aes(x=Hab_area, y=Species, color=factor(sim_nr))) + 
  geom_line() + 
  facet_grid(cols=vars(loss), rows=vars(mu), scales='free_x') + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  scale_color_viridis_d(direction=-1) + 
  xlab('Habitable area size (km2)') + 
  ylab('Number of species') + 
  theme_bw() + 
  theme(legend.position = 'none', 
        strip.placement = "outside", 
        strip.background = element_blank())
p

# 12 x 8:
tiff(filename = 'Fragmented-forest-communities/x64/Release/Figures/Figure 2.tif', 
     width = 12, height = 8, units = 'in', res = 600)
p
dev.off()

p <- ggplot(SAR[(SAR$Species > 0)&(SAR$sim_nr == 5),], aes(x=Hab_area, y=Species, color=loss)) + 
  geom_line() + 
  facet_wrap(vars(mu)) + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  scale_color_viridis_d(direction=-1, name='') + 
  xlab('Habitable area size (km2)') + 
  ylab('Number of species') + 
  theme_bw() + 
  theme(legend.position = 'top', 
        strip.placement = "outside", 
        strip.background = element_blank())
p

# now, let's take a look at the METE predictions of the species-area relation:
library(meteR)

loss       = 0
sim_nr     = 1
clustering = 1
sel = (SAR$f_loss == loss)&(SAR$clustering == clustering)&(SAR$sim_nr == sim_nr)

sar1 <- meteSAR(Amin=0.025, 
                A0=50, 
                S0=max(SAR$Species[(sel)&(SAR$Hab_area <= 51)]), 
                N0=50/0.025 * 1000)
plot(sar1$pred$A, sar1$pred$S)

SAR2 <- SAR[sel,]
SAR2$type <- 'Observed'

sar2 <- data.frame(Area       = NA,
                   Hab_area   = sar1$pred$A,
                   Species    = sar1$pred$S,
                   clustering = clustering,
                   sim_nr     = sim_nr,
                   f_loss     = loss,
                   type       = 'METE prediction')
sar2$loss <- paste(sar2$f_loss*100, '% habitat loss')
sar2$mu   <- paste('μ = ', sar2$clustering, sep='')

SAR2 <- rbind(SAR2, sar2)

ggplot(SAR2[SAR2$Hab_area <= 51,], aes(x=Hab_area, y=Species, color=type)) + 
  geom_line()  
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10')

