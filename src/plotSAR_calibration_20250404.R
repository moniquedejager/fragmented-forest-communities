library(ggplot2)
library(ggpubr)

# read in the community composition data and calculate part of the species 
# area curve. Fit the simulation results to the empirical data presented in 
# Ter Steege et al. 2020, to calibrate our model. 

GOFs <- data.frame(mut_rate  = vector(length=0),
                   n_species = vector(length=0),
                   H         = vector(length=0),
                   GOF       = vector(length=0)) 

dat            <- expand.grid(x=0:44, y=0:44)
uPatchID       <- paste(dat$x, dat$y, sep='-')
areas          <- c(10.725, 19.45)
spec_predicted <- c(4531, 5027)

filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition calibration/')

if (!("GOFs.txt" %in% filenames)){
  for (j in 1:length(filenames))
  {
    
    filename <- filenames[j]
    
    if (grepl("composition", filename))
    {
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition calibration/', filename, sep='')
      df <- read.table(filename2)
      
      names(df) <- c('x', 'y', 'species', 'n')
      df <- df[df$x != -1,]
      
      patchID <- paste(df$x, df$y, sep='-')
      uPatch2 <- unique(patchID)
      
      # For the two area sizes, find how many species are in the simulated
      # landscape:
      spec_observed <- vector(length=0)
      for (i in areas/0.025)
      {
        nspecies      <- length(unique(df$species[patchID %in% uPatchID[1:i]]))
        spec_observed <- c(spec_observed, nspecies)
      }
      
      # goodness-of-fit:
      GOF <- 1 - sum((spec_observed-spec_predicted)^2)/sum((spec_observed-mean(spec_predicted))^2)
      
      # get Smax, dmax, and Pmeta from the filename:
      dat        <- gsub("composition", "", filename)
      dat        <- unlist(strsplit(dat, '_'))
  
      clustering <- as.numeric(dat[1])
      sim_nr     <- as.numeric(dat[2])
      f_loss     <- as.numeric(dat[3])
      n_species  <- as.numeric(dat[4])
      mut_rate   <- as.numeric(dat[5])
      dat        <- gsub(".txt", "", dat[6])
      H          <- as.numeric(dat)
  
      GOFs2 <- data.frame(mut_rate = mut_rate,
                         n_species = n_species,
                         H         = H,
                         GOF       = GOF)
      GOFs <- rbind(GOFs, GOFs2)
    }
  }
} else {
  GOFs <- read.table('Fragmented-forest-communities/x64/Release/community composition calibration/GOFs.txt',
                     header = TRUE)
}

GOFs$Pmeta <- GOFs$mut_rate / 1000

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
