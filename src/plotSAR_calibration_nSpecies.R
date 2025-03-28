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
    f_loss     <- as.numeric(dat[3])
    dat        <- gsub(".txt", "", dat[4])
    n_species  <- as.numeric(dat)

    SAR2 <- data.frame(Area = df$Area,
                       Hab_area = df$Hab_area,
                       Species = df$Species,
                       clustering = clustering,
                       sim_nr     = sim_nr,
                       f_loss     = f_loss,
                       n_species  = n_species)
    SAR <- rbind(SAR, SAR2)
  }
}

table(SAR$f_loss, SAR$clustering)
SAR$clustering[SAR$clustering == 2.01] = 2

SAR$loss <- paste(SAR$f_loss*100, '% habitat loss')
SAR$loss <- factor(SAR$loss, levels=unique(SAR$loss))
SAR$mu   <- paste('μ = ', SAR$clustering, sep='')
SAR$mu   <- factor(SAR$mu, levels=unique(SAR$mu))

p <- ggplot(SAR[(SAR$Species > 0)&(SAR$sim_nr < 6),], aes(x=Area, y=Species, color=factor(n_species))) + 
  #geom_point() + 
  geom_line() + 
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

SARS <- data.frame(n_species = vector(length=0),
                   Hab_area = vector(length=0),
                   Species = vector(length=0),
                   type = vector(length=0)) 

maxArea = 60
for (n_species in unique(SAR$n_species))
{
  sel = (SAR$n_species == n_species)&(SAR$Hab_area <= maxArea)
  
  sar1 <- meteSAR(Amin=0.025, 
                  A0=max(SAR$Hab_area[sel])+5, 
                  S0=max(SAR$Species[sel]), 
                  N0=(max(SAR$Hab_area[sel])+5)/0.025 * 1000)
  plot(sar1$pred$A, sar1$pred$S)
  
  SAR2 <- SAR[sel,c(7, 2, 3, 8)]
  SAR2$type <- 'Observed'
  
  sar2 <- data.frame(n_species = n_species,
                     Hab_area   = sar1$pred$A,
                     Species    = sar1$pred$S,
                     #clustering = 1,
                     #sim_nr     = 1,
                     #f_loss     = 0,
                     type       = 'METE prediction')
  #sar2$loss <- paste(sar2$f_loss*100, '% habitat loss')
  #sar2$mu   <- paste('μ = ', sar2$clustering, sep='')
  
  SAR2 <- rbind(SAR2, sar2)
  
  ggplot(SAR2, aes(x=Hab_area, y=Species, color=type)) + 
    geom_line()  +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10')
  
  SARS <- rbind(SARS, SAR2)
}


ggplot(SARS, aes(x=Hab_area, y=Species, color=type)) +
  geom_point()  +
  facet_wrap(vars(n_species)) + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +
  xlab('Area size') + 
  ylab('Number of species') + 
  theme_bw()


filenames
filename <- filenames[1]

print(filename)
filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
df <- read.table(filename2)

sar1 <- meteSAR(spp = df$V3,
                abund = df$V4,
                row = df$V1+1,
                col = df$V2+1,
                Amin=1, 
                A0=2025)

plot(sar1)
warnings()

# anders moet ik het zelf maar doen, want zo werkt het niet... 
# beginnen bij het hele grid, en dan steeds halveren: 

dat <- expand.grid(x=0:44, y=0:44)
uPatchID <- paste(dat$x, dat$y, sep='-')

codes <- matrix('A', nrow = 2025, ncol= 11)

x <- dat$x
y <- dat$y

i <- 2

while (max(x) >= 2){
  # divide into two equal parts: 
  codes2 <- codes[,1]
  codes2[x < max(x)/2] <- 'B'
  codes[,i] <- codes2
  i <- i + 1
  
  # again: 
  codes2 <- codes[,1]
  codes2[y <= max(y)/2] <- 'B'
  codes[,i] <- codes2
  i <- i + 1
  
  # dubbelklappen: 
  x[x >= max(x)/2] <- x[x >= max(x)/2] - max(x)/2
  y[y >= max(y)/2] <- y[y >= max(y)/2] - max(y)/2
  
  max(x)
  max(y)
  i
}

# als het goed is, is nu iedere rij uniek:
codes2 <- paste(codes[,1], codes[,2], codes[,3], codes[,4], codes[,5],codes[,6],
                codes[,7], codes[,8], codes[,9], codes[,10], codes[,11])
length(unique(codes2))

# write the codes to a file, so we can use them: 
write.table(codes, 'results/SARgroupcodes.txt', 
            row.names = FALSE, col.names = FALSE, append = FALSE)


filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition/')

for (j in 1:length(filenames))
{
  
  filename <- filenames[j]
  
  dat <- expand.grid(x=0:44, y=0:44)
  uPatchID <- paste(dat$x, dat$y, sep='-')
  
  if (grepl("composition", filename))
  {
    
    filename_out = gsub("composition", "SAR2_", filename)
    if ((filename_out %in% filenames) == FALSE)
    {
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
      df <- read.table(filename2)
      
      # dan nu de kolommen afgaan om data te groeperen:
      patchID <- paste(df$V1, df$V2, sep='-')
      
      SAR <- data.frame(Area = length(uPatchID),
                      Hab_area = length(uPatchID),
                      Species = length(unique(df$V3)))
      
      codes2 <- codes[,1]
      for (i in 2:ncol(codes)){
        codes2 <- paste(codes2, codes[,i])
        uCodes <- unique(codes2)
        
        for (j2 in uCodes){
          SAR2 <- data.frame(Area = length(uPatchID[codes2 == j2]),
                            Hab_area = length(unique(patchID[patchID %in% uPatchID[codes2 == j2]])),
                            Species = length(unique(df$V3[patchID %in% uPatchID[codes2 == j2]])))
          SAR <- rbind(SAR, SAR2)
        }
      }
      
      ggplot(SAR, aes(x=Hab_area, y=Species)) + 
      geom_point()
      
      # write to file: 
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename_out, sep='')
      write.table(SAR, filename2, 
                row.names = FALSE, 
                col.names = TRUE, 
                append = FALSE)
    }
  }
}




# compare the SARs: 
# load in data from the SAR files:
filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition/')

SAR <- data.frame(Area = vector(length=0),
                  Hab_area = vector(length=0),
                  Species = vector(length=0),
                  clustering = vector(length=0),
                  sim_nr     = vector(length=0),
                  f_loss     = vector(length=0),
                  n_species  = vector(length=0),
                  type       = vector(length=0))

for (j in 1:length(filenames))
{
  
  filename <- filenames[j]
  
  if (grepl("SAR", filename))
  {
    print(filename)
    
    filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
    df <- read.table(filename2, header = TRUE)

    if (grepl("SAR2_", filename))
    {
      dat        <- gsub("SAR2_", "", filename)
      type       <- 'SAR2'
    } else { 
      dat        <- gsub("SAR", "", filename)
      type       <- 'SAR1'
    }
    
    dat        <- unlist(strsplit(dat, '_'))
    clustering <- as.numeric(dat[1])
    sim_nr     <- as.numeric(dat[2])
    f_loss     <- as.numeric(dat[3])
    dat        <- gsub(".txt", "", dat[4])
    n_species  <- as.numeric(dat)
    
    SAR2 <- data.frame(Area = df$Area,
                       Hab_area = df$Hab_area,
                       Species = df$Species,
                       clustering = clustering,
                       sim_nr     = sim_nr,
                       f_loss     = f_loss,
                       n_species  = n_species,
                       type       = type)
    SAR <- rbind(SAR, SAR2)
  }
}

SAR$Hab_area[SAR$type == 'SAR2'] <- SAR$Hab_area[SAR$type == 'SAR2'] *0.025

ggplot(SAR, aes(x=Hab_area, y=Species, color=type)) + 
  geom_point() + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  facet_wrap(vars(n_species)) + 
  theme_bw()
