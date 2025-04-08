library(ggplot2)

# read in the community composition data and calculate the species area curve. 
# write these data to a new file. 


filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition_calibration/')

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
filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition_calibration/')

SAR <- data.frame(Area = vector(length=0),
                  Hab_area = vector(length=0),
                  Species = vector(length=0),
                  clustering = vector(length=0),
                  sim_nr     = vector(length=0),
                  f_loss     = vector(length=0),
                  n_species  = vector(length=0),
                  mut_rate   = vector(length=0),
                  H          = vector(length=0))

for (j in 1:length(filenames))
{
  
  filename <- filenames[j]
  
  if (grepl("SAR", filename))
  {
    print(filename)
    filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition_calibration/', filename, sep='')
    df <- read.table(filename2, header = TRUE)
    
    dat        <- gsub("SAR", "", filename)
    dat        <- unlist(strsplit(dat, '_'))
    print(length(dat))
    
    clustering <- as.numeric(dat[1])
    sim_nr     <- as.numeric(dat[2])
    f_loss     <- as.numeric(dat[3])
    
    if (length(dat) == 4){
      dat        <- gsub(".txt", "", dat[4])
      n_species  <- as.numeric(dat)
      mut_rate   <- 0.0003*1000
      H          <- 5.5
    } else {
      n_species  <- as.numeric(dat[4])
      if (length(dat) == 5){
        dat        <- gsub(".txt", "", dat[5])
        mut_rate   <- as.numeric(dat)
        H          <- 5.5
      } else {
        mut_rate   <- as.numeric(dat[5])
        dat        <- gsub(".txt", "", dat[6])
        H          <- as.numeric(dat)
      }
    }

    SAR2 <- data.frame(Area = df$Area,
                       Hab_area = df$Hab_area,
                       Species = df$Species,
                       clustering = clustering,
                       sim_nr     = sim_nr,
                       f_loss     = f_loss,
                       n_species  = n_species,
                       mut_rate   = mut_rate,
                       H          = H)
    SAR <- rbind(SAR, SAR2)
  }
}

table(SAR$f_loss, SAR$clustering)
SAR$clustering[SAR$clustering == 2.01] = 2
SAR$mut_rate <- SAR$mut_rate / 1000

SAR$loss <- paste(SAR$f_loss*100, '% habitat loss')
SAR$loss <- factor(SAR$loss, levels=unique(SAR$loss))
SAR$mu   <- paste('μ = ', SAR$clustering, sep='')
SAR$mu   <- factor(SAR$mu, levels=unique(SAR$mu))

sel <- (SAR$Species > 0) #&(SAR$mut_rate < 0.01)
p <- ggplot(SAR[sel,], aes(x=Area, y=Species, color=factor(H))) + 
  #geom_point() + 
  geom_line() + 
  facet_grid(cols=vars(mut_rate), rows=vars(n_species), scales='free_y') +
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  ylim(c(0, 7000)) + 
  scale_color_viridis_d(direction=-1) + 
  xlab('Area size (km2)') + 
  ylab('Number of species') + 
  theme_bw() + 
  theme(legend.position = 'top', 
        strip.placement = "outside", 
        strip.background = element_blank())
p

# now, let's take a look at the METE predictions of the species-area relation:
library(meteR)

SARS <- data.frame(mut_rate = vector(length=0),
                   Hab_area = vector(length=0),
                   Species  = vector(length=0),
                   type     = vector(length=0),
                   H        = vector(length=0)) 

GOFs <- data.frame(mut_rate  = vector(length=0),
                   n_species = vector(length=0),
                   H         = vector(length=0),
                   GOF       = vector(length=0),
                   GOF_log   = vector(length=0)) 
  
maxArea = 60
for (mut_rate in unique(SAR$mut_rate))
{
  for (H in unique(SAR$H))
  {
    for (n_species in unique(SAR$n_species))
    {
      sel = (SAR$mut_rate == mut_rate)&(SAR$Hab_area <= maxArea)&(SAR$H == H)&(SAR$n_species == n_species)
      
      if (nrow(SAR[sel,]) > 0)
      {
        sar1 <- meteSAR(Amin=0.025, 
                        A0=max(SAR$Hab_area[sel]), 
                        S0=max(SAR$Species[sel]), 
                        N0=(max(SAR$Hab_area[sel]))/0.025 * 1000)
        plot(sar1$pred$A, sar1$pred$S)
        
        SAR2 <- SAR[sel,c(8, 2, 3, 7, 9)]
        SAR2$type <- 'Observed'
        
        sar2 <- data.frame(mut_rate = mut_rate,
                           Hab_area   = sar1$pred$A,
                           Species    = sar1$pred$S,
                           n_species  = n_species,
                           #clustering = 1,
                           #sim_nr     = 1,
                           #f_loss     = 0,
                           type       = 'METE prediction',
                           H          = H)
        #sar2$loss <- paste(sar2$f_loss*100, '% habitat loss')
        #sar2$mu   <- paste('μ = ', sar2$clustering, sep='')
        
        spec_observed  <- SAR2$Species[round(SAR2$Hab_area, 2) %in% round(sar2$Hab_area, 2)]
        spec_predicted <- sar2$Species 
        
        # goodness-of-fit:
        GOF <- 1 - sum((spec_observed-spec_predicted)^2)/sum((spec_observed-mean(spec_predicted))^2)
        
        # goodness-of-fit on a log scale:
        spec_observed  <- log(spec_observed)
        spec_predicted <- log(spec_predicted)
        GOF_log <- 1 - sum((spec_observed-spec_predicted)^2)/sum((spec_observed-mean(spec_predicted))^2)
        
        SAR2 <- rbind(SAR2, sar2)
        
        ggplot(SAR2, aes(x=Hab_area, y=Species, color=type)) + 
          geom_line()  +
          scale_x_continuous(trans='log10') + 
          scale_y_continuous(trans='log10')
        
        SARS <- rbind(SARS, SAR2)
        
        GOFs2 <- data.frame(mut_rate  = mut_rate,
                            n_species = n_species,
                            H         = H,
                            GOF       = GOF,
                            GOF_log   = GOF_log)
        
        GOFs <- rbind(GOFs, GOFs2)
      }
    }    
  }
}

# and how close are they to the 2019 survey, which held 5027 tree species in 
# 19.46 km2 (Ter Steege et al., 2020)
GOFS$dif <- 0
for (mut_rate in unique(SAR$mut_rate))
{
  for (H in unique(SAR$H))
  {
    for (n_species in unique(SAR$n_species))
    {
      sel = (SAR$mut_rate == mut_rate)&
        (SAR$Hab_area <= maxArea)&(SAR$H == H)&
        (SAR$n_species == n_species)
      if (nrow(SAR[sel,]) > 0)
      {
        SAR2 <- SAR[sel,c(8, 2, 3, 7, 9)]
        areas <- c(10.725, 19.45)
        
        spec_observed  <- SAR2$Species[SAR2$Hab_area %in% areas]
        spec_predicted <- c(4531, 5027)
        
        # goodness-of-fit:
        GOF <- 1 - sum((spec_observed-spec_predicted)^2)/sum((spec_observed-mean(spec_predicted))^2)
        
        GOFs$dif[(GOFs$mut_rate == mut_rate)&
                   (GOFs$n_species == n_species)&
                   (GOFs$H == H)] <- GOF
      }
    }    
  }
}
#GOFs$dif  <- abs(GOFs$dif)
#GOFs$best <- 1/GOFs$dif
#GOFs$best  <- GOFs$best / max(GOFs$best)

ggplot(SARS, aes(x=Hab_area, y=Species, color=type)) +
  geom_point()  +
  facet_grid(cols=vars(mut_rate), rows=vars(n_species), scales='free_y') + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +
  scale_color_discrete(name='') + 
  xlab('Area size') + 
  ylab('Number of species') + 
  #ylim(c(0, 6600)) + 
  theme_bw() + 
  theme(legend.position = 'top')

GOFs$fit <- GOFs$dif
GOFs$fit[GOFs$fit < 0] <- 0

sel <- (GOFs$H == 5.5) & 
  (GOFs$mut_rate %in% unique(GOFs$mut_rate[GOFs$n_species == 10000]))
sel <- GOFs$fit >= 0.4 
p1 <- ggplot(GOFs[sel,], aes(x=n_species, y=fit, color=as.factor(mut_rate))) + 
  geom_line(linewidth=1.1) + 
  geom_point(size=2) + 
  xlab(expression(italic(S[max]))) + 
  ylab('GoF to empirical data') + 
  facet_wrap(vars(H)) + 
  scale_color_discrete(name=expression(italic(P[meta]))) + 
  #scale_y_continuous(trans='log10', limits = c(0.0001,1)) + 
  #scale_x_continuous(trans='log10') + 
  theme_bw()  
  #theme(legend.position = c(0.8, 0.8))
p1

sel <- (GOFs$n_species == 6000)
ggplot(GOFs[sel,], aes(x=mut_rate, y=best, color=as.factor(H))) + 
  geom_line(linewidth=1.1) + 
  geom_point(size=2) + 
  xlab(expression(italic(P[meta]))) + 
  ylab('Relative fit to empirical data') + 
  scale_color_discrete(name=expression(italic(d[max]))) + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') + 
  theme_bw() + 
  theme(legend.position = c(0.8, 0.8))


# waardes uit de SAR van Hans: 
Species = c(4587.126966,
            4853.980899,
            5002.857303,
            5154.542697,
            5320.273034,
            5469.149438,
            5609.598876,
            5713.531461,
            5828.7,
            5915.778652,
            6002.857303,
            7000,
            15874)
Area = c(1007.840314,
         1308.887435,
         1557.578534,
         1845.536649,
         2277.473822,
         2643.965969,
         3206.793194,
         3481.662304,
         4083.756545,
         4515.693717,
         4973.808901,
         16000,
         1946/0.00035)
Area <- Area / 100

plot(Area, Species, xlim=c(1, max(Area)), 
     ylim=c(1, 16000), col='red', pch=16, type='b', log='xy')

sel <- Area < 16000
sar1 <- meteSAR(Amin=0.025, 
                A0=max(Area[sel]), 
                S0=max(Species[sel]), 
                N0=max(Area[sel]) * 100* 400)
plot(Area[sel], Species[sel], xlim=c(0.025, 50), 
     ylim=c(1, max(Species[sel])), col='red', type='l', lwd=2) #, log='xy')
lines(sar1$pred$A, sar1$pred$S, col='blue', lwd=2, lty=2)




sel <- (SARS$n_species == 50000)
df <- data.frame(n_species = 50000,
                 Hab_area = Area/100,
                 Species = Species, 
                 type = 'Estimation from Ter Steege et al. (2020)')
df <- rbind(df, SARS[sel,])

ggplot(df, aes(x=Hab_area, y=Species, color=type)) + geom_point()

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
                  mut_rate   = vector(length=0),
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
    n_species  <- as.numeric(dat[4])
    dat        <- gsub(".txt", "", dat[5])
    mut_rate   <- as.numeric(dat)
    
    SAR2 <- data.frame(Area = df$Area,
                       Hab_area = df$Hab_area,
                       Species = df$Species,
                       clustering = clustering,
                       sim_nr     = sim_nr,
                       f_loss     = f_loss,
                       n_species  = n_species,
                       mut_rate   = mut_rate,
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
