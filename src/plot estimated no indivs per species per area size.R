# Here, we examine the relation between the assumed n = n0 * A/A0, that is used 
# in METE to estimate the SAR, and the observed number of individuals n per 
# area size A. 

library(ggplot2)

##### Step 1:
# First, examine this relation in a continuous environment:
filename  <- 'composition1.00_1.00_0.00_5500.00_2.00_6.50.txt'

plot_relation <- function(filename)
{
  filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
  
  df        <- read.table(filename2)
  names(df) <- c('x', 'y', 'species', 'n')
  df        <- df[df$x != -1,]
  
  # per species, we need to know n0 (the number of individuals in A0):
  n0 <- tapply(df$n, as.factor(df$species), sum)
  species <- tapply(df$species, as.factor(df$species), mean)
  
  df$n0 <- 0
  for (spec in species){
    df$n0[df$species == spec] <- n0[species == spec]
  }
  
  patchID  <- unique(paste(df$x, df$y, sep='-'))
  df$patch <- paste(df$x, df$y, sep='-')
  A0       <- length(patchID)
  
  # create a dataframe where we increment the area size: 
  df2 <- data.frame(A  = 1,
                    n  = df$n[df$patch == patchID[1]],
                    n0 = df$n0[df$patch == patchID[1]])
  areas <- c(0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6)/0.025
  for (i in areas){
    sel  <- df$patch %in% patchID[1:i] 
    n    <- tapply(df$n[sel], as.factor(df$species[sel]), sum)
    spec <- tapply(df$species[sel], as.factor(df$species[sel]), mean)
    n0   <- tapply(df$n0[sel], as.factor(df$species[sel]), mean)
    # hier verder met het berekenen van n0 per species
    df2a <- data.frame(A = i,
                       n = n,
                       n0 = n0)
    df2 <- rbind(df2, df2a)
  }
  
  # calculate a goodness of fit:
  GOF1 <- 1 - sum((df2$n-(df2$n0*df2$A/2025))^2)/sum((df2$n-mean(df2$n0*df2$A/2025))^2)
  
  return(GOF1)
}
GOF <- plot_relation(filename)


##### Step 2: 
# Then, examine this relation in environments with different levels of 
# habitat loss and habitat fragmentation: 
filename <- 'composition5.00_1.00_0.70_5500.00_2.00_6.50.txt'

df2 <- plot_relation(filename)


plot_relation2 <- function(filename)
{
  patchID1 <- paste(sort(rep(0:44, 45)), rep(0:44, 45), sep='-')
  codes    <- read.table('results/SARgroupcodes.txt') 
  
  filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
  
  df        <- read.table(filename2)
  names(df) <- c('x', 'y', 'species', 'n')
  df        <- df[df$x != -1,]
  
  # per species, we need to know n0 (the number of individuals in A0):
  n0 <- tapply(df$n, as.factor(df$species), sum)
  species <- tapply(df$species, as.factor(df$species), mean)
  
  # let's continue with only the most abundant species:
  #df <- df[df$species == species[n0 == max(n0)],]
  #df$n0 <- max(n0)  
  
  patchID  <- unique(paste(df$x, df$y, sep='-'))
  df$patch <- paste(df$x, df$y, sep='-')

  # per species per patch:
  
  # hier moet ik verder mee, maar heb nu niet de mentale capaciteit :(
  # wat ik wil laten zien is de relatie tussen 1. a/A (hoeveelheid habitat 
  # per hoeveelheid gebied) en 2. n / (n0 * A/A0) (oftewel: het werkelijke 
  # aantal individuen per geschatte aantal individuen). En dit alles per soort, 
  # en voor 12 verschillende groepsgroottes (van 1 patch tot de helft van het
  # totale gebied bij elkaar). Dit levert een hele hoop datapunten op per 
  # simulatie, en dus nog veel meer als we alle simulaties bij elkaar zetten. 
  # Daar moet ik dus ook nog een oplossing voor vinden. Maar ja, mijn hoofd doet
  # het vandaag even niet... Veel succes hiermee, mezelf in de toekomst :)
  
  
  a <- A0 / 2025   # average habitable area size
  A <- 1
  n <- tapply(df$n[sel], as.factor(df$species[sel]), sum)
  n <- n / 2025
  n0<- tapply(df$n0[sel], as.factor(df$species[sel]), mean)
  
  
  
  
  # create a dataframe where we increment the area size: 
  df2 <- data.frame(a  = 1,
                    A  = (1:2025)[patchID1 == patchID[1]],
                    n  = df$n[df$patch == patchID[1]],
                    n0 = df$n0[df$patch == patchID[1]])
  areas <- c(0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6)/0.025
  for (i in areas){
    sel  <- df$patch %in% patchID[1:i] 
    n    <- tapply(df$n[sel], as.factor(df$species[sel]), sum)
    spec <- tapply(df$species[sel], as.factor(df$species[sel]), mean)
    n0   <- tapply(df$n0[sel], as.factor(df$species[sel]), mean)

    df2a <- data.frame(a = i,
                       A  = (1:2025)[patchID1 == patchID[i]],
                       n = n,
                       n0 = n0)
    df2 <- rbind(df2, df2a)
  }
  
  p1 <- ggplot(df2, aes(x=n0*A/A0, y=n)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + 
    geom_smooth(color='red') + 
    scale_x_continuous(trans='log10',
                       limits=c(1, max(c(df2$n, df2$n0*df2$A/A0)))*2) + 
    scale_y_continuous(trans='log10',
                       limits=c(1, max(c(df2$n, df2$n0*df2$A/A0)))*2) +  
    theme_bw()
  
  print(p1)
  return(df2)
}

filename <- 'composition2.01_1.00_0.70_5500.00_2.00_6.50.txt'
df2 <- plot_relation2(filename)

# for all parameter combinations:
filenames <- list.files('Fragmented-forest-communities/x64/Release/community composition/')

df <- data.frame(clustering = vector(length=0),
                 sim_nr     = vector(length=0),
                 f_loss     = vector(length=0),
                 A          = vector(length=0),
                 n          = vector(length=0),
                 n0         = vector(length=0))

df_gof <- data.frame(clustering = vector(length=0),
                     sim_nr     = vector(length=0),
                     f_loss     = vector(length=0),
                     GOF        = vector(length=0))

for (j in 1:length(filenames))
{
  filename <- filenames[j]
  
  if (grepl("composition", filename))
  {
    # get f_loss, clustering, and sim_nr from the filename:
    dat        <- gsub("composition", "", filename)
    dat        <- unlist(strsplit(dat, '_'))
    
    clustering <- as.numeric(dat[1])
    sim_nr     <- as.numeric(dat[2])
    f_loss     <- as.numeric(dat[3])
    
    group2 <- paste(clustering, sim_nr, f_loss, sep='-')
    if (!(group2 %in% group)){
      df2 <- plot_relation2(filename)
      
      df3 <- data.frame(clustering = clustering,
                        sim_nr     = sim_nr,
                        f_loss     = f_loss,
                        A          = df2$A,
                        n          = df2$n,
                        n0         = df2$n0)
      df <- rbind(df, df3)
      
      GOF     <- plot_relation(filename)
      df_gof2 <- data.frame(clustering = clustering,
                            sim_nr     = sim_nr,
                            f_loss     = f_loss,
                            GOF        = GOF)
      df_gof <- rbind(df_gof, df_gof2)
    }
  }
}

group <- paste(df$clustering, df$sim_nr, df$f_loss, sep='-')



A0 <- 2025
df$clustering[df$clustering == 2.01] <- 2
sel <- (df$f_loss %in% c(0.1, 0.3, 0.5, 0.7, 0.9))&
  (df$clustering %in% c(1, 2, 3, 4, 5))

df$mu   <- paste('µ = ', df$clustering, sep='')
df$loss <- paste(df$f_loss*100, '% Habitat loss', sep='')
ggplot(df[sel,], aes(x=n0*A/A0, y=n, color=as.factor(sim_nr))) + 
  geom_point(alpha=0.3) +
  geom_smooth(se=FALSE) + 
  facet_grid(cols=vars(mu), rows=vars(loss)) + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') +  
  geom_abline(slope=1, intercept=0) + 
  theme_bw() + 
  theme(legend.position='none')
  
df_gof$clustering[df_gof$clustering == 2.01] <- 2
ggplot(df_gof, aes(x=clustering, y=f_loss, fill=A)) + 
  geom_raster() + 
  facet_wrap(vars(sim_nr)) + 
  scale_fill_viridis_c(name='Goodness-of-fit') +
  theme_bw() + 
  theme(legend.position='top')

# calculate a mean GOF for all simulation runs per clustering and f_loss:
group <- paste(df_gof$clustering, df_gof$f_loss, sep='-')
df2 <- data.frame(clustering = tapply(df_gof$clustering, group, mean),
                  f_loss = tapply(df_gof$f_loss, group, mean),
                  GOF = tapply(df_gof$A, group, mean))
# add the data for f_loss = 1:
uCl  <- unique(df2$clustering)
df2a <- df2[(df2$clustering == 1)&(df2$f_loss == 0),] 
for (i in uCl[2:9])
{
  df2a$clustering <- i
  df2 <- rbind(df2, df2a)
}


p1 <- ggplot(df2, aes(x=clustering, y=f_loss*100, fill=GOF)) + 
  geom_raster() + 
  scale_fill_viridis_c(name='Goodness-of-fit', limits=c(0,1)) +
  xlab('Clustering of habitat patches (µ)') + 
  ylab('% Habitat loss') + 
  theme_bw() + 
  theme(legend.position='top',
        legend.key.width = unit(2, "cm"))
  
# 6 x 5:
tiff(filename = 'Fragmented-forest-communities/x64/Release/Figures/Figure gof estimated number of individuals per species per area size.tif', 
     width = 6, height = 5, units = 'in', res = 600)
p1
dev.off()

# write the data to a file so you don't have to run this again:
#write.table(df2, 'gof estimated number of individuals per species per area size.txt',
#            row.names = FALSE, col.names = FALSE, append = FALSE)



