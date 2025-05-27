# wat ik nu wil doen, is kijken naar de relatie tussen 
# 1. het verschil in het aantal soorten in een gebied van grootte A, tussen 
# de oude situatie (alles is habitat) en de nieuwe (met habitatsverlies en 
# fragmentatie) en 
# 2. de relatie tussen habitatgebied en gebied (a/A). 

library(ggplot2)
library(ggpubr)

# first, write the following data per continuous, unfragmented habitat to a file:
# sim_nr
# S per A (# species per area size)

df2 <- data.frame(simnr = 0,
                  A = 0,
                  S = 0)
for (simnr in 1:5){
  filename  <- paste('composition1.00_', simnr, '.00_0.00_5500.00_2.00_6.50.txt', sep='')
  filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
  
  df        <- read.table(filename2)
  names(df) <- c('x', 'y', 'species', 'n')
  df        <- df[df$x != -1,]
  
  patchID  <- unique(paste(df$x, df$y, sep='-'))
  df$patch <- paste(df$x, df$y, sep='-')
  A0       <- length(patchID)
  
  for (i in 1:2025){
    A <- i
    S <- length(unique(df$species[df$patch %in% patchID[1:i]]))
    df2 <- rbind(df2, list(simnr, A, S)) 
  }
}
df2 <- df2[df2$simnr > 0,]

ggplot(df2, aes(x=A, y=S, color=factor(simnr))) + geom_point()

filename3 <- 'Fragmented-forest-communities/x64/Release/results/SARinPristineLandscape.txt'
write.table(df2, filename3, row.names = FALSE, col.names = TRUE, append = FALSE)

# now, we can use the data from the pristine environment to calculate the 
# biodiversity loss (in terms of numbers of species) per area size, for all
# simulations:

df3 <- data.frame(simnr      = 0,
                  clustering = 0,
                  f_loss     = 0,
                  A          = 0,
                  S          = 0,
                  a          = 0,
                  s          = 0)

for (clust in c('1.00', '2.01', '5.00')){
  for (loss in c(0.3, 0.5, 0.7, 0.9)){
    for (simnr in 1:5){
      filename <- paste('composition', clust,'_', simnr, 
                        '.00_', loss, '0_5500.00_2.00_6.50.txt', sep='')
      
      #filename <- 'composition1.00_1.00_0.10_5500.00_2.00_6.50.txt'
      
      # derive the simulation number, clustering, and habitat loss from the filename:
      dat        <- gsub("composition", "", filename)
      dat        <- unlist(strsplit(dat, '_'))
      
      clustering <- as.numeric(dat[1])
      sim_nr     <- as.numeric(dat[2])
      f_loss     <- as.numeric(dat[3])
      
      df3a <- data.frame(simnr      = sim_nr,
                         clustering = clustering,
                         f_loss     = f_loss,
                         A          = df2$A[df2$simnr == sim_nr],
                         S          = df2$S[df2$simnr == sim_nr])
      
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
      df        <- read.table(filename2)
      names(df) <- c('x', 'y', 'species', 'n')
      df        <- df[df$x != -1,]
      
      df$patch <- paste(df$x, df$y, sep='-')
      
      a <- length(unique(df$patch[df$patch %in% patchID[1]]))
      s <- length(unique(df$species[df$patch %in% patchID[1]]))
      for (i in 2:2025){
        a <- c(a, length(unique(df$patch[df$patch %in% patchID[1:i]])))
        s <- c(s, length(unique(df$species[df$patch %in% patchID[1:i]])))
      }
      
      df3a$a <- a
      df3a$s <- s
      
      ggplot(df3a, aes(x=A, y=a)) + geom_point()
      ggplot(df3a, aes(x=S, y=s)) + geom_point()
      
      df3 <- rbind(df3, df3a)
    }
    
    
  }
}
df3 <- df3[df3$simnr > 0,]

p1 <- ggplot(df3, aes(x=A, y=a, color=factor(simnr))) + 
  geom_line() + 
  geom_abline(slope=1, intercept=0, linetype=3) + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  theme_bw()
p1

p2 <- ggplot(df3, aes(x=S, y=s, color=factor(simnr))) + 
  geom_line() + 
  geom_abline(slope=1, intercept=0, linetype=3) + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  theme_bw()
p2
ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE)

p1 <- ggplot(df3, aes(x=A, y=a/A, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p1

p2 <- ggplot(df3, aes(x=A, y=1 - s/S, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p2

ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE)

ggplot(df3, aes(x=a/A, y=s/S, color=factor(simnr))) + 
  geom_point(alpha=0.3) + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss))

# Calculate the Kullback-Leibner divergence per a / A graph:
group  <- paste(df3$simnr, df3$clustering, df3$f_loss, sep='-')
ugroup <- unique(group)

df4 <- data.frame(simnr      = 0,
                  clustering = 0,
                  f_loss     = 0,
                  A          = 0,
                  KL_A       = 0,
                  KL_S       = 0)

for (i in ugroup){
  KL_A = 0
  KL_S = 0
  
  for (j in 1:2025){
    sel <- (group == i)&(df3$A == j)
    
    KL_A       <- KL_A + df3$a[sel]/df3$A[sel]
    KL_S       <- KL_S + df3$s[sel]/df3$S[sel]
    
    simnr      <- unique(df3$simnr[group == i])
    clustering <- unique(df3$clustering[group == i])
    f_loss     <- unique(df3$f_loss[group == i])
    df4        <- rbind(df4, list(simnr, clustering, f_loss, j, KL_A, KL_S))
  }
}
df4 <- df4[df4$simnr > 0,]

p1 <- ggplot(df4, aes(x=A, y=KL_A, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p1

p2 <- ggplot(df4, aes(x=A, y=KL_S, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p2

ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE)

ggplot(df4, aes(x=KL_A, y=KL_S, color=A)) + geom_point()

mod <- lm(KL_S ~ KL_A, data=df4)
summary(mod)

df4$exp_KL_S <- mod$fitted.values

p2 <- ggplot(df4, aes(x=A, y=KL_S, color=factor(simnr))) + 
  geom_line() + 
  geom_line(aes(x=A, y=exp_KL_S, color=factor(simnr)), linetype=2) + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p2

# what if we only use part of the data (as the METE SAR function does as well)?
areas <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
df5 <- df3[df3$a %in% areas,]

group  <- paste(df5$simnr, df5$clustering, df5$f_loss, sep='-')
ugroup <- unique(group)

df6 <- data.frame(simnr      = 0,
                  clustering = 0,
                  f_loss     = 0,
                  A          = 0,
                  KL_A       = 0,
                  KL_S       = 0,
                  SSQ_A      = 0,
                  SSQ_S      = 0)

for (i in ugroup){
  KL_A = 0
  KL_S = 0
  SSQ_A = 0
  SSQ_S = 0
  for (j in areas){
    sel <- (group == i)&(df5$a == j)
    
    KL_A       <- KL_A + df5$a[sel][1]/df5$A[sel][1]
    KL_S       <- KL_S + df5$s[sel][1]/df5$S[sel][1]
    
    SSQ_A      <- SSQ_A + (df3$a[sel][1] - df5$A[sel][1])^2
    SSQ_S      <- SSQ_S + (df3$s[sel][1] - df5$S[sel][1])^2
    
    simnr      <- unique(df5$simnr[group == i])
    clustering <- unique(df5$clustering[group == i])
    f_loss     <- unique(df5$f_loss[group == i])
    A          <- df5$A[sel][1]
    df6        <- rbind(df6, list(simnr, clustering, f_loss, j, KL_A, KL_S, SSQ_A, SSQ_S))
  }
}
df6 <- df6[df6$simnr > 0,]

p1 <- ggplot(df6, aes(x=A, y=KL_A, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p1

p2 <- ggplot(df6, aes(x=A, y=KL_S, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p2

ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE)

ggplot(df6, aes(x=KL_A, y=KL_S, color=factor(clustering))) + geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10')

lKL_S <- log(df6$KL_S)
lKL_A <- log(df6$KL_A)
sel <- (!is.infinite(lKL_A))
mod <- lm(lKL_S ~ lKL_A)
summary(mod)

df6$exp_KL_S <- mod$fitted.values

p2 <- ggplot(df6, aes(x=A, y=KL_S, color=factor(simnr))) + 
  geom_line() + 
  geom_line(aes(x=A, y=exp_KL_S, color=factor(simnr)), linetype=2) + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()
p2





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
filename <- 'composition1.00_1.00_0.70_5500.00_2.00_6.50.txt'


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
  
  table(codes[,2])
  
  
  # eerst maar eens 1 soort:
  df2 <- data.frame(species = vector(length=0),
                    A       = vector(length=0),
                    RMSE_A  = vector(length=0),
                    RMSE_n  = vector(length=0))
  
  for (i in 1:length(species)){
    # we moeten rekening houden met alle plekken waar de soort niet in voorkomt, 
    # en ook welke plekken er wel / geen habitat zijn (dit is niet hetzelfde!)
    # Oftewel, even het hele landschap reconstrueren...
    landscape <- data.frame(patchID = patchID1,
                            habitat = 0,
                            n       = 0)
    
    landscape$habitat[landscape$patchID %in% patchID] <- 1
    landscape$n[landscape$patchID %in% df$patch[df$species == species[i]]] <- df$n[df$species == species[i]]
    
    # we kunnen RMSE gebruiken, zodat de nullen mee kunnen doen?
    n_exp  <- n0[i] * 1 / 2025
    RMSE_n <- sqrt(mean((landscape$n - n_exp)^2))
    RMSE_A <- sqrt(mean((landscape$habitat - 1)^2))
    
    df2a <- data.frame(species = species[i],
                       A       = 1,
                       RMSE_A  = RMSE_A,
                       RMSE_n  = RMSE_n)
    df2 <- rbind(df2, df2a)
    
    codes2 <- codes[,1]
    
    for (j in 2:11){
      codes2 <- paste(codes2, codes[,j])
      uCodes <- unique(codes2)
      
      landscape2 <- data.frame(patchID = sort(uCodes),
                               A       = tapply(codes2, codes2, length),
                               habitat = tapply(landscape$habitat, codes2, sum),
                               n       = tapply(landscape$n, codes2, sum))
      
      RMSE_n <- sqrt(mean((landscape2$n - n_exp)^2))
      RMSE_A <- sqrt(mean((landscape2$habitat - landscape2$A)^2))
      
      df2a <- data.frame(species = species[i],
                         A       = mean(landscape2$A),
                         RMSE_A  = RMSE_A,
                         RMSE_n  = RMSE_n)
      df2 <- rbind(df2, df2a)
    }
  }
  
  ggplot(df2, aes(x=RMSE_A, y=RMSE_n, color=factor(A))) + 
    geom_point() +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10')
  
  
  
  
  
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



