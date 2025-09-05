
for (clust in c('1.00', '1.50', '2.01','2.50', '3.00', '3.50', '4.00', '4.50', '5.00')){
  for (simnr in 1:5){
    for (loss in 1:9/10){
      filename <- paste('composition', clust,'_', simnr, 
                        '.00_', loss, '0_5500.00_2.00_6.50.txt', sep='')
      
      # derive the simulation number, clustering, and habitat loss from the filename:
      dat        <- gsub("composition", "", filename)
      dat        <- unlist(strsplit(dat, '_'))
      
      clustering <- as.numeric(dat[1])
      sim_nr     <- as.numeric(dat[2])
      f_loss     <- as.numeric(dat[3])
      
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
      df        <- read.table(filename2)
      names(df) <- c('x', 'y', 'species', 'n')
      df        <- df[df$x != -1,]

      # to make sure that there is at least 1 subcommunity present in de 10km2 area,
      # start with the first subcommunity in the list: 
      minx <- min(df$x[1], 25)
      miny <- min(df$y[1], 25)
      
      df        <- df[(df$x >= minx)&(df$x < (minx+20))&
                        (df$y >= miny)&(df$y < (miny+20)),]
      df$patch  <- paste(df$x, df$y, sep='-')
      
      uPatch    <- unique(df$patch)
      
      filename2  <- paste('composition1.00_', simnr, '.00_0.00_5500.00_2.00_6.50.txt', sep='')
      filename2  <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename2, sep='')
      df3        <- read.table(filename2)
      names(df3) <- c('x', 'y', 'species', 'n')
      df3        <- df3[df3$x != -1,]
      df3$patch  <- paste(df3$x, df3$y, sep='-')
      
      cover <- length(uPatch)/length(unique(df3$patch))
      
      # Delete all the subcommunities that were lost:
      df3 <- df3[df3$patch %in% uPatch,]
      
      species <- rep(df3$species, df3$n)
      
      n <- tapply(species, species, length)
      n <- sort(n, decreasing = TRUE)
      
      plot(n)
      
      n <- c(n, rep(0, 5500 - length(n)))
      
      data <- data.frame(type = 'static',
                         clustering = clustering, 
                         f_loss = f_loss, 
                         abundance = n,
                         rank = 1:length(n),
                         cover = cover)
      
      # After community dynamics have stabilized:
      species <- rep(df$species, df$n)
      
      n <- tapply(species, species, length)
      n <- sort(n, decreasing = TRUE)
      n <- c(n, rep(0, 5500 - length(n)))
      
      data2 <- data.frame(type = 'dynamic',
                         clustering = clustering, 
                         f_loss = f_loss, 
                         abundance = n,
                         rank = 1:length(n),
                         cover = cover)
      data <- rbind(data, data2)
      
      filename <- paste( 'Fragmented-forest-communities/x64/Release/results/S0 per simulation 20250902/simulated SADs 20250902_simnr=',
                         sim_nr, 'in10km2.txt', sep='')
      write.table(data, filename, append = TRUE, row.names = FALSE, 
                  col.names = FALSE)
    }
  }
}

sim_nr <- 1
filename <- paste( 'Fragmented-forest-communities/x64/Release/results/S0 per simulation 20250902/simulated SADs 20250902_simnr=',
                   sim_nr, 'in10km2.txt', sep='')
df <- read.table(filename)

for (sim_nr in 2:5){
  filename <- paste( 'Fragmented-forest-communities/x64/Release/results/S0 per simulation 20250902/simulated SADs 20250902_simnr=',
                     sim_nr, 'in10km2.txt', sep='')
  df <- rbind(df, read.table(filename))
}

names(df) <- c('type','clustering', 'f_loss', 'abundance', 'rank', 'cover')

# cover klopt niet!! 
df$cover <- df$cover * 5.0625

ggplot(df[df$type == 'dynamic',], aes(x=rank, color=round(clustering,1), y=abundance)) + 
  geom_point(alpha=0.3) + 
  scale_color_viridis_c(name='Degree of habitat clustering (µ)') + 
  ylab('Abundance') +
  xlab('Rank') + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  guides (color = guide_colourbar(barwidth = 25)) + 
  facet_wrap(vars(round(cover, 1))) + 
  theme_bw() + 
  theme(legend.position = 'top', 
        strip.placement = "outside", 
        strip.background = element_blank())

group <- paste(df$type, round(df$clustering, 1), round(df$cover, 1), df$rank, sep='-')
sdf   <- data.frame(type = tapply((df$type == 'static'), group, mean),
                    clustering = tapply(round(df$clustering, 1), group, mean),
                    cover = tapply(round(df$cover, 1), group, mean),
                    rank = tapply(df$rank, group, mean),
                    abundance = tapply(df$abundance, group, mean)) 

ggplot(sdf[sdf$type == 0,], aes(x=rank, color=round(clustering,1), y=abundance)) + 
  geom_point(alpha=0.2) + 
  scale_color_viridis_c(name='Degree of habitat clustering (µ)') + 
  ylab('Abundance') +
  xlab('Rank') + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  guides (color = guide_colourbar(barwidth = 25)) + 
  facet_wrap(vars(round(cover, 1))) + 
  theme_bw() + 
  theme(legend.position = 'top', 
        strip.placement = "outside", 
        strip.background = element_blank())

ggplot(sdf[sdf$type == 0,], aes(x=rank, color=cover*100, y=abundance)) + 
  geom_point(alpha=0.3) + 
  scale_color_viridis_c(name='% habitat cover') + 
  ylab('Abundance') +
  xlab('Rank') + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  guides(color = guide_colourbar(barwidth = 25)) + 
  facet_wrap(vars(clustering)) + 
  theme_bw() + 
  theme(legend.position = 'top', 
        strip.placement = "outside", 
        strip.background = element_blank())

# Let's calculate densities instead of abundances:
sdf$density <- NA
group <- paste(sdf$type, sdf$clustering, sdf$cover, sep='-')

for (i in unique(group)){
  sel   <- group == i
  N0    <- sum(sdf$abundance[sel])
  sdf$density[sel] <- sdf$abundance[sel] / N0
}

ggplot(sdf[sdf$type == 0,], aes(x=rank, color=cover*100, y=density)) + 
  geom_point(alpha=0.3) + 
  scale_color_viridis_c(name='% habitat cover') + 
  ylab('Abundance') +
  xlab('Rank') + 
  #scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  guides(color = guide_colourbar(barwidth = 25)) + 
  facet_wrap(vars(clustering)) + 
  theme_bw() + 
  theme(legend.position = 'top', 
        strip.placement = "outside", 
        strip.background = element_blank())


sel <- (sdf$type == 0)&(sdf$cover == 0.3)&(sdf$clustering == 4)&(sdf$abundance > 0)
ggplot(sdf[sel,], aes(x=log(rank), y=log(density))) + geom_point() + 
  geom_smooth(method = "gam")

y = log(sdf$density[sel])

# guesses:
a2 = 1.5
a3 = 0.25

fitFun <- function(y, a2, a3){
  y  = sort(y, decreasing = TRUE)
  a1 = max(y) + a2
  x  = log(1:length(y))
  y2 = a1 - a2*exp(a3*x)
  fit <- 1 - sum((y - y2)^2)/sum((y - mean(y2))^2)
  return(fit)
}

fitFun(y, a2, a3)

plot(x,y, log='')
points(x,y2, col='red')

#bbmle-likelihood fit
library(stats)
library(bbmle)

m0= mle2(fitFun(y, a2,a3),data=list(y=y),method="L-BFGS-B",
         lower=c(a2=0, a3=0),
         start=list(a2=1.5, a3=0.25))

coef(m0)
