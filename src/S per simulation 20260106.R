# obtain the total number of species in the simulated area per simulation, 
# so we can make a similar figure to figure 3 of our previous paper... 
library(ggplot2)
library(ggpubr)
source('./src/summarySE.R')

get_SAR <- function(clust, loss, simnr){
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
  df$patch  <- paste(df$x, df$y, sep='-')
  uPatch    <- unique(df$patch)
  
  filename2  <- paste('composition1.00_', simnr, '.00_0.00_5500.00_2.00_6.50.txt', sep='')
  filename2  <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename2, sep='')
  df3        <- read.table(filename2)
  names(df3) <- c('x', 'y', 'species', 'n')
  df3        <- df3[df3$x != -1,]
  df3$patch  <- paste(df3$x, df3$y, sep='-')
  
  # Delete all the subcommunities that were lost:
  df3 <- df3[df3$patch %in% uPatch,]
  
  length(unique(df3$species))
  length(unique(df$species))
  
  df2a <- data.frame(clustering = rep(clustering, 2),
                     area_size = rep(2025*(1-f_loss), 2),
                     n_species = c(length(unique(df3$species)),
                                   length(unique(df$species))),
                     static_dynamic = c('Static', 'Dynamic'))
  
  write.table(df2a, 'Fragmented-forest-communities/x64/Release/results/S per simulation 20260106.txt', row.names = FALSE, col.names = FALSE,
              append=TRUE)
}

# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
dat <- expand.grid(clust = c('1.00', '1.50', '2.01','2.50', '3.00', '3.50', '4.00', '4.50', '5.00'), 
                   loss = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
                   simnr = 1:5)

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:nrow(dat), function(i) {
  get_SAR (dat$clust[i], 
           dat$loss[i], 
           dat$simnr[i])
}, future.seed = TRUE)


# plot the results: 
df        <- read.table('Fragmented-forest-communities/x64/Release/results/S per simulation 20260106.txt')
names(df) <- c('clustering', 'area_size', 'n_species', 'static_dynamic')

# we need to add the number of species in case of no habitat loss:

for (simnr in 1:5)
{
  filename2  <- paste('composition1.00_', simnr, '.00_0.00_5500.00_2.00_6.50.txt', sep='')
  filename2  <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename2, sep='')
  df3        <- read.table(filename2)
  names(df3) <- c('x', 'y', 'species', 'n')
  df3        <- df3[df3$x != -1,]
  df3$patch  <- paste(df3$x, df3$y, sep='-')
  
  length(unique(df3$species))
  
  dfa <-  data.frame(clustering = rep(unique(df$clustering), 2),
                     area_size = 2025,
                     n_species = length(unique(df3$species)),
                     static_dynamic = c('Static', 'Dynamic'))
  
  df <- rbind(df, dfa)
  
  # also, add the backward SAR: 
  uPatch <- unique(df3$patch)
  area_size <- sort(unique(df$area_size))
  
  dfa <- data.frame(clustering = 'Backwards SAR',
                    area_size = 2025,
                    n_species = length(unique(df3$species)),
                    static_dynamic = c('Static', 'Dynamic'))
  for (i in area_size){
    S <- length(unique(df3$species[df3$patch %in% uPatch[1:round(i)]]))
    
    dfa2 <- data.frame(clustering = 'Backwards SAR',
                       area_size = i,
                       n_species = S,
                       static_dynamic = c('Static', 'Dynamic'))
    dfa <- rbind(dfa, dfa2)
  }
  
  #ggplot(dfa, aes(x=area_size, y=n_species)) + geom_point()
  
  df <- rbind(df, dfa)
}


sdf <- summarySE(df, measurevar="n_species", 
                 groupvars=c("clustering","area_size", "static_dynamic"))

sdf <- sdf[sdf$clustering %in% c(1, 3, 5, 'Backwards SAR'),]

sdf$mu2 <- 'Random'
sdf$mu2[sdf$clustering == 3] <- 'Fractal'
sdf$mu2[sdf$clustering == 5] <- 'Clustered'
sdf$mu2[sdf$clustering == 'Backwards SAR'] <- 'Backwards SAR'

sdf$static_dynamic2 <- 'Immediately after habitat loss'
sdf$static_dynamic2[sdf$static_dynamic == 'Dynamic'] <- 'After stabilization'
sdf$static_dynamic2 <- factor(sdf$static_dynamic2, levels=c('Immediately after habitat loss', 'After stabilization'))
sdf$hab_loss <- (1 - sdf$area_size/2025)*100

# figures for presentation: 
sel <- (sdf$static_dynamic == 'Static')
ggplot(sdf[sel,], aes(x=hab_loss, y=n_species, color=mu2)) + 
  geom_errorbar(aes(ymin=n_species-sd, ymax=n_species+sd), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_color_manual(values = c('grey40', "#F8766D", "#00BA38", "#619CFF")) + 
  xlab('% Habitat loss') + 
  ylab('Total # species') +
  ylim(c(1500, 5500)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        legend.title=element_blank())

sel <- (sdf$static_dynamic == 'Dynamic')
ggplot(sdf[sel,], aes(x=hab_loss, y=n_species, color=mu2)) + 
  geom_errorbar(aes(ymin=n_species-sd, ymax=n_species+sd), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_color_manual(values = c('grey40', "#F8766D", "#00BA38", "#619CFF")) + 
  xlab('% Habitat loss') + 
  ylab('Total # species') +
  ylim(c(1500, 5500)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        legend.title=element_blank())


dat <- data.frame(Rank = 1:6,
                  Abundance = c(10, 7, 3, 2, 1, 1))

cols = c('grey60', "#F8766D", "#00BA38", "#619CFF", "plum3", "gold")
ggplot(dat, aes(x=Rank, y=Abundance)) + 
  geom_col(col='black', fill=cols) + 
  ylim(c(0, 10)) + 
  theme_bw()

sum(dat$Abundance)
p <- dat$Abundance / 24

N0 = 24
A0 = 100
dat2 <- data.frame(S = 6, A = 100)

for (A in seq(0, 0.95, 0.05)*A0){
  S    <- sum(1 - (1 - p)^(A*(N0/A0)))
  dat2 <- rbind(dat2, list(S = S, A = A)) 
  
}
dat2 <- dat2[dat2$A < 100,]

ggplot(dat2, aes(x=A, y=S)) + geom_line() + 
  geom_point() + 
  xlab('Area size (A)') + 
  ylab('Number of species (S)') + 
  theme_bw()

dat3 <- data.frame(A = 1/c(1, 2, 4, 8, 16, 32, 64),
                   S = c(6, 4, 4, 4, 2, 1, 1))
ggplot(dat3, aes(x=A, y=S)) + geom_point() + 
  geom_line() + 
  xlab('Area size (A)') + 
  ylab('Number of species (S)') + 
  theme_bw()

df <- read.table('Fragmented-forest-communities/x64/Release/results/SAR per simulation 20250709.txt')
names(df) <- c('S', 'A', 'type', 'simnr', 'f_loss', 'clustering')

sdf <- summarySE(df, measurevar="S", 
                 groupvars=c("clustering","A", "type", "f_loss"))

sdf <- sdf[sdf$clustering %in% c(1, 3, 5),]

sdf$mu2 <- 'Random'
sdf$mu2[sdf$clustering == 3] <- 'Fractal'
sdf$mu2[sdf$clustering == 5] <- 'Clustered'

sdf$f_loss2 <- paste(sdf$f_loss*100, '% habitat loss', sep='') 
sel <- (sdf$type == 'dynamic')&(sdf$f_loss %in% c(0.2, 0.5, 0.8))
ggplot(sdf[sel,], aes(x=A, y=S, color=mu2)) + 
  geom_errorbar(aes(ymin=S-sd, ymax=S+sd), width=0) + 
  geom_line() + 
  geom_point() +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) + 
  xlab('Area size (A)') + 
  ylab('Number of species (S)') + 
  facet_grid(rows=vars(f_loss2)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())
