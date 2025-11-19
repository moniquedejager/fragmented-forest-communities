# From SAD to SAR: 

#libraries:
library(ggplot2)
library(meteR)

simnr <- 1
filename2  <- paste('composition1.00_', simnr, '.00_0.00_5500.00_2.00_6.50.txt', sep='')
filename2  <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename2, sep='')
df         <- read.table(filename2)
names(df)  <- c('x', 'y', 'species', 'n')

df         <- df[df$x != -1,]
df$patch   <- paste(df$x, df$y, sep='-')

species    <- rep(df$species, df$n)

n <- tapply(species, species, length)
n <- sort(n, decreasing = TRUE)

# SAD:
plot(n, log='xy')

# shannon diversity:
p <- n/sum(n)
E <- -1 * sum(p * log(p))

esf1 <- meteESF(S0 = length(n), N0 = sum(n))
sad1 <- sad(esf1)

plot(sad1$d(1:length(n))*sum(n), log='xy', xlab='Rank', ylab='Abundance')
p_mete <- sad1$d(1:length(n))
E_mete <- -1 * sum(p_mete * log(p_mete))

# Mete's Shannon entropy is a bit larger than that of our simulation data... 
E_mete - E

# Now, calculate SAR from SAD:
areas <- 2^(1:20)

SAR <- data.frame(A      = 1, 
                  S      = 1, 
                  type   = c('Using p from simulations', 'Using p from METE'))

for (i in areas){
  SAR <- rbind(SAR, list(i, sum(1 - (1 - p)^i), 'Using p from simulations'))
  SAR <- rbind(SAR, list(i, sum(1 - (1 - p_mete)^i), 'Using p from METE'))
}

# SAR made with meteR function:
sar <- meteSAR(Amin=2, 
               A0=sum(n), 
               S0=length(n), 
               N0=sum(n))
SAR <- rbind(SAR, list(sar$pred$A, 
                       sar$pred$S, 
                       rep('Using meteR SAR function', length(sar$pred$A))))

# SAR made with the simulation data: 
for (j in areas){
  S <- vector(length=0)
  for (i in 1:10){
    spec <- sample(species, j, replace = FALSE)
    S    <- c(S, length(unique(spec)))
    SAR  <- rbind(SAR, list(j, length(unique(spec)), 'Using simulation data'))
  }
}

ggplot(SAR, aes(x=A/1000*0.025, y=S, color=type)) + 
  geom_point() + 
  geom_line() + 
  xlab('Area (km2)') + 
  ylab('Species') + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  theme_bw()
 



sel <- SAR$type == 'Using p from METE' 
ggplot(SAR[sel,], aes(x=A/1000*0.025, y=S)) + 
  geom_line() +
  geom_point() + 
  xlab('Area (km2)') + 
  ylab('Species') + 
  theme_bw()

# Biodiversity loss directly after habitat loss: 
# p per species decreases with habitat loss: p_new = p * cover

SAR2 <- SAR[SAR$A == -1,]

cover <- 1
p_new <- p * cover

for (i in areas){
  SAR2 <- rbind(SAR2, list(A=i, S=sum(1 - (1 - p_new)^i), type=paste(cover*100, '% Habitat cover')))
}

ggplot(SAR2, aes(x=A/1000*0.025, y=S, color=type)) + 
  geom_point() + 
  geom_line() + 
  xlab('Area (km2)') + 
  ylab('Species') + 
  scale_color_discrete(name='') + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()







