# In this R-script, I am going to try to estimate the SAD for a large area, 
# in which species are NOT homogeneously distributed. 

# For each simulated environment, we obtain the actual SAD for the entire area.
# We also calculate the mean S0 per habitat patch. 
# With this S0 (and N0 = 1000 individuals per habitat patch), we can estimate 
# METEs Lagrange multiplier lambda1. But this is not needed yet,

# even mijn gedachten ordenen, want het is een zooitje in mijn hoofd... 
# Voor een habitat patch kun je met METE de SAD goed schatten, omdat er 
# geen ruimtelijke structuur in zit en dus alle individuen van alle soorten
# homogeen over de (niet-bestaande) ruimte verdeeld zijn. 
# Als er meerdere patches zijn, dan zijn er ook meerdere SADs, waarbij de kans 
# dat een soort veel voorkomt in alle patches weer afhangt van hoe goed de 
# patches met elkaar verbonden zijn. Bij hoge connectiviteit (zoals wanneer 
# er 100% cover is), is de kans groter dat een soort in ongeveer dezelfde 
# aantallen voorkomt in de patches (dit is dan soort van Poisson verdeeld). 
# Maar als er weinig connectiviteit is, dan is de kans even groot dat een soort 
# veel of weinig voorkomt. De kans bestaat ook dat de soort helemaal niet in 
# andere patches voorkomt. 

# mijn plan is om per simulatie, per habitat patch, en per soort, na te gaan
# wat zijn ICFD- of Quantile Q(n) waarde is. Daarna maak ik een QS(n) voor de 
# abundantie van deze soort in de andere patches. Door deze logit te 
# transformeren, kan ik met een lineair model de vorm ervan schatten. Hiervan
# sla ik dan ook de slope en intercept op. OF: ik kan anders ook eerst alle QS(n)
# van de abundantie per soort opslaan, en dan de relatie bepalen tussen deze
# QS(n), de Quantile waarde van de soort in de patch, en de abundantie in
# andere patches. 

# Dan moet ik alleen nog een manier vinden om dit ook nog aan connectiviteit te 
# koppelen. Het wordt wel een beetje een rommelig zooitje zo... 
# Maar laten we eerst maar even kijken of het voor 1 simulatie werkt, voordat 
# we de connectiviteit er ook bij gaan betrekken. 

clust <- '1.00'  
simnr <- 1 
loss  <- '0.0' 

filename <- paste('composition', clust,'_', simnr, 
                  '.00_', loss, '0_5500.00_2.00_6.50.txt', sep='')
filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
df        <- read.table(filename2)
names(df) <- c('x', 'y', 'species', 'n')
df        <- df[df$x != -1,]
df$patch  <- paste(df$x, df$y, sep='-')
uPatch    <- unique(df$patch)

# collect the data per habitat patch:
df$Q <- NA
for (iPatch in 1:length(uPatch)){
  print(iPatch / length(uPatch) * 100)
  patch  <- df[df$patch == uPatch[iPatch],]
  patch$n <- patch$n + runif(length(patch$n), -0.1, 0.1)
  # create the quantile function for this patch:
  n       <- sort(patch$n)
  
  patchN <- n
  # include all species with zero abundance as well: 
  #patchN <- c(rep(0, 5500 - length(n)), n)
  patchQ  <- (length(patchN):1)/length(patchN)
  #plot(patchN,patchQ, type='l')
  
  # collect the data per species:
  for (iSpec in 1:nrow(patch)){
    nSpec <- patch$n[iSpec]
    qSpec <- max(patchQ[patchN == nSpec])
    patch$Q[iSpec] <- qSpec
  }
  df$Q[df$patch == uPatch[iPatch]] <- patch$Q
}

library(ggplot2)

df$logitQ <- log(df$Q/(1-df$Q))
ggplot(df, aes(x=Q)) + geom_histogram(bins=20)

mean(df$Q[df$n == 1])
h <- hist(df$Q, breaks=20)

df2 <- data.frame(Q1   = vector(length=0),
                   Q2   = vector(length=0), 
                   P_Q2 = vector(length=0),
                   P0   = vector(length=0))

for (i2 in 1:20){
  # select all species with a certain Q(n):
  spec <- df$species[(df$Q > h$breaks[i2])&(df$Q <= h$breaks[i2+1])]
  w    <- tapply(spec, factor(spec), length)
  spec <- tapply(spec, factor(spec), mean)
  
  # because it concerns a lot of species, we will have to do random sampling
  # to speed things up... 
  sample_Qs <- vector(length=0)
  for (i in 1:length(spec)){
    q <- df$Q[df$species == spec[i]]
    # add the patches in which the species is absent (n = 0, Q = 1):
    q <- c(q, rep(1, length(uPatch) - length(q))) 
    sq <- sample(q, 10*w[i], replace=TRUE)
    sample_Qs <- c(sample_Qs, sq)
  }
  P0 <- length(sample_Qs[sample_Qs == 1]) / length(sample_Qs) # probability of zero individuals on a patch
  h2 <- hist(sample_Qs[sample_Qs < 1], breaks=seq(0, 1, 0.05))
  plot(h2$mids, h2$counts/sum(h2$counts))
  
  df2a <- data.frame(Q1   = h$mids[i2],
                     Q2   = h2$mids, 
                     P_Q2 = h2$counts/sum(h2$counts),
                     P0   = P0)
  df2 <- rbind(df2, df2a)
}
ggplot(df2, aes(x=Q1, y=Q2, fill=P_Q2)) + 
  geom_raster() +
  scale_fill_viridis_c() + 
  theme_bw()
  #ggplot(df, aes(x=n, y=Q)) + geom_point()

ggplot(df2, aes(x=log(Q1/(1-Q1)), y=P_Q2, color=factor(Q2))) + 
  geom_line() +
  scale_color_viridis_d() + 
  theme_bw()

# CFD per Q1: 
df2$CFD <- NA
for (i in unique(df2$Q1)){
  df3 <- df2[df2$Q1 == i,]
  CFD <- df3$P_Q2[1]
  for (j in 2:20){
    CFD <- c(CFD, CFD[j-1] + df3$P_Q2[j])
  }
  df2$CFD[df2$Q1 == i] <- CFD
}
ggplot(df2, aes(x=Q1, y=CFD, color=factor(Q2))) + 
  geom_line() + 
  scale_color_viridis_d() + 
  theme_bw()

ggplot(df2, aes(x=Q2, y=CFD, color=factor(Q1))) + 
  geom_line() + 
  scale_color_viridis_d() + 
  theme_bw()

mean = 1
rate = 10
shape = mean * rate
x <- rgamma(1000, shape, rate)
hist(x, breaks=100)

# collect the data per species:
df$intercept <- NA
df$slope     <- NA
df$r2        <- NA
uSpec <- unique(df$species)
for (iSpec in 1:length(uSpec)){
  print(iSpec)
  df2  <- df[(df$species == uSpec[iSpec]),]
  
  n        <- sort(df2$n)
  #addPatch <- length(uPatch[(uPatch %in% df2$patch) == FALSE])
  #n        <- c(rep(0, addPatch), n)
  #Q        <- (length(n):1)/length(n)  
  
  q <- sort(df2$Q)
  addPatch <- length(uPatch[(uPatch %in% df2$patch) == FALSE])
  q        <- c(q, rep(1, addPatch))
  hist(q)
  
  # for each patch, what are the differences in q between this patch and others?
  
  ix <- 1:length(uPatch)
  dq <- vector(length=0)
  for (iPatch in ix){
    dq <- c(dq, q[iPatch] - q[ix != iPatch])
  }
  hist(dq)
  
  
  Q        <- (1:length(q))/length(q)
  Q        <- Q[q < 1]
  q        <- q[q < 1]
  
  logitq <- log(q/(1-q))
  logitQ <- log(Q/(1-Q))
  plot(logitq, logitQ)
  
  sel <- (!is.infinite(logitQ))
  if (length(q[sel]) > 1){
    mod  <- lm(logitQ[sel]~logitq[sel])
    smod <- summary(mod)
    
    df$intercept[df$species == uSpec[iSpec]] <- mod$coefficients[1]
    df$slope[df$species == uSpec[iSpec]]     <- mod$coefficients[2]
    df$r2[df$species == uSpec[iSpec]]        <- smod$adj.r.squared
  } else {
    df$intercept[df$species == uSpec[iSpec]] <- logitQ
    df$slope[df$species == uSpec[iSpec]]     <- NA
    df$r2[df$species == uSpec[iSpec]]        <- NA
  }
}

df$logitQ <- log(df$Q / (1 - df$Q))
library(ggplot2)
ggplot(df, aes(x=logitQ, y=slope, color=intercept)) + 
  geom_point() + 
  scale_y_continuous(trans='log10') + 
  scale_color_viridis_c(trans='log10')

ggplot(df, aes(x=intercept, y=slope, color=logitQ)) + 
  geom_point() + 
  scale_color_viridis_c()

ggplot(df, aes(x=logitQ, y=intercept, color=slope)) + 
  geom_point() + 
  scale_color_viridis_c()

mod <- lm(log(intercept)~Q,data=df[df$r2 > 0.8,], weights=r2)
summary(mod)

# summarize per species
sdf <- data.frame(species   = tapply(df$species, factor(df$species), mean),
                  Q         = tapply(df$Q, factor(df$species), mean),
                  intercept = tapply(df$intercept, factor(df$species), mean),
                  slope     = tapply(df$slope, factor(df$species), mean),
                  r2        = tapply(df$r2, factor(df$species), mean))
sdf$Q[sdf$Q == 1] <- 0.999
sdf$logitQ <- log(sdf$Q/(1-sdf$Q))

ggplot(sdf[sdf$r2 > 0.0,], aes(x=Q, y=-slope, color=intercept)) + 
  geom_point() + 
  scale_color_viridis_c() +
  scale_y_continuous(trans='log10') 


ggplot(sdf[sdf$r2 > 0.0,], aes(x=logitQ, y=1/slope, color=r2)) + 
  geom_point() + 
  scale_color_viridis_c() 

ggplot(sdf[sdf$r2 > 0.0,], aes(x=logitQ, y=intercept, color=r2)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'R2') + 
  xlab('logit Q(n)') + 
  ylab('Intercept (α)') + 
  theme_bw() + 
  theme(legend.position = 'top',
        legend.key.width= unit(1, 'cm'))

ggplot(sdf[sdf$r2 > 0.0,], aes(x=intercept, y=1/slope, color=r2)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'R2') + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  xlab('Intercept (α)') + 
  ylab('Slope (β)') + 
  theme_bw() + 
  theme(legend.position = 'top',
        legend.key.width= unit(1, 'cm'))

ggplot(sdf[sdf$r2 > 0.0,], aes(y=1/slope, x=logitQ, color=r2)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'r2') + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  xlab('logit Q(n)') + 
  ylab('Slope (β)') + 
  theme_bw() + 
  theme(legend.position = 'top',
        legend.key.width= unit(1, 'cm'))

smod_a <- summary(lm(intercept~logitQ,data=sdf,weights=r2))
smod_b <- summary(lm(1/slope~logitQ, data=sdf,weights=r2))

# plot it:
df2 <- expand.grid(Q = seq(0.001, 0.999, 0.001),
                   n = 0:300)
df2$logitQ <- log(df2$Q/(1-df2$Q))
df2$a <- smod_a$coefficients[1,1] + smod_a$coefficients[2,1]*df2$logitQ
df2$b <- 1/(smod_b$coefficients[1,1] + smod_b$coefficients[2,1]*df2$logitQ)
df2$b[df2$b > 0] <- 1000
df2$y <- df2$a + df2$b*df2$n
df2$y <- exp(df2$y)/(1+ exp(df2$y))
df2$y[is.na(df2$y)] <- 0

#ggplot(df2[df2$b > 0,], aes(x=Q, y=b)) + geom_point()

ggplot(df2, aes(x=n, y=Q, fill=y)) +
  geom_raster() + 
  scale_fill_viridis_c(limits=c(0,1)) + 
  theme_bw()+ 
  theme(legend.position = 'top',
        legend.key.width= unit(1, 'cm'))

ggplot(df2[df2$Q %in% seq(0, 0.8, 0.05),], aes(x=n, y=y, color=factor(Q))) + 
  geom_line() + 
  scale_color_viridis_d(name = 'Q(n)') + 
  scale_y_continuous() + 
  theme_bw() + 
  theme(legend.position = 'right')

# per species, you randomly select a quantile. For this quantile, 
# you randomly select from the second quantile function to sample the abundances
# per patch. You add all abundances together to get an overall abundance per 
# species. 

S0       <- length(unique(df$species))
nPatches <- length(uPatch)

# per species, calculate the estimated abundance per patch and subsequently
# calculate the total estimated abundance for the entire area:
a1 <- smod_a$coefficients[1,1]
a2 <- smod_a$coefficients[2,1]
b1 <- smod_b$coefficients[1,1]
b2 <- smod_b$coefficients[2,1]

SAD <- expand.grid(rank = 1:S0,
                  abundance = 1)
for (i in 1:S0){
  Q1       <- 1/S0*i 
  Q1[Q1 == 1] <- 0.999
  logitQ1  <- log(Q1/(1 - Q1))
  a        <- a1 + a2*logitQ1
  b        <- 1/(b1 + b2*logitQ1)
  b[b > 0] <- 1000
  
  Q2      <- 1/nPatches * (1:nPatches)
  Q2[Q2 == 1] <- 0.999
  logitQ2 <- log(Q2 / (1 - Q2))
  n       <- round((logitQ2 - a)/b) 
  
  SAD$abundance[i] <- sum(n)
}

ggplot(SAD, aes(x=rank, y=abundance)) + 
  geom_point()

n <- tapply(df$n, df$species, sum)
n <- sort(n, decreasing = TRUE)

SAD2 <- data.frame(rank = 1:S0,
                   abundance = n,
                   type = 'Simulation data')
SAD$type <- 'Estimate'
SAD <- rbind(SAD, SAD2)

sum(n)
sum(SAD$abundance[SAD$type == 'Estimate'])

ggplot(SAD, aes(x=rank, y=abundance, color=type)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10')


# what if we use random sampling?
SAD2 <- expand.grid(rank = 1:S0,
                    abundance = 1)
for (i in 1:S0){
  Q1 <- runif(1)
  Q1[Q1 == 1] <- 0.999
  logitQ1  <- log(Q1/(1 - Q1))
  a        <- a1 + a2*logitQ1
  b        <- 1/(b1 + b2*logitQ1)
  b[b > 0] <- 1000
  
  Q2      <- runif(nPatches)
  Q2[Q2 == 1] <- 0.999
  logitQ2 <- log(Q2 / (1 - Q2))
  n       <- round((logitQ2 - a)/b) 
  
  SAD2$abundance[i] <- sum(n)
}

ix   <- order(SAD2$abundance, decreasing = TRUE)
SAD2 <- SAD2[ix,]
SAD2$rank <- 1:nrow(SAD2)
ggplot(SAD2, aes(x=rank, y=abundance)) + 
  geom_point()
SAD2$type <- 'Random sampling'

SAD <- rbind(SAD, SAD2)

ggplot(SAD, aes(x=rank, y=abundance, color=type)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10')

##############################
# vanaf hier de SAR berekenen en kijken hoe goed hij fit!

# try it out:
# first, create the SAD per patch:
S0 <- mean(tapply(df$species, df$patch, length))
N0 <- 1000

# we have to find the value of lambda1, the Lagrange multiplier we need
# to calculate the SAD:
# lambda1 is equivalent to fisher's alpha / N0
a  <- S0
f  <- a * log(1 + N0/a) - S0
af <- data.frame(a = S0, f = f)
while(f > 0.01){
  a <- a - 1
  f <- a * log(1 + N0/a) - S0  
  af <- rbind(af, list(a = a, f = f))
}
# smaller steps:
a  <- a + 1
f  <- a * log(1 + N0/a) - S0
af <- data.frame(a = a, f = f)
while(f > 0.01){
  a <- a - 0.01
  f <- a * log(1 + N0/a) - S0  
  af <- rbind(af, list(a = a, f = f))
}
lambda1 <- a / N0

# estimate the number of species with abundance n:
n <- 1:(N0 - S0 + 1)
P <- 1 / log(lambda1^-1) * exp(-lambda1 * n)/n
P <- P/sum(P)

# create an ICFD / a quantile function of this: 
x    <- 1
Q    <- rep(0, length(n))
for (i in 1:length(n)){
  x       <- x - P[i]
  Q[i]    <- x
}
n <- n[Q > 0]
Q <- Q[Q > 0]
plot(n, Q, type='l', log='xy')

logitQ <- log(Q/(1-Q))
a <- smod_a$coefficients[1,1] + smod_a$coefficients[2,1]*logitQ
b <- 1/(smod_b$coefficients[1,1] + smod_b$coefficients[2,1]*logitQ)
y <- a + b*n

df2 <- data.frame(n = n,
                  Q = Q,
                  y = y)

ggplot(df2, aes(x=n, y=Q, color=y))  





        