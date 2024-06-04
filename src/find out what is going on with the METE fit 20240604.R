# find out what is going on with the METE fit:
source('./src/est_RAI.R')

filename <- 'results/community_composition/1_4_0.05_0_0.05.txt'
m <- read.table(filename)

table(m$V6)
#m <- m[m$V6 == 250,]
m <- m[1:380, 7:1006]
m <- as.matrix(m)

i <- 3
species <- m[i,]
n <- tapply(species, species, length)
n <- sort(n, decreasing = TRUE)
plot(n)

difference <- est_RAI(n,'plot')


# in reality, only a subsample of all individuals is taken, which results in 
# more singles than there actually are...


fits <- rep(0, 1000)
for (ix in 100:1000){
  spec <- sample(species, ix, replace=FALSE)
  n <- tapply(spec, spec, length)
  n <- sort(n, decreasing = TRUE)
  #est_RAI(n,'plot')
  
  difference <- est_RAI(n,'differences')
  fit        <- 1 - 
    sum(difference^2)/sum((n - mean(n))^2)
  
  fits[ix] <- fit
}
plot(fits)

# the smaller the subsample, the better the fit! 
# is dit ook zo bij echte data? 

df <- data.frame(n = 100:1000,
                 fit = fits[100:1000])
ggplot(df, aes(x = n, y=fit)) + 
  geom_point() + 
  xlab('Sample size') + 
  ylab('METE fit to SAD') + 
  theme_bw()


filename <- 'data/CSIRO_PermanentPlots_TreeMeasurementData.csv'
df <- read.table(filename, header=TRUE, sep=',')

df$species <- paste(df$family, df$genus, sep=' ')

species <- df$species[(df$year > 2012)&((df$epNumber == 'ep40')|(df$epNumber == 'ep19'))]
n <- tapply(species, species, length)
n <- sort(n, decreasing = TRUE)
plot(n)

difference <- est_RAI(n,'plot')

fits <- rep(0, length(species))
for (ix in 100:length(species)){
  spec <- sample(species, ix, replace=FALSE)
  n <- tapply(spec, spec, length)
  n <- sort(n, decreasing = TRUE)
  #est_RAI(n,'plot')
  
  difference <- est_RAI(n,'differences')
  fit        <- 1 - 
    sum(difference^2)/sum((n - mean(n))^2)
  
  fits[ix] <- fit
}
plot(fits)

df2 <- data.frame(n = 100:length(species),
                 fit = fits[100:length(species)])
ggplot(df2, aes(x = n, y=fit)) + 
  geom_point(color = 'red') + 
  xlab('Sample size') + 
  ylab('METE fit to SAD') + 
  xlim(c(100, 1000)) + 
  ylim(c(0, 1)) + 
  theme_bw()

table(df$epNumber, df$year)

 

group <- paste(df$epNumber, df$year, sep='-')
sapply(group, function(x) {
  spec <- df$species[group == x]
  length(unique(spec))
} )



