# find out what is going on with the METE fit:
source('./src/est_RAI.R')

df <- data.frame(n      = vector(length=0),
                 fit    = vector(length=0),
                 f_loss = vector(length=0))

for (f_loss in  seq(0.05, 0.95, 0.05)){
  filename <- paste('results/community_composition/1_2_', f_loss, '_0_0.05.txt', 
                    sep='')
  m <- read.table(filename)
  
  table(m$V6)
  n_communities <- min(tapply(m$V6, m$V6, length))
  m <- m[m$V6 == max(m$V6),]
  m <- m[1:n_communities, 7:1006]
  
  for (j in 1:10){
    i <- sample(1:n_communities, 1)
    species <- as.matrix(m[i,])
    n <- tapply(species, species, length)
    n <- sort(n, decreasing = TRUE)
    
    est_RAI(n,'plot')
    # in reality, only a subsample of all individuals is taken, which results in 
    # more singles than there actually are...
    
    fits <- vector(length=0)
    for (ix in (1:10)*100){
      spec <- sample(species, ix, replace=FALSE)
      n <- tapply(spec, spec, length)
      n <- sort(n, decreasing = TRUE)
      #est_RAI(n,'plot')
      
      difference <- est_RAI(n,'differences')
      fit        <- 1 - 
        sum(difference^2)/sum((n - mean(n))^2)
      
      fits <- c(fits, fit)
    }
    # the smaller the subsample, the better the fit! 
    # is dit ook zo bij echte data? 
    
    df2 <- data.frame(n      = (1:10)*100,
                      fit    = fits,
                      f_loss = f_loss)
    p1 <- ggplot(df2, aes(x = n, y=fit)) + 
      geom_point() + 
      xlab('Sample size') + 
      ylab('METE fit to SAD') + 
      theme_bw()
    print(p1)
    
    df <- rbind(df, df2)
  }
  
}

df$good_fit <- df$fit >= 0.5
df$sample_size <- paste('sample size = ', df$n)
df$sample_size <- factor(df$sample_size, levels = unique(df$sample_size))
ggplot(df, aes(x = f_loss, y=fit, color=good_fit)) + 
  geom_point() + 
  facet_wrap(vars(sample_size)) + 
  xlab('Fraction of habitat lost') + 
  ylab('METE fit to SAD') + 
  theme_bw()

ggplot(df[df$n == 200,], aes(x = f_loss, y=fit)) + 
  geom_point() + 
  xlab('Fraction of habitat lost') + 
  ylab('METE fit to SAD') + 
  theme_bw()


#filename <- 'data/CSIRO_PermanentPlots_TreeMeasurementData.csv'
#df <- read.table(filename, header=TRUE, sep=',')
#df$species <- paste(df$family, df$genus, sep=' ')
#species <- df$species[(df$year > 2012)&((df$epNumber == 'ep40')|(df$epNumber == 'ep19'))]


filename <- 'data/Raevel2012_count_data.txt'
df <- read.table(filename, header=TRUE)
m <- as.matrix(df)
m <- m[,2:98]

rowSums(m)

df <- data.frame(n = vector(length=0),
                  fit = vector(length=0),
                 location = vector(length=0))

for (i in 1:52){
  #species <- m[46,]
  species <- m[i,] #colSums(m[1,])
  spec <- lapply(1:length(species), function(ix){
    rep(ix, species[ix])
  })
  species <- unlist(spec)
  
  n <- tapply(species, species, length)
  n <- sort(n, decreasing = TRUE)
  plot(n)
  
  difference <- est_RAI(n,'plot')
  
  if (length(species) > 100){
    fits <- vector(length=0)
    ns <- vector(length = 0)
    for (ix in 100:length(species)){
      spec <- sample(species, ix, replace=FALSE)
      n <- tapply(spec, spec, length)
      n <- sort(n, decreasing = TRUE)
      
      difference <- est_RAI(n,'differences')
      fit        <- 1 - 
        sum(difference^2)/sum((n - mean(n))^2)
      
      fits <- c(fits, fit)
      ns <- c(ns, ix)
    }
    plot(fits)
    
    df2 <- data.frame(n = ns,
                      fit = fits,
                      location = i)
    p1 <- ggplot(df2, aes(x = n, y=fit)) + 
      geom_point(color = 'red') + 
      xlab('Sample size') + 
      ylab('METE fit to SAD') + 
      xlim(c(100, 1000)) + 
      ylim(c(0, 1)) + 
      theme_bw()
    print(p1)
    
    df <- rbind(df, df2)
  }
}

ggplot(df, aes(x = n, y=fit, color=as.factor(location))) + 
  geom_point(alpha = 0.1) + 
  xlab('Sample size') + 
  ylab('METE fit to SAD') + 
  xlim(c(100, 1000)) + 
  ylim(c(0, 1)) + 
  theme_bw() + 
  theme(legend.position = 'none')

 

group <- paste(df$epNumber, df$year, sep='-')
sapply(group, function(x) {
  spec <- df$species[group == x]
  length(unique(spec))
} )



