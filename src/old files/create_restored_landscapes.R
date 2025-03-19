# Create the different spatial configurations of habitat restoration. 
# Each landscape file that is loaded here, consists of 20 landscapes, 
# differing in their % habitat cover (from 100 to 5%). Per landscape, 
# habitat is restored, with % restored habitat cover = 5, 10, 15% etc., 
# until the landscape is once more completely covered. 
#
# Hence, per landscape file, this results in: 
# 190 restoration landscapes per spatial configuration (+ full cover). 
# Using random, fractal, and clustered restoration, this totals to
# 570 restoration landscapes per landscape file. 
# There are 10 x 3 = 30 landscape files; hence the total number of restoration 
# landscapes will be: 
570 * 30
# 17100! Unless we will only use one simulation per parameter combination, 
# then it will 'just' be 1710 restoration landscapes to simulate... 

library(ggplot2)

n <- 20
x <- rep(1:n, n)
y <- sort(x)
id <- 1:400

for (sim_nr in 1:10){
  for (clustering in c(1, 3, 5)){
    filename   <- paste('./results/landscapes/landscapes_', sim_nr, '_', 
                        clustering, '.txt', sep='')
    m <- read.table(filename)
    
    df <- data.frame(x      = vector(length=0),
                     y      = vector(length=0), 
                     hab    = vector(length=0),
                     f_loss = vector(length=0),
                     hab_cover = vector(length=0),
                     clustering = vector(length=0))
    
    for (i in 2:20){
      for (mu in c(1, 3, 5))
      {
        hab <- m[,i]
        
        # add weights to the patches depending on distance to habitat patches and 
        # the amount of clustering of restoration:
        
        # weights for sampling patch destruction:
        w <- rep(0, n*n)
        
        f_loss = 1 - mean(hab)
        
        # 2D pareto distribution of weights around the destroyed patch:
        lmin <- 1
        lmax <- 50
        
        for (j in id[hab == 1]){
          dist <- sqrt((x - x[j])^2 + (y - y[j])^2)
          w    <- w + 1/(2*pi) * ((2 - mu) / (lmax^(2-mu) - lmin^(2-mu))) * dist^(-1*mu)
          w[is.na(w)] <- 0  # focal patch
        }
        w[hab == 1] <- 0
        
        # now, add restored sites:
        hab_cover <- mean(hab)
        for (f in seq(hab_cover, 0.95, 0.05)){
          print(paste('f = ', f))
          
          n_patches <- f * n * n  # number of habitat patches that should remain
          n_restore <- n_patches - sum(hab) # number of patches that needs to be restored
          restored    <- sample(1:(n*n), n_restore, prob = w)
          
          hab[restored] <- 1
          w[hab == 1] <- 0
          
          # write the data to a dataframe and then to a file:
          df2 <- data.frame(x      = x,
                            y      = y, 
                            hab    = hab,
                            f_loss = f_loss,
                            hab_cover  = round(mean(hab), 2),
                            clustering = mu)
          df <- rbind(df, df2)
        }
        
      }
    }
    
    #sel <- round(df$f_loss, 2) == 0.80
    #sel <- df$clustering == 5
    #ggplot(df[sel,], aes(x=x, y=y, fill=as.factor(hab))) + 
    #  geom_raster() + 
    #  scale_fill_manual(values=c('white', 'darkolivegreen4')) + 
    #  facet_grid(cols=vars(hab_cover), rows=vars(f_loss))
    
    # write to a restoration landscape file:
    filename <- paste('./results/landscapes/restored_landscapes_', sim_nr, '_', 
                      clustering, '.txt', sep='')
    write.table(df, filename, col.names = TRUE, row.names = FALSE, append = FALSE)
    
  }
}

