# create landscapes with different habitat covers and fragmentation levels, 
# starting from 1 randomly placed patch. 

# each patch contains 1000 trees. Assuming a tree density of 400 trees per ha 
# (Ter Steege et al. 2003), patch size equals 2.5 ha. 
# The landscape is 45 x 45 = 2025 patches, which thus equals 50 km^2 (appr. 
# 7 x 7 km).

# I can create the model in C++ to increase computation speed? 
# Yes, this does indeed work!  

# Patch creation using weighted sampling: 

for (sim_nr in 1:10){
  for (mu in seq(1, 5, 2)){
    
    print(paste('mu = ', mu))
    
    # start with a uniform, continuous habitat and fragment it...
    if (mu == 2){ mu <- 2.01}
    n <- 45
    x <- rep(1:n, n)
    y <- sort(x)
    h <- rep(0, n*n)
    
    mat_h <- h
    
    # first habitat patch:
    hab_patch    <- sample(1:(n*n), 1)
    h[hab_patch] <- 1
    
    # weights for sampling patch creation:
    w <- rep(0, n*n)
    
    # 2D pareto distribution of weights around the habitat patch:
    lmin <- 1
    lmax <- 50
    
    for (f in seq(0.05, 0.95, 0.05)){
      print(paste('f = ', f))
      
      n_patches <- f * n * n  # number of habitat patches that should be created
      while (sum(h) < n_patches){
        dist <- sqrt((x - x[hab_patch])^2 + (y - y[hab_patch])^2)
        w    <- w + 1/(2*pi) * ((2 - mu) / (lmax^(2-mu) - lmin^(2-mu))) * dist^(-1*mu)
        w    <- w * (h == 0)
        w[hab_patch] <- 0  # focal patch
        
        hab_patch    <- sample(1:(n*n), 1, prob = w)
        h[hab_patch] <- 1
      }
      
      # write the data to a dataframe and then to a file:
      mat_h <- cbind(mat_h, h)
      
    }
    mu[mu == 2.01] <- 2
    
    data <- cbind(x, y, mat_h)
    
    filename <- paste('./results/landscapes/landscapes_', sim_nr, '_', mu, '.txt', sep='')
    write.table(data, filename, append=FALSE, col.names = FALSE, row.names = FALSE)
  }
}

library(ggplot2)
df <- data.frame(x=x,
                 y=y,
                 hab=mat_h[,1],
                 cover=0)

f <- round(seq(0.05, 0.95, 0.05), 2)
for (i in 1:20){
  df2 <- data.frame(x=x,
                    y=y, 
                    hab=mat_h[,i],
                    cover=f[i])
  df <- rbind(df, df2)
}

ggplot(df, aes(x=x, y=y, fill=hab)) + geom_raster() + facet_wrap(vars(cover))
