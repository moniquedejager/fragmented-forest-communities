# Edwin's comment: but this probably depends on the amount of fragmentation, 
# I can imagine there is like a nested relationship that even for spatial 
# randomly distributed species at some point this will also become a 
# negative effect.


df <- data.frame(mu = vector(length=0),
                 n_ind = vector(length=0),
                 hab_cover = vector(length=0))

for (sim_nr in 1:10){
  
  X <- sample(1:20, 50, replace=TRUE)
  Y <- sample(1:20, 50, replace=TRUE)
  species_loc <- paste(X, Y, sep='-')
  plot(X, Y)
  
  for (mu in seq(1, 5, 0.1)){
    
    print(paste('mu = ', mu))
    
    # start with a uniform, continuous habitat and fragment it...
    if (mu == 2){ mu <- 2.01}
    n <- 20
    x <- rep(1:n, n)
    y <- sort(x)
    h <- rep(1, n*n)
    
    # first destroyed patch:
    destroyed    <- sample(1:(n*n), 1)
    h[destroyed] <- 0
    
    # weights for sampling patch destruction:
    w <- rep(0, n*n)
    
    # 2D pareto distribution of weights around the destroyed patch:
    lmin <- 1
    lmax <- 50
    
    for (f in seq(0.9, 0.1, -0.1)){
      print(paste('f = ', f))
      
      n_patches <- f * n * n  # number of habitat patches that should remain
      while (sum(h) > n_patches){
        dist <- sqrt((x - x[destroyed])^2 + (y - y[destroyed])^2)
        w    <- w + 1/(2*pi) * ((2 - mu) / (lmax^(2-mu) - lmin^(2-mu))) * dist^(-1*mu)
        w    <- w * h
        w[is.na(w)] <- 0  # focal patch
        
        destroyed    <- sample(1:(n*n), 1, prob = w)
        h[destroyed] <- 0
      }
      
      # write the data to a dataframe and then to a file:
      habitat <- paste(x[h == 1], y[h == 1], sep='-')
      df2 <- data.frame(mu = mu,
                       n_ind = sum(species_loc %in% habitat),
                       hab_cover = f)
      df <- rbind(df, df2)
    }
    mu[mu == 2.01] <- 2
  }
}

df$hab_cover2 <- paste(df$hab_cover*100, ' % Habitat cover')
ggplot(df, aes(x=mu, y=n_ind*2)) + 
  geom_point() + 
  facet_wrap(vars(hab_cover2)) + 
  xlab(expression(paste('Fragmentation level (', mu, ')'))) + 
  ylab('% Individuals left in population') +
  theme_bw() +
  theme(strip.placement = "outside", 
        strip.background = element_blank())
