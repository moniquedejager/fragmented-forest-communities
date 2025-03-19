# Patch destruction using weighted sampling. 

for (sim_nr in 1:10){
  for (mu in seq(1, 5, 2)){
    
    print(paste('mu = ', mu))
    
    # start with a uniform, continuous habitat and fragment it...
    if (mu == 2){ mu <- 2.01}
    n <- 20
    x <- rep(1:n, n)
    y <- sort(x)
    h <- rep(1, n*n)
    
    mat_h <- h
    
    # first destroyed patch:
    destroyed    <- sample(1:(n*n), 1)
    h[destroyed] <- 0
    
    # weights for sampling patch destruction:
    w <- rep(0, n*n)
    
    # 2D pareto distribution of weights around the destroyed patch:
    lmin <- 1
    lmax <- 50
    
    for (f in seq(0.95, 0.05, -0.05)){
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
      mat_h <- cbind(mat_h, h)
      
    }
    mu[mu == 2.01] <- 2
    filename <- paste('./results/simulation results 03-2024/landscapes/landscapes_', sim_nr, '_', mu, '.txt', sep='')
    write.table(mat_h, filename, append=FALSE, col.names = FALSE, row.names = FALSE)
  }
}
