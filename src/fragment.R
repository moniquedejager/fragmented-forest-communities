# Patch destruction using weighted sampling. 

fragment <- function(rv){

  # first destroyed patch:
  destroyed    <- sample(1:rv$n, 1)
  rv$comm_type[destroyed] <- 'fragmented'
  
  # weights for sampling patch destruction:
  w <- rep(0, rv$n)
  
  # 2D pareto distribution of weights around the destroyed patch:
  lmin <- 1
  lmax <- 50
  
  n_patches <- (1 - rv$f_loss) * rv$n  # number of habitat patches that should remain
  while (length(rv$comm_type[rv$comm_type == 'sub']) > n_patches){
    dist <- sqrt((rv$x - rv$x[destroyed])^2 + (rv$y - rv$y[destroyed])^2)
    w    <- w + 1/(2*pi) * 
      ((2 - rv$clustering) / (lmax^(2-rv$clustering) - 
                                lmin^(2-rv$clustering))) * 
      dist^(-1*rv$clustering)
    w[rv$comm_type == 'fragmented'] <- 0  # focal patch
    
    destroyed    <- sample(1:rv$n, 1, prob = w)
    rv$comm_type[destroyed] <- 'fragmented'
  }
  
  rv$comm_type2 <- unlist(lapply(rv$comm_type, function(x) 
    rep(x, rv$n_ind)))
  rv$species[rv$comm_type2 == 'fragmented'] <- 0
  rv$Pm[rv$comm_type2 == 'fragmented'] <- 0
  
  return(rv)
}


    

