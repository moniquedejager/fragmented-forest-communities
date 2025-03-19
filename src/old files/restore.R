restore <- function(rv){
  w <- rep(0, rv$n)
  
  # 2D pareto distribution of weights around the destroyed patch:
  lmin <- 1
  lmax <- 50
  
  for (j in rv$comm_ID[rv$comm_type == 'sub']){
    dist <- sqrt((rv$x - rv$x[j])^2 + (rv$y - rv$y[j])^2)
    w    <- w + 1/(2*pi) * 
      ((2 - rv$clustering_restored) / 
         (lmax^(2-rv$clustering_restored) - 
            lmin^(2-rv$clustering_restored))) * 
      dist^(-1*rv$clustering_restored)
    w[is.na(w)] <- 0  # focal patch
  }
  w[rv$comm_type == 'sub'] <- 0
  
  n_patches <- rv$hab_cover * rv$n  # number of habitat patches that should remain
  n_restore <- n_patches - length(rv$comm_type[rv$comm_type == 'sub']) # number of patches that needs to be restored
  restored    <- sample(1:rv$n, n_restore, prob = w)
    
  rv$comm_type[restored] <- 'sub'
  rv$comm_type2[rv$comm_ID2 %in% restored] <- 'sub'
  
  return(rv)
}
