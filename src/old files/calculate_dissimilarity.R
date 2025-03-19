calc_dissimilarity <- function(rv){
  # Calculate the Bray-Curtis dissimilarity between subcommunities
  df2 <- matrix(rv$species, rv$n_ind,)
  
  # select 20 random subcommunities:
  sel <- rv$comm_ID %in% sample(rv$comm_ID[rv$comm_type == 'sub'], 20, 
                                 replace = FALSE)
  
  species   <- df2[,sel]
  
  # Euclidian distance between the subpopulations:
  x <- rv$x[sel]
  y <- rv$y[sel]
  
  n <- length(x)
  i  <- rep(1, (n-1))
  j  <- 2:n
  for (ix in 2:(n-1)){
    i <- c(i, rep(ix, n-ix))
    j <- c(j, (ix+1):n)
  }
  
  d <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
  
  dissimilarity <-  sapply(1:length(i), function(ix){
    
    spec1 <- species[,i[ix]]
    spec2 <- species[,j[ix]]
    
    max_spec <- max(spec1, spec2)
    
    h1 <- hist(spec1, breaks=0:max_spec + 0.5, plot=FALSE)
    h2 <- hist(spec2, breaks=0:max_spec + 0.5, plot=FALSE)
    
    # Bray-Curtis dissimilarity (0..1, 0 = equal, 1 = different)
    dissimilarity <- bray_curtis(h1$counts, h2$counts)
    return(dissimilarity)
  })
  
  return(cbind(d, dissimilarity))
}
