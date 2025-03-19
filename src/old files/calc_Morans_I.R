calc_Morans_I <- function(rv){
  # calculate Moran's I, mean patch size, and number of patches
  
  # libraries:
  library(GGally)    
  source('./src/calculate_patch_size.R')
  
  n <- rv$nx
  x <- rep(1:n, n)
  y <- sort(x)
  h <- rv$comm_type == 'sub'
  
  patch_size <- calc_patch_size(x, y, h)
          
  # create vectors of the habitat values of each neighbor:
  zmean      <- mean(h)
  d0         <- h - zmean
  
  x2         <- x - 1
  x2[x2 < 1] <- n
  ID2        <- (y - 1)*n + x2
  n1         <- h[ID2]
  d1         <- n1 - zmean
  
  x2         <- x + 1
  x2[x2 > n] <- 1
  ID2        <- (y - 1)*n + x2
  n2         <- h[ID2]
  d2         <- n2 - zmean
  
  y2         <- y - 1
  y2[y2 < 1] <- n
  ID2        <- (y2 - 1)*n + x
  n3         <- h[ID2]
  d3         <- n3 - zmean
  
  y2         <- y + 1
  y2[y2 > n] <- 1
  ID2        <- (y2 - 1)*n + x
  n4         <- h[ID2]
  d4         <- n4 - zmean
  
  # Moran's I = 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2))  
  I <- 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2)) 
  
  df <- data.frame( morans_I     = I,
                    m_patch_size = mean(patch_size),
                    n_patches    = length(patch_size))
  return(df)
}