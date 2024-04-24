calc_patch_size <- function(x, y, h){
  # function to calculate patch sizes with. 
  # Input: x and y coordinates of all cells, 
  # and whether they are habitat cells or not (h = 1 / h = 0)
  
  mat_x <- matrix(rep(x[h == 1], length(x[h==1])), 
                  length(x[h==1]), length(x[h==1]))
  mat_y <- matrix(rep(y[h == 1], length(x[h==1])), 
                  length(x[h==1]), length(x[h==1]))
  
  dist     <- sqrt((mat_x - t(mat_x))^2 + (mat_y - t(mat_y))^2)
  patch_nr <- 1:ncol(mat_x)
  together <- (dist <= 1)*patch_nr
  
  group_nr <- rep(0, ncol(mat_x))
  
  nr <- 0
  while (sum(group_nr == 0) > 0){
    nr <- nr + 1
    group_members <- patch_nr[group_nr == 0][1]
    l             <- 0
    
    while (length(group_members) > l){
      l <- length(group_members)
      group_members <- together[,group_members]
      group_members <- unique(group_members[group_members > 0])
    }
    
    group_nr[group_members] <- nr
  }
  patch_size <- tapply(group_nr, group_nr, length)
  return(patch_size)
}