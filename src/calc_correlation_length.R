calc_correlation_length <- function(rv){
  # calculate spatial correlation length per dispersal capacity
  C    <- vector(length = 0)
  v_Pm <- unique(rv$Pm) 
  for (jPm in v_Pm){
    # what if all individuals dispersed using dispersal strategy Pm?
    Pm <- rep(jPm, length(rv$species))
    Pm[rv$species == 0] <- 0
    
    # Per Pm value, select all individuals with this value and randomly 
    # select the subcommunity to reproduce and disperse to and the random
    # individual to replace in this subcommunity
    r_comm <- rep(61, rv$tot2) # 61 is the focal subcommunity
    r_ind  <- sample(1:rv$n_ind, rv$tot2, replace=TRUE)
    for (iPm in unique(Pm[Pm > 0])) {
      sel <- Pm == iPm
      Pm2 <- (2*pi*iPm^2)^-1 * exp(-1*rv$dist / iPm) 
      Pm2 <- Pm2 / sum(Pm2)
      r_comm[sel] <- sample(1:121, sum(sel), replace=TRUE, prob=Pm2)
    }
    x_comm   <- rv$x2 + rv$dx[r_comm]
    y_comm   <- rv$y2 + rv$dy[r_comm]

    x_comm <- x_comm - rv$nx*(x_comm > rv$nx) + rv$nx*(x_comm < 1)
    y_comm <- y_comm - rv$ny*(y_comm > rv$ny) + rv$ny*(y_comm < 1)
    r_commID <- (x_comm - 1)*rv$ny + y_comm
    
    # and in a random order (otherwise, the last individual always replaces 
    # an earlier dispersing individual)
    a <- sample((1:rv$tot2), rv$tot2, replace=F)
    
    ind_id   <- ((r_commID - 1) * rv$n_ind + r_ind)      
    
    origin_ID <- sort(rep(rv$comm_ID, rv$n_ind)) 
    origin_ID2 <- origin_ID
    origin_ID[ind_id[a]] <- 
      origin_ID2[a]
    
    comm_ID   <- rv$comm_ID2[(rv$comm_type2 == 'sub')]
    origin_ID <- origin_ID[(rv$comm_type2 == 'sub')]
    
    notSub <- rv$comm_ID[rv$comm_type != 'sub'] 
    origin_ID[origin_ID %in% notSub] <- comm_ID[origin_ID %in% notSub]
    u_ID      <- sort(unique(c(comm_ID, origin_ID)))
    
    x   <- rv$x[rv$comm_ID %in% u_ID]
    y   <- rv$y[rv$comm_ID %in% u_ID]
    IDs <- rv$comm_ID[rv$comm_ID %in% u_ID]
    n   <- length(x)
    
    connected <- matrix(0, n, n)
    for (i in IDs){
      j = unique(origin_ID[comm_ID %in% i])
      connected[(IDs == i),(IDs %in% j)] <- 1
      connected[(IDs == i),(IDs == i)] <- 1
    }
    
    patch_nr <- 1:n
    together <- connected*patch_nr
    group_nr <- rep(0, n)
    
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
    cluster_nr <- group_nr
    
    R <- vector(length = 0)
    n <- vector(length = 0)
    for (i in unique(cluster_nr))
    {
      x2 <- x[cluster_nr == i]
      y2 <- y[cluster_nr == i]
      R <- c(R, mean(sqrt((x2 - mean(x2))^2 + (y2 - mean(y2))^2)))
      n <- c(n, length(x))
    }
    C <- c(C, sum(n * R)/sum(n))
  }
  return(C)
}
