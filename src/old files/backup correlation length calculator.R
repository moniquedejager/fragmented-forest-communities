sim_nr <- 10
clustering <- 1
Pm_range <- 0.5
n_ind <- 1000

calc_corr_length <- function(sim_nr, clustering, Pm_range, n_ind){
  filename  <- paste('./results/simulation results 03-2024/landscapes/landscapes_',
                           sim_nr, '_', clustering, '.txt', sep='')
  landscape <- read.table(filename)
  
  rv <- list(n          = nrow(landscape),
             n_ind      = n_ind, 
             Pm_range   = Pm_range,
             nx         = sqrt(nrow(landscape)),
             ny         = sqrt(nrow(landscape)), 
             clustering = clustering) 
  
  # position the cells in the lattice:
  rv$x <- rep(1:rv$nx, rv$ny)
  rv$y <- sort(rep(1:rv$ny, rv$nx))
  
  # create the meta community around the forest (at 1-5 cells outside 
  # of the edge subcommunities):
  y_meta                                 <- rep(-4:(rv$ny+5), 
                                                length(-4:(rv$nx+5)))
  x_meta                                 <- sort(rep(-4:(rv$nx+5), 
                                                     length(-4:(rv$ny+5))))
  group_meta                             <- paste(x_meta, y_meta, sep='-')
  group                                  <- paste(rv$x, rv$y, sep='-')
  rv$comm_type                           <- rep('sub', length(x_meta))
  rv$comm_type[!(group_meta %in% group)] <- 'meta'
  
  rv$x      <- x_meta
  rv$y      <- y_meta
  
  rv$x2 <- unlist(lapply(rv$x, function(x) 
    rep(x, rv$n_ind)))
  rv$y2 <- unlist(lapply(rv$y, function(x) 
    rep(x, rv$n_ind)))
  
  rv$comm_ID    <- 1:length(rv$x)
  rv$comm_type2 <- unlist(lapply(rv$comm_type, function(x) 
    rep(x, rv$n_ind)))
  rv$comm_ID2   <- unlist(lapply(rv$comm_ID, function(x) 
    rep(x, rv$n_ind)))
  
  # start with an initial community with 1 individual per species per subcommunity 
  rv$species                         <- rep(1:rv$n_ind, length(rv$x))
  
  # we furthermore need to define the local neighborhood per cell
  rv$dx   <- rep(-5:5, 11)
  rv$dy   <- sort(rv$dx)
  rv$dist <- sqrt(rv$dx^2 + rv$dy^2)
  
  calc_n_spec <- function(x) {
    spec <- rv$species[x*rv$n_ind - (1:rv$n_ind) + 1]
    return(length(unique(spec[spec > 0])))
  } 
  rv$nspecies <- sapply(rv$comm_ID, calc_n_spec)
  rv$tot      <- rv$n_ind * rv$n
  rv$tot2     <- length(rv$x2)
  
  # We fragment the environment with steps of 5% of all patches. Patch 
  # destruction is already done with the r-script 'create_landscapes.R'. 
  # We only have to select the right column from the landscape file to 
  # read in which patches are still habitat and which are destructed. 
  for (iFrag in 1:20){
    frag <- seq(0, 95, 5)[iFrag]
    
    filename <- paste('./results/simulation results 03-2024/community_composition/',
                            rv$clustering,'_',sim_nr, '_',frag, '.txt', sep='')
    
    if (file.exists(filename)){
      a        <- as.vector(as.matrix(read.table(filename)))
    
      hab      <- landscape[,iFrag]
      comm     <- rv$comm_ID[rv$comm_type != 'meta'][hab == 1]
      rv$species[rv$comm_ID2 %in% comm] <- a
      
      rv$comm_type[rv$comm_type != 'meta'][!hab] <- 'fragmented'
      rv$comm_type2 <- unlist(lapply(rv$comm_type, function(x) 
            rep(x, rv$n_ind)))
      rv$species[rv$comm_ID2 %in% rv$comm_ID[rv$comm_type == 'fragmented']] <- 0
      rv$nspecies <- sapply(rv$comm_ID, calc_n_spec)
      rv$frag     <- frag
      rv$sim_nr   <- sim_nr
          
      # calculate spatial correlation length per dispersal capacity
      C    <- vector(length = 0)
      v_Pm <- unique((floor((rv$species / (rv$n_ind + 1)) * 9) / 10 + 0.1) * rv$Pm_range) 
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
        r_commID <- (x_comm - min(rv$x))*(rv$ny+10) + (y_comm - min(rv$y)) + 1
        
        # only replace individuals in the subcommunities:
        sel      <- r_commID %in% rv$comm_ID[rv$comm_type == 'sub']
        
        # and in a random order (otherwise, the last individual always replaces 
        # an earlier dispersing individual)
        a <- sample(1:sum(sel), sum(sel), replace=F)
        
        spec     <- rv$species
        spec[(r_commID[sel] - 1) * rv$n_ind + r_ind[sel]][a] <- rv$species[sel][a]
        
        origin_ID <- sort(rep(rv$comm_ID, rv$n_ind)) 
        origin_ID[(r_commID[sel] - 1) * rv$n_ind + r_ind[sel]][a] <- 
          origin_ID[sel][a]
        
        comm_ID   <- rv$comm_ID2[(rv$comm_type2 == 'sub')]
        origin_ID <- origin_ID[(rv$comm_type2 == 'sub')]
        u_ID      <- sort(unique(c(comm_ID, origin_ID)))
        
        x   <- rv$x[rv$comm_ID %in% u_ID]
        y   <- rv$y[rv$comm_ID %in% u_ID]
        IDs <- rv$comm_ID[rv$comm_ID %in% u_ID]
        n   <- length(x)
        
        connected <- matrix(0, n, n)
        for (i in IDs){
          j = unique(origin_ID[comm_ID %in% i])
          connected[(IDs == i),(IDs %in% j)] <- 1
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
      
      df <- data.frame(mu         = rv$clustering,
                       f_hab_loss = rv$frag,
                       sim_nr     = rv$sim_nr,
                       dispersal_capacity = unique((floor((rv$species / 
                                                             (rv$n_ind + 1)) * 9) / 
                                                      10 + 0.1) * rv$Pm_range),
                       corr_length= C)
      filename <- paste('./results/simulation results 03-2024/correlation_length/correlation_length_data_',
                        rv$sim_nr, '.txt', sep='')
      if (file.exists(filename)){
        write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
      } else {
        write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
      }
    }
  }
}

# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
n_ind        <- 1000
Pm_range     <- 0.5
clustering   <- 1:4

v_n_ind      <- rep(n_ind, 40)
v_Pm_range   <- rep(Pm_range, 40)
v_clustering <- sort(rep(clustering, 10))
v_sim_nr     <- rep(1:10, 4)

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:length(v_n_ind), function(i) {
  calc_corr_length(v_sim_nr[i], v_clustering[i], v_Pm_range[i], v_n_ind[i])
})

