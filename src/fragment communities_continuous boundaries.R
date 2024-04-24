# Fragment communities.R

fragment_community <- function(n_ind, Pm_range,clustering, sim_nr){
  library(ggplot2)
  source('./src 03-2024/simulate_community_dynamics.R')
  source('./src 03-2024/write_data_to_files.R')
  
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

  startup <- TRUE
  for (iFrag in 1:20){
    frag <- seq(0, 95, 5)[iFrag]
    
    filename <- paste('./results/simulation results 03-2024/community_composition/',
                      rv$clustering,'_',sim_nr, '_',frag, '.txt', sep='')
    if (!file.exists(filename)){

     if ((frag > 0)&(startup == TRUE)){
        filename <- paste('./results/simulation results 03-2024/community_composition/',
                          rv$clustering,'_',sim_nr, '_',frag-5, '.txt', sep='')
        a        <- as.vector(as.matrix(read.table(filename)))
        hab      <- landscape[,iFrag-1]
        comm     <- rv$comm_ID[rv$comm_type != 'meta'][hab == 1]
        rv$species[rv$comm_ID2 %in% comm] <- a
     }
      startup <- FALSE
      hab  <- landscape[,iFrag]
      rv$comm_type[rv$comm_type != 'meta'][!hab] <- 'fragmented'
      rv$comm_type2 <- unlist(lapply(rv$comm_type, function(x) 
        rep(x, rv$n_ind)))
      rv$species[rv$comm_ID2 %in% rv$comm_ID[rv$comm_type == 'fragmented']] <- 0
      rv$nspecies <- sapply(rv$comm_ID, calc_n_spec)
      rv          <- simulate_community_dynamics(rv)
      rv$frag     <- frag
      rv$sim_nr   <- sim_nr
      write_data_to_files(rv)
    }
  }
}

fragment_community(1000, 0.5, 4, 5)

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
  fragment_community(v_n_ind[i], v_Pm_range[i], 
                     v_clustering[i], v_sim_nr[i])
}, future.seed = TRUE)

