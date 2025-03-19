# Fragment communities.R

fragment_community <- function(n_ind, Pm_range,clustering, sim_nr, mutation_rate, max_mutation){
  library(ggplot2)
  source('src/simulate_community_dynamics.R')
  source('./src/write_data_to_files.R')
  
  filename  <- paste('./results/landscapes/landscapes_',
                     sim_nr, '_', clustering, '.txt', sep='')
  landscape <- as.matrix(read.table(filename))
  
  rv <- list(n          = nrow(landscape),
             n_ind      = n_ind, 
             Pm_range   = Pm_range,
             nx         = sqrt(nrow(landscape)),
             ny         = sqrt(nrow(landscape)), 
             clustering = clustering,
             mutation_rate = mutation_rate,
             max_mutation = max_mutation) 
  
  # position the cells in the lattice:
  rv$x <- rep(1:rv$nx, rv$ny)
  rv$y <- sort(rep(1:rv$ny, rv$nx))
  
  # create an empty landscape around the forest (at 1-5 cells outside 
  # of the edge subcommunities):

  rv$comm_type                           <- rep('sub', length(rv$x))
    
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
  rv$species                          <- rep(1:rv$n_ind, length(rv$x))
  Pm      <- floor((rv$species / (rv$n_ind + 1)) * 9) / 10 + 0.1
  Pm      <- Pm * rv$Pm_range
  rv$Pm <- Pm

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
    
    filename <- paste('./results/community_composition/',
                      rv$clustering,'_',sim_nr, '_',frag, 
                      '_', rv$mutation_rate,'_', rv$max_mutation,'.txt', sep='')
    if (!file.exists(filename)){

     if ((frag > 0)&(startup == TRUE)){
       filename <- paste('./results/community_composition/',
                         rv$clustering,'_',sim_nr, '_',frag - 5, 
                         '_', rv$mutation_rate,'_', rv$max_mutation,'.txt', sep='')
        a        <- as.vector(as.matrix(read.table(filename)))
        a2        <- as.numeric(unlist(strsplit(a, '-'))[(1:(length(a)))*2 - 1])
        hab      <- landscape[,iFrag-1]
        comm     <- rv$comm_ID[hab == 1]
        rv$species[rv$comm_ID2 %in% comm] <- a2
        rv$Pm[rv$comm_ID2 %in% comm] <- as.numeric(unlist(strsplit(a, '-'))[(1:(length(a)))*2])
     }
      startup <- FALSE
      hab  <- landscape[,iFrag]
      rv$comm_type[hab == 0] <- 'fragmented'
      rv$comm_type2 <- unlist(lapply(rv$comm_type, function(x) 
        rep(x, rv$n_ind)))
      rv$species[rv$comm_type2 == 'fragmented'] <- 0
      rv$Pm[rv$comm_type2 == 'fragmented'] <- 0
      rv$nspecies <- sapply(rv$comm_ID, calc_n_spec)
      rv          <- simulate_community_dynamics(rv)
      rv$frag     <- frag
      rv$sim_nr   <- sim_nr
      write_data_to_files(rv)
    }
  }
}

#fragment_community(1000, 0.5, 5, 1, 0.0001, 3)

# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
n_ind        <- 1000
Pm_range     <- 0.5
clustering   <- c(1, 3, 5)
mutation_rate <- 0.0001 #c(0.00001, 0.0001, 0.001, 0.01)
max_mutation <- 0.05 # c(0, 0.05, 0.1, 0.15, 0.2)
sim_nr       <- 1:10

dat <- expand.grid(n_ind = n_ind, 
                   Pm_range = Pm_range, 
                   clustering = clustering, 
                   mutation_rate = mutation_rate, 
                   max_mutation = max_mutation, 
                   sim_nr = sim_nr)

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:length(dat$n_ind), function(i) {
  fragment_community(dat$n_ind[i], dat$Pm_range[i], 
                     dat$clustering[i], dat$sim_nr[i], dat$mutation_rate[i],
                     dat$max_mutation[i])
}, future.seed = TRUE)

