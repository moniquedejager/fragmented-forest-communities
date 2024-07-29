fragment_community <- function(n_ind, 
                              Pm_range,
                              clustering, 
                              sim_nr, 
                              mutation_rate, 
                              max_mutation,
                              f_loss,
                              dispersal){
  # We simulated subcommunity dynamics in a 2-dimensional, semi-spatial, 
  # near-neutral, individual-based model. The environment consists of 20 x 20 
  # cells, each cell having a subcommunity of 1,000 individuals. Initially, 
  # all cells are habitable (i.e. there is no fragmentation) and, contain a 
  # subcommunity. Each subcommunity that during an initialization phase is 
  # populated with the highest diversity theoretically possible during an 
  # initialization phase and is subsequently allowed to stabilize, after which 
  # we run the simulations with fragmentationhabitat loss. We simulated 20 
  # different levels of habitat loss, from 0 to 95% habitat destruction, at 
  # three fragmentation configurations: (i) clustered, (ii) fractal, and (iii) 
  # random habitat destruction. For full details on the model, see De Jager et
  # al. (in prep.). 
  
  # Libraries and custom functions:
  library(ggplot2)
  source('src/simulate_community_dynamics.R')
  source('src/fragment.R')
  source('./src/write_data_to_files_fragmentation.R')
  
  rv <- list(n          = 400,
             n_ind      = n_ind, 
             Pm_range   = Pm_range,
             nx         = 20,
             ny         = 20,
             sim_nr     = sim_nr,
             clustering = clustering,
             mutation_rate = mutation_rate,
             max_mutation  = max_mutation,
             f_loss        = f_loss,
             dispersal     = dispersal)
             
  # position the cells in the lattice:
  rv$x <- rep(1:rv$nx, rv$ny)
  rv$y <- sort(rep(1:rv$ny, rv$nx))
  rv$comm_type <- rep('sub', length(rv$x))
    
  rv$x2 <- unlist(lapply(rv$x, function(x) 
    rep(x, rv$n_ind)))
  rv$y2 <- unlist(lapply(rv$y, function(x) 
    rep(x, rv$n_ind)))
  
  rv$comm_ID    <- 1:length(rv$x)
  rv$comm_type2 <- unlist(lapply(rv$comm_type, function(x) 
    rep(x, rv$n_ind)))
  rv$comm_ID2   <- unlist(lapply(rv$comm_ID, function(x) 
    rep(x, rv$n_ind)))
  
  # start with the pristine, unfragmented landscape, with the initial community:
  if (dispersal == 'similar'){
    m <- read.table('results/community_composition/initial_community_sim.txt')
    m <- as.vector(t(as.matrix(m)))
    rv$species <- m
    rv$Pm      <- rep(0.8, rv$n*rv$n_ind)
  } else {
    m <- read.table('results/community_composition/initial_community_dif.txt')
    m <- as.vector(t(as.matrix(m)))
    m2  <- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2 - 1])
    rv$Pm     <- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2])
    rv$species <- m2
  }
  
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
  
  # now that we have a starting community, we can fragment the environment:
  if (rv$f_loss > 0){
    rv  <- fragment(rv)
  }
  rv$simulation_type <- 'Fragmentation'
  #rv$simulation_type <- 'Initialization'
  rv  <- simulate_community_dynamics(rv)
  
  write_data_to_files_fragmentation(rv)
}

# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
dat <- expand.grid(n_ind = 1000, 
                   Pm_range = 1, 
                   clustering = c(1,3,5), 
                   mutation_rate = 0.0003, 
                   max_mutation = 0,  
                   sim_nr = 1:10,
                   f_loss = round(seq(0, 0.95, 0.05), 2),
                   dispersal = c('similar', 'different'))

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:length(dat$n_ind), function(i) {
  fragment_community(dat$n_ind[i], 
                    dat$Pm_range[i], 
                    dat$clustering[i], 
                    dat$sim_nr[i], 
                    dat$mutation_rate[i],
                    dat$max_mutation[i],
                    dat$f_loss[i],
                    dat$dispersal[i])
}, future.seed = TRUE)

