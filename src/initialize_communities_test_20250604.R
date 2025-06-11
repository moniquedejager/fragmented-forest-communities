# To examine how starting with a different distribution of dispersal strategies 
# will affect this distribution at the end of the initialization phase, we now 
# performed two additional simulations which start with the following initial 
# distributions of dispersal strategies: (i) a peak around λ = 0.1 (using a 
# logit back-transform from a normal distribution with mean = -1.5 and sd = 1), 
# and (ii) a peak around λ = 0.9 (using a logit back-transform from a normal 
# distribution with mean = 1.5 and sd = 1). 


# use logit function to get the normal distribution within the boundaries of lambda:

x <- rnorm(1000, -1.5, 1)
hist(x)

# logit = log(x/(1-x))
# exp(logit) = x / (1-x)
# exp(logit) - exp(logit)x = x
# exp(logit)x + x = exp(logit)
# x(exp(logit) + 1) = exp(logit)
# x = exp(logit)/(exp(logit) + 1)

y = exp(x)/(exp(x) + 1) * 1000
hist(y, breaks=0:1000 + 0.5)

x = (1:1000) /1000
y = sample(x, 1000, replace=TRUE, prob=x)

initialize_community <- function(n_ind, 
                               Pm_range,
                               disp_distribution,
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
  source('./src/write_data_to_files_initialization_test_20250604.R')
  
  rv <- list(n          = 400,
             n_ind      = n_ind, 
             Pm_range   = Pm_range,
             disp_distribution = disp_distribution,
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
  if (disp_distribution == 1){
    # all species are equally represented
    rv$species <- rep(1:rv$n_ind, rv$n) 
  } else { 
    if (disp_distribution == 2){
      # peak at lambda = 0.1
      x          <- rnorm(length(rv$x2), -1.5, 1)
      rv$species <- round(exp(x)/(exp(x) + 1) * 1000)
    } else {
      # peak at lambda = 0.9
      x          <- rnorm(length(rv$x2), 1.5, 1)
      rv$species <- round(exp(x)/(exp(x) + 1) * 1000)
    }
  }
  rv$Pm      <- ceiling(rv$species/(rv$n_ind/10))/10 * rv$Pm_range
  
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
  
  rv$simulation_type <- 'Initialization'
  rv  <- simulate_community_dynamics(rv)
  
  write_data_to_files_initialization_test_20250604(rv)
}

# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
dat <- expand.grid(n_ind = 1000, 
                   Pm_range = 1, 
                   disp_distribution = 1:3,
                   clustering = 1, 
                   mutation_rate = 0.0003, 
                   max_mutation = 0,  
                   sim_nr = 1:5,
                   f_loss = 0,
                   dispersal = 'different')

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:length(dat$n_ind), function(i) {
  initialize_community(dat$n_ind[i], 
                     dat$Pm_range[i], 
                     dat$disp_distribution[i],
                     dat$clustering[i], 
                     dat$sim_nr[i], 
                     dat$mutation_rate[i],
                     dat$max_mutation[i],
                     dat$f_loss[i],
                     dat$dispersal[i])
}, future.seed = TRUE)

