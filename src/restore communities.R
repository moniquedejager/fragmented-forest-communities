# Restore communities.R
# this script is not finished! It is a copy of fragment communities.R

# nog in te verwerken: tijd tussen fragmentatie en restoratie
# wegschrijven data per 50 iteraties

n_ind = 1000
Pm_range = 0.5
clustering = 1
sim_nr = 4
mutation_rate = 0 #0.0001
max_mutation = 0.05
f_loss = 0.95
hab_cover = 0.95
clustering_restored = 1
n_iterations = 50

restore_community <- function(n_ind, 
                              Pm_range,
                              clustering, 
                              sim_nr, 
                              mutation_rate, 
                              max_mutation,
                              f_loss,
                              hab_cover,
                              clustering_restored,
                              n_iterations){
  library(ggplot2)
  source('src/simulate_community_dynamics.R')
  source('src/fragment.R')
  source('src/restore.R')
  source('./src/write_data_to_files_restoration.R')
  
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
             hab_cover     = hab_cover,
             clustering_restored = clustering_restored,
             n_iterations  = n_iterations) 
  
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
  
  # start with the pristine, unfragmented landscape, with 1 individual of 
  # species 1 to 1000 per subcommunity:
  rv$species <- rep(1:rv$n_ind, rv$n)
  rv$Pm      <- ceiling(rv$species/100)/10 * rv$Pm_range
  
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
  
  rv$simulation_type <- 'Initializing'
  rv  <- simulate_community_dynamics(rv)

  # now that we have a starting community, we can fragment the environment:
  rv  <- fragment(rv)
  rv$simulation_type <- 'Fragmentation'
  rv  <- simulate_community_dynamics(rv)
  
  # restore some of the communities:
  rv <- restore(rv)
  rv$simulation_type <- 'Restoration'
  rv  <- simulate_community_dynamics(rv)
  
  #write_data_to_files_restoration(rv)
}




# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
dat <- expand.grid(n_ind = 1000, 
                   Pm_range = 0.5, 
                   clustering = c(1,3,5), 
                   mutation_rate = 0.0001, 
                   max_mutation = 0.05, 
                   sim_nr = 1,
                   f_loss = round(seq(0.05, 0.95, 0.05), 2),
                   hab_cover = round(seq(0.05, 0.95, 0.05), 2),
                   clustering_restored = c(1, 3, 5))

dat <- dat[round((1 - dat$f_loss),2) <= dat$hab_cover,]

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:length(dat$n_ind), function(i) {
  restore_community(dat$n_ind[i], 
                    dat$Pm_range[i], 
                    dat$clustering[i], 
                    dat$sim_nr[i], 
                    dat$mutation_rate[i],
                    dat$max_mutation[i],
                    dat$f_loss[i],
                    dat$hab_cover[i],
                    dat$clustering_restored[i])
}, future.seed = TRUE)

