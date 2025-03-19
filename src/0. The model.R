# to change: 
# - max_mutation doet niets, deze mag eruit
# - fragment.R moet anders
# - welke data willen we wegschrijven? Volgens mij alle individuen en hun locatie?
# - er is nog geen start community, dus dat moet ook anders
# - omschrijven zodat we per landschapsgrootte de berekensnelheid kunnen achterhalen. 

# hoe groter het model, hoe langzamer hij wordt (heel logisch natuurlijk). Is
# het een optie om het model toch maar om te schrijven naar C++? 

fragment_community <- function(n_ind, 
                              Pm_range,
                              clustering, 
                              sim_nr, 
                              mutation_rate, 
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
  source('src/1. simulate_community_dynamics.R')
  #source('src/fragment.R')
  source('./src/2. write_species_to_files.R')
  
  for (nx in (1:10)*10){
    start_time1   <- Sys.time()
    rv <- list(n          = nx*nx,
               n_ind      = n_ind, 
               Pm_range   = Pm_range,
               nx         = nx,
               ny         = nx,
               sim_nr     = sim_nr,
               clustering = clustering,
               mutation_rate = mutation_rate,
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
    
    rv$species <- rep(1:rv$n_ind, rv$n);
    rv$Pm      <- ceiling(rv$species/(rv$n_ind/10))/10 * rv$Pm_range
    
    # start with the pristine, unfragmented landscape, with the initial community:
    #if (dispersal == 'similar'){
    #  m <- read.table('results/community_composition/initial_community_sim.txt')
    #  m <- as.vector(t(as.matrix(m)))
    #  rv$species <- m
    #  rv$Pm      <- rep(0.8, rv$n*rv$n_ind)
    #} else {
    #  m <- read.table('results/community_composition/initial_community_dif.txt')
    #  m <- as.vector(t(as.matrix(m)))
    #  m2  <- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2 - 1])
    #  rv$Pm     <- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2])
    #  rv$species <- m2
    #}
    
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
    
    write_species_to_files(rv)
    end_time1 <- Sys.time()
    print(end_time1 - start_time1)
  }
}

fragment_community(1000, 1, 1, 0.0003, 1, 0, 'different')

# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
dat <- expand.grid(n_ind = 1000, 
                   Pm_range = 1, 
                   clustering = c(1,3,5), 
                   mutation_rate = 0.0003, 
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
                    dat$f_loss[i],
                    dat$dispersal[i])
}, future.seed = TRUE)


df <- data.frame(x = 1:10*10,
                 landscape_size  = 0.1*0.25*(1:10*10)^2, # in km^2
                 simulation_time = c(4.7, 20.1, 43.3, 1.28*60,
                                     2.03*60, 3.01*60, 4.09*60, 5.48*60,
                                     7.72*60, 10.50*60),
                 sim_time_cpp = c(0.8, 3.4, 7.6, 15.4,
                                  20.6, 30.0, 40.8, 67.9, 95.3, 113.4)) # per 100 iterations

df2 <- data.frame(landscape_size = rep(df$landscape_size, 2), 
                  simulation_time = c(df$simulation_time, df$sim_time_cpp),
                  code = c(rep("R", nrow(df)), rep("C++", nrow(df))))
library(ggplot2)

ggplot(df2, aes(x=landscape_size, y=simulation_time/100, color=code)) + 
  geom_point() + 
  geom_line() + 
  xlab('Landscape size (km2)') + 
  ylab('Simulation time (seconds per iteration)') + 
  theme_bw()

df2$total_sim_time = df2$simulation_time*1.67/24

ggplot(df2, aes(x=landscape_size, y=simulation_time*1.67/24, color=code)) + 
  geom_point(aes(shape=code), size=3) + 
  geom_line(aes(linetype=code), linewidth=1.2) + 
  xlab('Landscape size (km2)') + 
  ylab('Simulation time (days per simulation)') + 
  scale_color_discrete(name='') + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  theme_bw() + 
  theme(legend.position = 'top')

df2
# 4x4


# dus, als we een gebied van ongeveer 50 km2 willen simuleren, dan hebben we
# 45 x 45 = 2025 patches nodig. In C++ zou dat dan ongeveer 1.2 dagen duren, 
# en in R 7 dagen. 

# 23*1.67/24 = 1.6 dagen
# 207 sec. voor 1000 iteraties



