# As an example, we will create a smaller, 1-D environment, consisting of n_cell 
# cells. f_hab % of these cells are habitable and contain n_indiv individuals. 
# There are 4 different dispersal strategies: 
# (1) local dispersal only, 
# (2) 80% local and 20% to an adjacent cell, 
# (3) 60% local, 30% to an adjacent cell, and 10% to a cell one cell further, 
# and (4) equal chances of dispersing to a cell at distance 0, 1, and 2. 
# Each individual disperses n_seeds seeds. Only one seed can replace an 
# individual; which seed will replace the individual (if any) will be randomly 
# selected. 

# assumptions: all species produce an equal amount of seeds. There is no 
# difference in the cost of producing long-distance or short-distance dispersing
# seeds. 

seed_dispersal <- function(par){
  n_cell  <- par$n_cell
  n_indiv <- par$n_indiv
  n_seeds <- par$n_seeds
  f_hab   <- par$f_hab
  simnr   <- par$simnr
  
  n_gens  <- 1000
  
  habitat <- sample(1:n_cell, f_hab*n_cell, replace = FALSE)
  individuals <- data.frame(x = sort(rep(1:n_cell, n_indiv)),
                            y = rep(1:n_indiv, n_cell),
                            strategy = rep(1:4, n_cell*n_indiv / 4))
  
  individuals <- individuals[individuals$x %in% habitat,]
  
  probs <- matrix(c(1, 0, 0, 
                    0.8, 0.2, 0,
                    0.6, 0.3, 0.1,
                    0.33, 0.33, 0.33), 3, 4)
  
  k <- 0
  while ((k < n_gens)&(length(unique(individuals$strategy)) > 1)){
    print(k)
    k <- k + 1
    
    seeds <- data.frame(x = vector(length=0),
                        y = vector(length=0),
                        strategy = vector(length=0))
    
    for (i in 1:nrow(individuals)){
      print(i)
      dx  <- sample(0:2, n_seeds, replace = TRUE, prob=probs[,individuals$strategy[i]])
      dir <- sample(c(-1, 1), n_seeds, replace = TRUE)
      y   <- sample(1:n_indiv, n_seeds, replace = TRUE)
      seeds2 <- data.frame(x = individuals$x[i] + dx*dir,
                           y = y,
                           strategy = individuals$strategy[i])
      seeds <- rbind(seeds, seeds2)
    }
    seeds$x[seeds$x < 1] <- n_cell + seeds$x[seeds$x < 1]
    seeds$x[seeds$x > n_cell] <- seeds$x[seeds$x > n_cell] - n_cell 
    
    rand_order <- sample(1:nrow(seeds), nrow(seeds), replace = FALSE)
    
    for (i in rand_order){
      individuals[(individuals$x == seeds$x[i])&(individuals$y == seeds$y[i]),] <- seeds[i,]
    }
    individuals <- individuals[individuals$x %in% habitat,]
  }
  
  hist(individuals$strategy, breaks=0:4+0.5)
  
  data          <- par
  data$strategy <- mean(individuals$strategy) 
  data$generation <- k
  
  write.table(data, 'results/example of seed dispersal.txt', row.names = FALSE,
              col.names = FALSE, append = TRUE)
}


# Install and load the future package
# install.packages("future")
library(future)
library(future.apply)

# Create input vectors/lists
params <- expand.grid(n_cell  = c(10, 50, 100),
                      n_indiv = c(4, 40, 400),
                      n_seeds = c(1, 10, 100),
                      f_hab   = c(1, 0.9, 0.7),
                      simnr   = 1:5)

# discard the parameter combinations that lead to too great simulations to run:
param_size <- params$n_cell*params$n_indiv*params$n_seeds
params <- params[param_size < 40000,]

# Set up parallel processing with future
plan(multisession, workers = 10)  # Adjust the number of workers based on your system

result_parallel <- future.apply::future_lapply(1:nrow(params), function(i) {
  seed_dispersal(params[i,])
}, future.seed = TRUE)
