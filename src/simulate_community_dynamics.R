
simulate_community_dynamics <- function(rv){
  ###################################################################
  # second, replace individuals in the local communities until 
  # the number of species stabilizes:
  ###################################################################
    
  source('src/record.R')
  print(table(rv$Pm))
  
  calc_n_spec <- function(x) {
    spec <- rv$species[x*rv$n_ind - (1:rv$n_ind) + 1]
    length(unique(spec[spec > 0]))
  } 

  rv$static_n_species <- length(unique(rv$species))
  mean_nspecies   <- vector(length=0)
  total_species   <- vector(length=0)
  rv$iteration_nr <- 0

  while (rv$iteration_nr < 501){
    start_time   <- Sys.time()
    rv$origin_ID_t50 <- rv$comm_ID2 
    for (i in 1:50){
      record(rv)
      # print(table(rv$Pm))
      
      # we let each individual reproduce and disperse to a subcommunity:
      # dispersal depends on the dispersal parameter, which corresponds to 
      # the species number. To be able to code this efficiently, we make 
      # 10 dispersal parameter values instead of the entire range (0.1 - 0.9):
      Pm <- rv$Pm
      
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
      r_commID <- (y_comm - 1)*rv$nx + x_comm
      
      # and in a random order (otherwise, the last individual always replaces 
      # an earlier dispersing individual)
      sel <- (rv$species > 0)
      a <- sample((1:rv$tot2)[sel], sum(sel), replace=F)

      spec     <- rv$species
      ind_id   <- ((r_commID - 1) * rv$n_ind + r_ind)      
      spec[ind_id[a]] <- rv$species[a]
      Pm       <- rv$Pm
      Pm[ind_id[a]] <- rv$Pm[a]
      rv$species <- spec
      rv$Pm      <- Pm
      
      # add mutations: (which are now incoming seeds from the metapopulation!)
      mutated                  <- runif(length(rv$species), 0, 1) < rv$mutation_rate
      mutated[rv$species == 0] <- 0
      if (sum(mutated) > 0){
        new_species         <- sample(1:rv$n_ind, 
                                      sum(mutated), 
                                      replace=TRUE, 
                                      prob=(ceiling((1:rv$n_ind)/(rv$n_ind/10))/10 * rv$Pm_range))
        rv$species[mutated==1] <- new_species
        
        if (rv$dispersal == 'different'){
          new_Pm              <- ceiling(new_species/(rv$n_ind/10))/10 * rv$Pm_range
          rv$Pm[mutated==1]      <- new_Pm
        }
      }
      rv$species[rv$comm_type2 != 'sub'] <- 0
      rv$Pm[rv$comm_type2 != 'sub'] <- 0
      
      origin_ID_t50 <- rv$origin_ID_t50
      rv$origin_ID_t50[ind_id[a]] <- 
        origin_ID_t50[a]
      
      origin_ID_t1 <- rv$comm_ID2
      rv$origin_ID_t1[ind_id[a]] <- 
        origin_ID_t1[a]
      
      rv$nspecies     <- sapply(rv$comm_ID, calc_n_spec)
      rv$iteration_nr <- rv$iteration_nr + 1
      
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    
    record(rv)
    
    mn            <- mean(rv$nspecies[rv$comm_type =='sub'])
    mean_nspecies <- c(mean_nspecies, mn)
    total_species <- c(total_species, length(unique(rv$species)))
    
    #if (length(mean_nspecies) %in% ((1:10)*5)){
      #record(rv)
    #}
  }
  return(rv)
} 
