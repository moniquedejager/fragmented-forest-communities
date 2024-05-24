
simulate_community_dynamics <- function(rv){
  ###################################################################
  # second, replace individuals in the local communities until 
  # the number of species stabilizes:
  ###################################################################
    
  # every Pm-value needs a constant to correct for the surface area under the
  # graph: 
  subcom <- 1
  dx <- abs(rv$x - rv$x[rv$comm_ID == subcom])
  dy <- abs(rv$y - rv$y[rv$comm_ID == subcom])
  # continuous boundaries:
  dx[dx > rv$nx/2] <- rv$nx - dx[dx > rv$nx/2]
  dy[dy > rv$ny/2] <- rv$ny - dy[dy > rv$ny/2]
  dist <- sqrt(dx^2 + dy^2)
  sel  <- (dist < 5)
  sum_Pm <- vector(length = 0)
  for (Pm in seq(0.05, 1, 0.05)){
    Pm_weight  <- (2*pi*Pm^2)^-1 * exp(-1*dist[sel]/Pm) 
    sum_Pm <- c(sum_Pm, sum(Pm_weight))
  }
  rv$Pm_C <- 1/sum_Pm
  
  print(table(rv$Pm))
  
  calc_n_spec <- function(x) {
    spec <- rv$species[x*rv$n_ind - (1:rv$n_ind) + 1]
    length(unique(spec[spec > 0]))
  } 
  
  # per cell, select one of the subcommunities (or the metacommunity) 
  # closest to this distance
  mean_nspecies        <- vector(length=0)
  rv$iteration_nr <- 0
  p               <- 0
  while (p < 50){
    start_time   <- Sys.time()
    rv$origin_ID_t50 <- rv$comm_ID2 
    for (i in 1:50){
      #print(round(table(rv$Pm)/(rv$n*rv$n_ind), 2))
      #print(paste('Working on iteration ', i))
      #df <- data.frame(x=rv$x,
      #                 y=rv$y,
      #                 nspec = rv$nspecies,
      #                 mPm   = tapply(rv$Pm, rv$comm_ID2, mean))
      #print(ggplot(df, aes(x=x, y=y, fill=nspec)) + geom_raster())
      
      # for each subcommunity, we will perform a weighted random sampling of
      # all individuals to replace the old generation. Weights are based on the
      # dispersal capacity of the individuals and their distance to the subcom-
      # munity. 
      id_replacement <- sapply(rv$comm_ID[rv$comm_type == 'sub'], 
                               function(subcom){
                                 dx <- abs(rv$x - rv$x[rv$comm_ID == subcom])
                                 dy <- abs(rv$y - rv$y[rv$comm_ID == subcom])
                                 # continuous boundaries:
                                 dx[dx > rv$nx/2] <- rv$nx - dx[dx > rv$nx/2]
                                 dy[dy > rv$ny/2] <- rv$ny - dy[dy > rv$ny/2]
                                 dist <- sqrt(dx^2 + dy^2)
                                 dist <- unlist(lapply(dist, function(x) 
                                   rep(x, rv$n_ind)))
                                 
                                 sel <- (rv$species > 0)&(dist < 5)
                                 ids <- (1:length(rv$x2))[sel]
                                 Pm  <- (2*pi*rv$Pm[sel]^2)^-1 * 
                                   exp(-1*dist[sel]/rv$Pm[sel])*rv$Pm_C[rv$Pm[sel]/0.05] 
                                 
                                 id_replacement <- sample(ids, 
                                                          rv$n_ind, 
                                                          prob = Pm, 
                                                          replace = TRUE)
                                 return(id_replacement)
                               })
      id_replacement_all <- as.vector(id_replacement)
      
      new_species  <- rv$species[id_replacement_all]
      new_Pm       <- rv$Pm[id_replacement_all]
      origin_ID_t1 <- rv$comm_ID2[id_replacement_all]
      origin_ID_t50<- rv$origin_ID_t50[id_replacement_all]
      
      rv$species[rv$comm_type2 == 'sub']       <- new_species
      rv$Pm[rv$comm_type2 == 'sub']            <- new_Pm
      rv$origin_ID_t1[rv$comm_type2 == 'sub']  <- origin_ID_t1
      rv$origin_ID_t50[rv$comm_type2 == 'sub'] <- origin_ID_t50
 
      # add mutations:
      mutated             <- (runif(length(rv$species), 0, 1) < rv$mutation_rate)*1
      mutated[rv$species == 0] <- 0
      if (sum(mutated) > 0){
        new_species         <- max(rv$species) + 1:sum(mutated)
        #new_Pm              <- rv$Pm[mutated == 1] + sample(seq((-1*rv$max_mutation),rv$max_mutation, 0.05), sum(mutated), replace=TRUE) 
        new_Pm              <- rv$Pm[mutated == 1] + runif(sum(mutated), -1*rv$max_mutation,rv$max_mutation) 
        new_Pm[new_Pm < 0.05] <- 0.05
        
        rv$species[mutated == 1] <- new_species
        rv$Pm[mutated == 1]      <- new_Pm
      }
      
      rv$nspecies   <- sapply(rv$comm_ID, calc_n_spec)
      mn            <- mean(rv$nspecies[rv$comm_type =='sub'])
      mean_nspecies <- c(mean_nspecies, mn)
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    print(round(hist(rv$Pm, breaks = seq(-0.025, 1.025, 0.05))$counts / (rv$n*rv$n_ind), 3))
    
    x <- 1:50
    s <- summary(lm(mean_nspecies[length(mean_nspecies) - 49:0]~x))
    if (s$coefficients[2,4] > 0.05) { p <- p + 1 }
    rv$iteration_nr <- rv$iteration_nr + 50
  }
  return(rv)
} 
