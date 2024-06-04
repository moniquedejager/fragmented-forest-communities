
write_data_to_files_fragmentation <- function(rv){
  #library(ggplot2)
  library(meteR) # rominger and merow, 2016
  library(benthos)
  source('./src/est_RAI.R')
  source('./src/calculate_dissimilarity.R')

  # Per metacommunity, record: 
  # - the total number of species
  # - the frequency distribution of dispersal strategies
  # 
  # Per subcommunity, record:
  # - the number of species
  # - % ancestors from elsewhere
  # - METE's fit to SAD
  # 
  # For a subset of subcommunities, record: 
  # - the Bray-Curtis dissimilarity. 

  # calculate distance from t-50 ancestor:   
  Pm <- tapply(rv$comm_ID2 == rv$origin_ID_t50, rv$comm_ID2, mean)
  Pm <- Pm[rv$comm_type == 'sub']
  
  # Calculate dissimilarity between subcommunity pairs:
  dissimilarity <- calc_dissimilarity(rv)
  
  # calculate the METE fit:
  fit <-  sapply(rv$comm_ID[rv$comm_type == 'sub'], 
               function(j){
                 species   <- rv$species[j*rv$n_ind - (1:rv$n_ind) + 1]
                 species   <- species[species > 0]
                 abundance <- sort(tapply(species, as.factor(species), length),
                                   decreasing = TRUE)
                 
                 if (length(abundance) >= 10){
                   difference <- est_RAI(abundance,'differences')
                   fit        <- 1 - 
                     sum(difference^2)/sum((abundance - mean(abundance))^2)
                 }else { 
                   fit <- 0
                 }
                 fit
               })

  # Record per subcommunity, mu, % habitat loss, sim_nr, x, y, n_species, METE fit, and Pm
  df <- data.frame(mu            = rv$clustering,
                   sim_nr        = rv$sim_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation,
                   f_loss     = rv$f_loss,
                   n_iteration = rv$iteration_nr,
                   x          = rv$x[rv$comm_type == 'sub'],
                   y          = rv$y[rv$comm_type == 'sub'],
                   n_species  = rv$nspecies[rv$comm_type == 'sub'],
                   METE_fit   = fit,
                   Pm         = Pm)
  
  filename <- paste('./results/subcommunity_data/fragmented_subcommunity_data_',
                    rv$f_loss, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }

  # Record, per subcommunity pair, mu, % habitat loss, sim_nr, distance between
  # the subcommunities, and Bray-Curtis dissimilarity. 
  df <- data.frame(mu            = rv$clustering,
                   f_hab_loss    = rv$f_loss,
                   sim_nr        = rv$sim_nr,
                   n_iterations  = rv$iteration_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation,
                   distance      = dissimilarity[,1],
                   dissimilarity = dissimilarity[,2])
  filename <- paste('./results/dissimilarity/dissimilarity_data_',
                    rv$sim_nr, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Record, per simulation:
  df <- data.frame(mu            = rv$clustering,
                   sim_nr        = rv$sim_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation,
                   f_loss        = rv$f_loss,
                   m_nspecies    = mean(rv$nspecies[rv$comm_type == 'sub']),
                   m_METE_fit    = mean(fit),
                   n_species     = length(unique(rv$species[(rv$comm_type2 == 'sub')&(rv$species > 0)])),
                   n_iterations  = rv$iteration_nr)
  
  filename <- paste('./results/simulation_data/fragmented_simulation_data_',
                    rv$f_loss, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Record, per simulation, mu, % habitat loss, sim_nr, and % individuals per 
  # dispersal strategy
  disp_cap <- rv$Pm[rv$species > 0]
  p_disp   <- tapply(disp_cap, disp_cap, length)/length(disp_cap)*100
  disp     <- tapply(disp_cap, disp_cap, mean)
  df       <- data.frame(mu            = rv$clustering,
                         sim_nr        = rv$sim_nr,
                         mutation_rate = rv$mutation_rate,
                         max_mutation  = rv$max_mutation,
                         f_loss        = rv$f_loss,
                         n_iterations  = rv$iteration_nr,
                         disp_cap      = disp,
                         perc_indiv    = p_disp)
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_',
                    rv$f_loss, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Record the composition of the subcommunities per simulation for use in 
  # the restoration model. (matrix of species numbers of all subcommunities 
  # x 1000 individuals).
  df <- t(matrix(paste(rv$species, rv$Pm, sep='-'), rv$n_ind, rv$n))
  df <- df[rv$comm_type == 'sub',]
  df2 <- cbind(rv$clustering, 
               rv$sim_nr, 
               rv$x[rv$comm_type == 'sub'], 
               rv$y[rv$comm_type == 'sub'], 
               rv$f_loss,
               rv$iteration_nr)
  df <- cbind(df2,
              df)
  
  filename <- paste('./results/community_composition/',
                    rv$clustering,'_',rv$sim_nr, '_',rv$f_loss, 
                    '_', rv$mutation_rate,'_', rv$max_mutation,'.txt', sep='')
  write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
}