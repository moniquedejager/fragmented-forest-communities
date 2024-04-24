
write_data_to_files <- function(rv){
  #library(ggplot2)
  library(meteR) # rominger and merow, 2016
  library(benthos)
  source('./src/est_RAI.R')
  #source('./src/calc_correlation_length.R')
  #source('./src/calculate_dissimilarity.R')
  #source('./src/calc_Morans_I.R')
  
  # Per community, calculate the number of species, METE fit, Pm, m_dist, 
  # and shannon index.
  # Per species, calculate correlation length.
  # Per subcommunity pair (of the 125 subcomm. left at 95% habitat loss), 
  # the Bray-Curtis dissimilarity and distance between the subcommunities.
  
  # Record, for the 125 subcommunities still present at 95% habitat loss, 
  # mu, % habitat loss, sim_nr, x, y, n_species, METE fit, Pm, m_dist, and
  # Shannon index H. 
  # Record, per simulation, mu, % habitat loss, sim_nr, Moran's I, 
  # mean patch size, n_patches, mean n_species, mean METE fit, mean Shannon
  # Index H, mean correlation length, and % individuals per dispersal strategy
  # Record, per species, mu, % habitat loss, sim_nr, species_nr, and 
  # correlation length. 
  # Record, per subcommunity pair, mu, % habitat loss, sim_nr, distance between
  # the subcommunities, and Bray-Curtis dissimilarity. 
  # Record the composition of the subcommunities per simulation for use in 
  # the restoration model. (matrix of species numbers of all subcommunities 
  # x 1000 individuals).
  
  filename  <- paste('./results/landscapes/landscapes_',
                     rv$sim_nr, '_', rv$clustering, '.txt', sep='')
  landscape <- read.table(filename)
  
  # calculate distance from t-50 ancestor:                ################## here, we do not yet take into account that the boundary conditions are continuous!!!!
  #y_origin         <- rv$y[rv$origin_ID_t50]
  #x_origin         <- rv$x[rv$origin_ID_t50]
  #dist_from_origin <- sqrt((rv$x2 - x_origin)^2 + (rv$y2 - y_origin)^2)
  #group            <- sort(rep(rv$comm_ID[rv$comm_type == 'sub'], rv$n_ind))
  #sel              <- rv$comm_type[rv$comm_ID2] == 'sub'
  #dist_from_origin <- dist_from_origin[sel]
  #m_dist           <- tapply(dist_from_origin, group, mean)
  #Pm               <- tapply(dist_from_origin > 0, group, mean)
  
  # calculate the METE fit:
  # fit <-  sapply(rv$comm_ID[rv$comm_type == 'sub'], 
  #              function(j){
  #                species   <- rv$species[j*rv$n_ind - (1:rv$n_ind) + 1]
  #                abundance <- sort(tapply(species, as.factor(species), length),
  #                                  decreasing = TRUE)
  #                
  #                if (length(abundance) >= 10){
  #                  difference <- est_RAI(abundance,'differences')
  #                  fit        <- 1 - 
  #                    sum(difference^2)/sum((abundance - mean(abundance))^2)
  #                }else { 
  #                  fit <- 0
  #                }
  #                fit
  #              })

  # calculate Shannon Index per subcommunity:
  #df2 <- matrix(rv$species, rv$n_ind, )
  #df2 <- df2[,rv$comm_type == 'sub']
  
  #H <- sapply(1:ncol(df2), function(i){
  #  species <- df2[,i]
  #  n <- tapply(species, as.factor(species), length)
  #  p <- n / sum(n)
  #  
  #  H <- -1*sum(p * log(p))
  #  return(H)
  #})
  
  # calculate spatial correlation length per species:
  #C <- calc_correlation_length(rv)                                              
  
  # Calculate dissimilarity between subcommunity pairs:
  #dissimilarity <- calc_dissimilarity(rv)
  
  # calculate Moran's I, mean patch size, and number of patches:
  #morans_I <- calc_Morans_I(rv)
  
  # Record, for the 125 subcommunities still present at 95% habitat loss, 
  # mu, % habitat loss, sim_nr, x, y, n_species, METE fit, Pm, m_dist, and
  # Shannon index H. 
  iFrag <- (1:20)[seq(0, 95, 5) == rv$frag]
  df <- data.frame(mu         = rv$clustering,
                   f_hab_loss = rv$frag,
                   sim_nr     = rv$sim_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation,
                   x          = rv$x[landscape[,20]==1],
                   y          = rv$y[landscape[,20]==1],
                   n_species  = rv$nspecies[landscape[,20]==1])
                   
  filename <- paste('./results/subcommunity_data/subcommunity_data_',
                    rv$sim_nr, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }

  # Record, per simulation, mu, % habitat loss, sim_nr, Moran's I,            
  # mean patch size, n_patches, mean n_species, mean METE fit, mean Shannon
  # Index H, and mean correlation length
  df <- data.frame(mu         = rv$clustering,
                   f_hab_loss = rv$frag,
                   sim_nr     = rv$sim_nr,
                   #Morans_I   = morans_I[1],
                   #m_patch_size = morans_I[2],
                   #n_patches    = morans_I[3],
                   m_nspecies   = mean(rv$nspecies[rv$comm_type == 'sub']),
                   #m_METE_fit   = mean(fit),
                   #m_Shannon    = mean(H),
                   #m_corr_length = mean(C),
                   n_species     = length(unique(rv$species[rv$comm_type2 == 'sub'])),
                   n_iterations  = rv$iteration_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation)

  filename <- paste('./results/simulation_data/simulation_data_',
                    rv$sim_nr, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Record, per simulation, mu, % habitat loss, sim_nr, and % individuals per 
  # dispersal strategy
  # disp_cap <- rv$Pm
  # p_disp   <- tapply(disp_cap, disp_cap, length)/length(disp_cap)*100
  # disp     <- tapply(disp_cap, disp_cap, mean)
  # df       <- data.frame(mu         = rv$clustering,
  #                        f_hab_loss = rv$frag,
  #                        sim_nr     = rv$sim_nr,
  #                        mutation_rate = rv$mutation_rate,
  #                        max_mutation  = rv$max_mutation,
  #                        disp_cap   = disp,
  #                        perc_indiv = p_disp)
  # filename <- paste('./results/dispersal_capacity/dispersal capacity_data_',
  #                   rv$sim_nr, '.txt', sep='')
  # if (file.exists(filename)){
  #   write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  # } else {
  #   write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  # }
  
  # Record, per dispersal capacity, mu, % habitat loss, sim_nr, species_nr, and 
  # correlation length. 
  
  # df <- data.frame(mu         = rv$clustering,
  #                  f_hab_loss = rv$frag,
  #                  sim_nr     = rv$sim_nr,
  #                  mutation_rate = rv$mutation_rate,
  #                  max_mutation  = rv$max_mutation,
  #                  dispersal_capacity = unique(rv$Pm),
  #                  corr_length= C)
  # filename <- paste('./results/correlation_length/correlation_length_data_',
  #                   rv$sim_nr, '.txt', sep='')
  # if (file.exists(filename)){
  #   write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  # } else {
  #   write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  # }
  # 
  # Record, per subcommunity pair, mu, % habitat loss, sim_nr, distance between
  # the subcommunities, and Bray-Curtis dissimilarity. 
  # df <- data.frame(mu            = rv$clustering,
  #                  f_hab_loss    = rv$frag,
  #                  sim_nr        = rv$sim_nr,
  #                  mutation_rate = rv$mutation_rate,
  #                  max_mutation  = rv$max_mutation,
  #                  distance      = dissimilarity[,1],
  #                  dissimilarity = dissimilarity[,2])
  # filename <- paste('./results/dissimilarity/dissimilarity_data_',
  #                   rv$sim_nr, '.txt', sep='')
  # if (file.exists(filename)){
  #   write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  # } else {
  #   write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  # }
  
  # Record the composition of the subcommunities per simulation for use in 
  # the restoration model. (matrix of species numbers of all subcommunities 
  # x 1000 individuals).
  # df <- matrix(paste(rv$species, rv$Pm, sep='-'), rv$n_ind, )
  # df <- df[,rv$comm_type == 'sub']
  # filename <- paste('./results/community_composition/',
  #                   rv$clustering,'_',rv$sim_nr, '_',rv$frag, 
  #                   '_', rv$mutation_rate,'_', rv$max_mutation,'.txt', sep='')
  # write.table(df, filename, append=FALSE, col.names = FALSE, row.names = FALSE)
}