
write_data_to_files_restoration <- function(rv){
  #library(ggplot2)
  library(meteR) # rominger and merow, 2016
  library(benthos)
  source('./src/est_RAI.R')
  source('./src/calc_correlation_length.R')
  source('./src/calculate_dissimilarity.R')
  source('./src/calc_Morans_I.R')
  
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
  
  
  # calculate distance from t-50 ancestor:               
  y_origin         <- rv$y[rv$origin_ID_t50]
  x_origin         <- rv$x[rv$origin_ID_t50]
  dist_from_origin <- sqrt((rv$x2 - x_origin)^2 + (rv$y2 - y_origin)^2)
  group            <- sort(rep(rv$comm_ID[rv$comm_type == 'sub'], rv$n_ind))
  sel              <- rv$comm_type[rv$comm_ID2] == 'sub'
  dist_from_origin <- dist_from_origin[sel]
  Pm               <- tapply(dist_from_origin > 0, group, mean)
  
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

  # calculate Shannon Index per subcommunity:
  df2 <- matrix(rv$species, rv$n_ind, )
  df2 <- df2[,rv$comm_type == 'sub']
  
  H <- sapply(1:ncol(df2), function(i){
    species <- df2[,i]
    n <- tapply(species, as.factor(species), length)
    p <- n / sum(n)
    
    H <- -1*sum(p * log(p))
    return(H)
  })
  
  # Record, for the 20 subcommunities still present at 95% habitat loss, 
  # mu, % habitat loss, sim_nr, x, y, n_species, METE fit, Pm, m_dist, and
  # Shannon index H. 
  filename   <- paste('./results/landscapes/restored_landscapes_',
                      rv$sim_nr, '_', rv$clustering, '.txt', sep='')
  landscape  <- read.table(filename, header=TRUE)
  
  landscape <- landscape[(round(landscape$f_loss, 2) == 0.95)&
                  (landscape$hab_cover == 0.05)&
                  (landscape$clustering == rv$clustering_restored),]
  
  fit2 <- rep(0, rv$n)
  fit2[rv$comm_type == 'sub'] <- fit
  Pm2  <- rep(0, rv$n)
  Pm2[rv$comm_type == 'sub'] <- Pm
  H2   <- rep(0, rv$n)
  H2[rv$comm_type == 'sub'] <- H
  
  df <- data.frame(mu            = rv$clustering,
                   sim_nr        = rv$sim_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation,
                   f_loss     = rv$f_loss,
                   hab_cover  = rv$hab_cover,
                   clustering_restored = rv$clustering_restored,
                   x          = rv$x[landscape$hab==1],
                   y          = rv$y[landscape$hab==1],
                   n_species  = rv$nspecies[landscape$hab==1],
                   METE_fit   = fit2[landscape$hab==1],
                   Pm         = Pm2[landscape$hab==1],
                   Shannon    = H2[landscape$hab==1])
  
  filename <- paste('./results/subcommunity_data/restored_subcommunity_data_',
                    rv$sim_nr, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }

  # Record, per simulation:
  df <- data.frame(mu         = rv$clustering,
                   sim_nr     = rv$sim_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation,
                   f_loss     = rv$f_loss,
                   hab_cover  = rv$hab_cover,
                   clustering_restored = rv$clustering_restored,
                   m_nspecies_all_patches = mean(rv$nspecies[rv$comm_type == 'sub']),
                   m_nspecies_last_20 = mean(rv$nspecies[landscape$hab==1]),
                   m_METE_fit_all_patches = mean(fit),
                   m_METE_fit_last_20 = mean(fit2[landscape$hab==1]),
                   m_Shannon_all_patches = mean(H),
                   m_Shannon_last_20 = mean(H2[landscape$hab==1]),
                   n_species     = length(unique(rv$species[(rv$comm_type2 == 'sub')&(rv$species > 0)])),
                   n_iterations  = rv$iteration_nr,
                   mutation_rate = rv$mutation_rate,
                   max_mutation  = rv$max_mutation)
  
  filename <- paste('./results/simulation_data/restored_simulation_data_',
                    rv$sim_nr, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Record, per simulation, mu, % habitat loss, sim_nr, and % individuals per 
  # dispersal strategy
  disp_cap <- rv$Pm
  p_disp   <- tapply(disp_cap, disp_cap, length)/length(disp_cap)*100
  disp     <- tapply(disp_cap, disp_cap, mean)
  df       <- data.frame(mu         = rv$clustering,
                         sim_nr     = rv$sim_nr,
                         mutation_rate = rv$mutation_rate,
                         max_mutation  = rv$max_mutation,
                         f_loss     = rv$f_loss,
                         hab_cover  = rv$hab_cover,
                         clustering_restored = rv$clustering_restored,
                         disp_cap   = disp,
                         perc_indiv = p_disp)
  filename <- paste('./results/dispersal_capacity/restored_dispersal capacity_data_',
                    rv$sim_nr, '.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
}