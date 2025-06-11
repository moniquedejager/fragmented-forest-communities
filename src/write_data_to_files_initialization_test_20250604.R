
write_data_to_files_initialization_test_20250604 <- function(rv){
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

  # Record, per simulation, mu, % habitat loss, sim_nr, dispersal distribution, 
  # and % individuals per dispersal strategy
  disp_cap <- rv$Pm[rv$species > 0]
  p_disp   <- tapply(disp_cap, disp_cap, length)/length(disp_cap)*100
  disp     <- tapply(disp_cap, disp_cap, mean)
  df       <- data.frame(mu            = rv$clustering,
                         sim_nr        = rv$sim_nr,
                         disp_distribution = rv$disp_distribution,
                         mutation_rate = rv$mutation_rate,
                         max_mutation  = rv$max_mutation,
                         f_loss        = rv$f_loss,
                         n_iterations  = rv$iteration_nr,
                         dispersal     = rv$dispersal,
                         disp_cap      = disp,
                         perc_indiv    = p_disp)
  filename <- paste('./results/dispersal_capacity/test_initial_dispersal_distribution',
                    rv$disp_distribution,'_', rv$sim_nr,'.txt', sep='')
  if (file.exists(filename)){
    write.table(df, filename, append=TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(df, filename, append=FALSE, col.names = TRUE, row.names = FALSE)
  }
  
}