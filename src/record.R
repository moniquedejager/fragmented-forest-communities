record <- function(rv){
  # record, per 50 iterations:
  # - average number of species per subcommunity 
  # - total number of (unique) species
  # - fraction of t-1 ancestors from another subcommunity
  # - fraction of t-50 ancestors from another subcommunity
  # - average dispersal capacity (Pm) 
  
  data <- data.frame(sim_nr              <- rv$sim_nr,
                     clustering          <- rv$clustering,
                     mutation_rate       <- rv$mutation_rate,
                     max_mutation        <- rv$max_mutation,
                     f_loss              <- rv$f_loss,
                     hab_cover           <- rv$hab_cover,
                     clustering_restored <- rv$clustering_restored,
                     n_iterations        <- rv$n_iterations, 
                     simulation_type     <- rv$simulation_type,
                     iteration_nr        <- rv$iteration_nr,
                     m_species           <- mean(rv$nspecies),
                     total_species       <- length(unique(rv$species)),
                     f_t1_same_subcom    <- mean(rv$comm_ID2 == rv$origin_ID_t1),
                     f_t50_same_subcom   <- mean(rv$comm_ID2 == rv$origin_ID_t50),
                     m_Pm                <- mean(rv$Pm))
  
  
  filename <- paste('results/data per iteration/sim_nr=', rv$sim_nr, 
                    'clustering=', rv$clustering,
                    'f_loss=', rv$f_loss,
                    "hab_cover=", rv$hab_cover,
                    'clustering_restored=', rv$clustering_restored, '.txt', sep='')
  if (file.exists(filename)){
    write.table(data, filename, append = TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(data, filename, col.names = TRUE, row.names = FALSE)
  }
}
