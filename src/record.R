record <- function(rv){
  # record, per 50 iterations:
  # - average number of species per subcommunity 
  # - total number of (unique) species
  # - fraction of t-1 ancestors from another subcommunity
  # - fraction of t-50 ancestors from another subcommunity
  # - average dispersal capacity (Pm) 
  # - average METE fit to SAD
  
  library(meteR) # rominger and merow, 2016
  source('./src/est_RAI.R')
  
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
  
  
  
  data <- data.frame(sim_nr              = rv$sim_nr,
                     clustering          = rv$clustering,
                     mutation_rate       = rv$mutation_rate,
                     max_mutation        = rv$max_mutation,
                     f_loss              = rv$f_loss,
                     simulation_type     = rv$simulation_type,
                     iteration_nr        = rv$iteration_nr,
                     dispersal           = rv$dispersal,
                     hab_cover           = rv$hab_cover,
                     clustering_restored = rv$clustering_restored,
                     m_species           = mean(rv$nspecies[rv$comm_type == 'sub']),
                     m_METE_fit          = mean(fit),
                     total_species       = length(unique(rv$species[rv$species > 0])),
                     f_t50_same_subcom   = mean(rv$comm_ID2[rv$comm_type2 == 'sub'] 
                                                == rv$origin_ID_t50[rv$comm_type2 == 'sub']),
                     m_Pm                = mean(rv$Pm[(rv$comm_type2 == 'sub')&(rv$Pm > 0)]))
  
  filename <- paste('results/data per iteration/sim_nr=', rv$sim_nr, 
                    'clustering=', rv$clustering,
                    'f_loss=', rv$f_loss, 
                    'mutation_rate=',rv$mutation_rate,
                    'restoration', rv$hab_cover,
                    'b.txt', sep='')
  if (file.exists(filename)){
    write.table(data, filename, append = TRUE, col.names = FALSE, row.names = FALSE)
  } else {
    write.table(data, filename, col.names = TRUE, row.names = FALSE)
  }
}
