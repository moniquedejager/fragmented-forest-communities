write_species_to_files <- function(rv){
 m <- matrix(rv$species, rv$n_ind, )
 filename <- paste('./results/METE_in_fragmented_landscapes/community_composition/',
                  rv$clustering,'_',rv$sim_nr, '_',rv$frag, 
                  '_', rv$mutation_rate,'_','.txt', sep='')
 write.table(m, filename, append=FALSE, col.names = FALSE, row.names = FALSE)
}