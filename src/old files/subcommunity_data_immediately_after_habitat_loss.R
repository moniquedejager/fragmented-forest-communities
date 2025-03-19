# Find out what the immediate effects of habitat fragmentation and loss
# are on biodiversity at the subcommunity level:

f_loss <- 0
filename <- paste('./results/subcommunity_data/fragmented_subcommunity_data_',
                  f_loss, '.txt', sep='')
m1 <- read.table(filename, header=TRUE)
m1$loc <- paste(m1$x, m1$y, sep='-')

for (f_loss in 1:19/20){
  filename <- paste('./results/subcommunity_data/fragmented_subcommunity_data_',
                    f_loss, '.txt', sep='')
  m2 <- read.table(filename, header=TRUE)
  m2$loc <- paste(m2$x, m2$y, sep='-')
  
  df <- m2[m2$mu == 10,]
  
  for (sim_nr in 1:10){
    for (mu in c(1, 3, 5)){
      for (sim_type in c('different', 'similar')){
        sel <- (m2$mu == mu)&(m2$sim_nr == sim_nr)&(m2$dispersal == sim_type)
        m3 <- m2[sel,]
        
        sel <- (m1$mu == mu)&(m1$sim_nr == sim_nr)&(m1$dispersal == sim_type)
        m0 <- m1[sel,]
        
        m3$n_species <- m0$n_species[m0$loc %in% m3$loc]
        m3$Pm <- m0$Pm[m0$loc %in% m3$loc]
        
        df <- rbind(df, m3)
      }
    }
  }
  
  filename <- paste('./results/subcommunity_data/fragmented_subcommunity_data_',
                    f_loss, 'immediate_after_hab_loss.txt', sep='')
  write.table(df, filename, append=FALSE, col.names = TRUE, row.names=FALSE)
}




