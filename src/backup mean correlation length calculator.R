# calculate mean correlation length based on the occurrences of the dispersal 
# strategies per simulation:

filename <- './results/simulation results 03-2024/dispersal_capacity/dispersal capacity_data_1.txt'
df       <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('./results/simulation results 03-2024/dispersal_capacity/dispersal capacity_data_', i, '.txt', sep='')
  df2      <- read.table(filename, header = TRUE)
  df       <- rbind(df, df2)
}

filename <- './results/simulation results 03-2024/correlation_length/correlation_length_data_1.txt'
df3       <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('./results/simulation results 03-2024/correlation_length/correlation_length_data_', i, '.txt', sep='')
  df2      <- read.table(filename, header = TRUE)
  df3       <- rbind(df3, df2)
}

df4 <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df4) <- c('mu', 'f_hab_loss', 'sim_nr', 'm_corr_length')

df$corr_length <- NA
group <- paste(df$mu, df$f_hab_loss, df$sim_nr, df$disp_cap, sep='-')
group2 <- paste(df3$mu, df3$f_hab_loss, df3$sim_nr, df3$dispersal_capacity, sep='-')
for (i in unique(group)){
  df$corr_length[group == i] <- df3$corr_length[group2 == i]
}

group <- paste(df$mu, df$f_hab_loss, df$sim_nr, sep='-')

df4 <- data.frame(mu = tapply(df$mu, group, mean),
                  f_hab_loss = tapply(df$f_hab_loss, group, mean),
                  sim_nr     = tapply(df$sim_nr, group, mean),
                  m_corr_length = tapply(df$perc_indiv*df$corr_length / 100, group, sum))

# these need to be added to the data per simulation:
for (i in 1:10){
  filename <- paste('results/simulation results 03-2024/simulation_data/simulation_data_',
                    i, '.txt', sep='')
  df      <- read.table(filename, header=T)

  group <- paste(df$mu, df$f_hab_loss, df$sim_nr, sep='-')
  group2 <- paste(df4$mu, df4$f_hab_loss, df4$sim_nr, sep='-')
  
  for (j in unique(group)){
    df$m_corr_length[group == j] <- df4$m_corr_length[group2 == j]
  }
  
  filename <- paste('results/simulation results 03-2024/new simulation_data/simulation_data_',
                    i, '.txt', sep='')
  write.table(df, filename, row.names = FALSE, col.names = TRUE)
}


