# Simulation output for Fleur!

filenames <- list.files('Fragmented-forest-communities/x64/Release/record one subcommunity/')

filename2 <- paste('Fragmented-forest-communities/x64/Release/record one subcommunity/', filenames[1], sep='')
df        <- read.table(filename2)
df$sim_nr <- 1

for (i in 2:5){
  filename2 <- paste('Fragmented-forest-communities/x64/Release/record one subcommunity/', filenames[i], sep='')
  df2        <- read.table(filename2)
  df2$sim_nr <- i
  df <- rbind(df, df2)
}

names(df) <- c('t', 'species', 'n', 'sim_nr')

df$t <- df$t - 9899

# let's try out some species:
sel <- (df$species < 10)&(df$sim_nr == 1)
ggplot(df[sel,], aes(x=t, y=n, color=factor(species))) + geom_line()


sel <- (df$sim_nr == 1)&(df$species > 3000)
ggplot(df[sel,], aes(x=t, y=n, color=factor(species))) + geom_line() + 
  theme(legend.position = 'none')


# write the data to a file for Fleur:
write.table(df, 
            'results/simulation results for Fleur.txt',
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE)

df[df$n == max(df$n),]
