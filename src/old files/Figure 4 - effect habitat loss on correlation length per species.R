library(ggplot2)

filename <- './results/simulation results 04-2024/correlation_length/correlation_length_data_1.txt'
df       <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('./results/simulation results 04-2024/correlation_length/correlation_length_data_', i, '.txt', sep='')
  df2      <- read.table(filename, header = TRUE)
  df       <- rbind(df, df2)
}


ggplot(df, aes(x=f_hab_loss, y=corr_length)) + 
  geom_point(alpha = 0.3, color='salmon') + 
  facet_grid(rows=vars(dispersal_capacity), cols=vars(mu)) + 
  theme_classic()

df2 <- df[(df$mu != 2)&(df$dispersal_capacity %in% c(0.2, 0.3, 0.4)),]
  
df2$clustering <- 'Random'
df2$clustering[df2$mu == 3] <- 'Fractal'
df2$clustering[df2$mu == 5] <- 'Clustered'
#df2$clustering <- factor(df2$clustering, levels=unique(df2$clustering))
df2$disp <- "λ = 0.2"
df2$disp[df2$dispersal_capacity == 0.3] <- "λ = 0.3"
df2$disp[df2$dispersal_capacity == 0.4] <- "λ = 0.4"
ggplot(df2, aes(x=f_hab_loss, y=corr_length)) + 
  geom_point(alpha = 0.3, color='darkcyan') + 
  facet_grid(rows=vars(disp), cols=vars(clustering)) + 
  scale_y_continuous(trans='log10') + 
  theme_bw() + 
  xlab('% Habitat loss') + 
  ylab('Correlation length (C)') + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

group <- paste(df2$f_hab_loss, df2$mu, df2$dispersal_capacity,sep='-')
df3   <- data.frame(loss = tapply(df2$f_hab_loss, group, mean),
                  mu   = tapply(df2$mu, group, mean),
                  C    = tapply(df2$corr_length, group, median), 
                  disp = tapply(df2$dispersal_capacity, group, mean))

df3$clustering <- 'Random'
df3$clustering[df3$mu == 3] <- 'Fractal'
df3$clustering[df3$mu == 5] <- 'Clustered'
df3$clustering <- factor(df3$clustering, levels=unique(df3$clustering))
df3$disp2 <- "λ = 0.2"
df3$disp2[df3$disp == 0.3] <- "λ = 0.3"
df3$disp2[df3$disp == 0.4] <- "λ = 0.4"
ggplot(df3, aes(x=loss, y=C, color=clustering)) + 
  geom_point(alpha = 0.3) + 
  geom_line() + 
  facet_grid(cols=vars(disp2)) + 
  #scale_y_continuous(trans='log10') + 
  theme_bw() + 
  xlab('% Habitat loss') + 
  ylab('Correlation length (C)')
