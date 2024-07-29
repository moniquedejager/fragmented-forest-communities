# Species-area curves:
library(ggplot2)
library(ggpubr)
source('./src/summarySE.R')

filename <- 'results/simulation_data/fragmented_simulation_data_0.txt'
df <- read.table(filename, header=T)

for (i in seq(0, 0.95, 0.05)){
  filename <- paste('results/simulation_data/fragmented_simulation_data_',
                    i,'.txt', sep='')
  df2 <- read.table(filename, header=T)
  df  <- rbind(df, df2)
  
  #df <- df2[df2$dispersal == 'different',]
  #write.table(df, filename, append = FALSE, row.names = FALSE, col.names = TRUE)
}

df <- data.frame(clustering = rep(df$mu, 2),
                 area_size = rep(400*(1-df$f_loss), 2),
                 n_species = c(df$n_species, df$static_n_species),
                 dispersal_type = rep(df$dispersal),
                 static_dynamic = c(rep('Dynamic', length(df$mu)), 
                                    rep('Static', length(df$mu))))

for (sim_nr in 1:10){
  for (mu in c(1, 3, 5)){
    for (disp in c('different', 'similar')){
      filename <- paste('results/community_composition/', 
                        disp, mu, '_', sim_nr, '_0_3e-04_0.txt', sep='')
      if (file.exists(filename)){
        m <- read.table(filename)
        landscape <- as.matrix(m)
        for (i in 1:20){
          # backwards SAR:
          area_size <- i * 20
          n_species <- length(unique(as.vector(landscape[1:area_size,])))
          
          df2 <- data.frame(clustering = 0,
                            area_size = area_size,
                            n_species = n_species,
                            dispersal_type = disp,
                            static_dynamic = c('Dynamic'))
          df <- rbind(df, df2)
        }
      }
    }
  }
}

filenames <- c('results/community_composition/initial_community_sim.txt',
               'results/community_composition/initial_community_dif.txt')
disp <- c('similar', 'different')
for (j in 1:2){
  m <- read.table(filenames[j])
  landscape <- as.matrix(m)
  for (i in 1:20){
    # backwards SAR:
    area_size <- i * 20
    n_species <- length(unique(as.vector(landscape[1:area_size,])))
    
    df2 <- data.frame(clustering = 0,
                      area_size = area_size,
                      n_species = n_species,
                      dispersal_type = disp[j],
                      static_dynamic = c('Static'))
    df <- rbind(df, df2)
  }
}

sdf <- summarySE(df, measurevar="n_species", 
                 groupvars=c("clustering","area_size", "dispersal_type", "static_dynamic"))

sdf$mu2 <- 'Random'
sdf$mu2[sdf$clustering == 3] <- 'Fractal'
sdf$mu2[sdf$clustering == 5] <- 'Clustered'
sdf$mu2[sdf$clustering == 0] <- 'Backwards SAR'

sdf$static_dynamic2 <- 'Static'
sdf$static_dynamic2[sdf$static_dynamic == 'Dynamic'] <- 'Dynamic'
sdf$static_dynamic2 <- factor(sdf$static_dynamic2, levels=c('Static', 'Dynamic'))

sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal_type == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

pd <- position_dodge(0.01) 
ggplot(sdf, aes(x=area_size, y=n_species, color=mu2)) + 
  geom_errorbar(aes(ymin=n_species-sd, ymax=n_species+sd), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  scale_color_manual(values = c('grey40', "#F8766D", "#00BA38", "#619CFF")) + 
  facet_grid(cols=vars(disp_type), rows=vars(static_dynamic2), scales = 'free_y') + 
  xlab('Area size') + 
  ylab('Total # species') +
  theme_bw() + 
  theme(legend.position = 'top',
      strip.placement = "outside", 
      strip.background = element_blank(),
      legend.title=element_blank())

pd <- position_dodge(0.1) 
sdf$hab_loss <- (1 - sdf$area_size/400)*100
p1 <- ggplot(sdf, aes(x=hab_loss, y=n_species, color=mu2)) + 
  geom_errorbar(aes(ymin=n_species-sd, ymax=n_species+sd), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  #scale_x_continuous(trans='log10') + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values = c('grey40', "#F8766D", "#00BA38", "#619CFF")) + 
  facet_grid(cols=vars(disp_type), rows=vars(static_dynamic2)) + 
  xlab('% Habitat loss') + 
  ylab('Total # species') +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())

# tiff file 600 dpi:
tiff(filename = 'figures/Figure 2.tif', 
     width = 5, height = 5, units = 'in', res = 600)
p1
dev.off()



range(sdf$area_size)

ggplot(sdf, aes(x=area_size, y=n_species, color=dispersal_type)) + 
  geom_errorbar(aes(ymin=n_species-ci, ymax=n_species+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  facet_grid(cols=vars(mu2), rows=vars(static_dynamic2), scales = 'free_y') + 
  xlab('Area size') + 
  ylab('Total # species') +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())


m <- read.table('results/community_composition/initial_community_dif.txt')
m <- as.vector(t(as.matrix(m)))
m2<- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2 - 1])
Pm <- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2])

hist(Pm)
mean(Pm)


m <- read.table('results/community_composition/different1_1_0_3e-04_0.txt')
m <- t(as.matrix(m))
m2 <- m[7:1006,]
m <- as.vector(m2)

m2<- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2 - 1])
Pm <- as.numeric(unlist(strsplit(m, '-'))[(1:(length(m)))*2])

hist(Pm)
mean(Pm)
