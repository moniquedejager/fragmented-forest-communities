# Species-area curves:
library(ggplot2)
library(ggpubr)
source('./src 03-2024/summarySE.R')

df <- data.frame(clustering = vector(length=0),
                 sim_nr = vector(length=0),
                 frag = vector(length=0),
                 mutation_rate = vector(length=0),
                 max_mutation = vector(length=0),
                 area_size = vector(length=0),
                 n_species = vector(length=0),
                 dispersal_type = vector(length=0),
                 static_dynamic = vector(length=0))

mutation_rate <- 0.0001
for (clustering in c(1, 3, 5)){
  for (sim_nr in 1:10){
    for (frag in seq(0, 95, 5)){
      for (i_disp_type in 1:2){
        if (i_disp_type == 1){
          disp_type <-  'Same dispersal strategy'
          max_mutation <- 0
          filename <- paste('./results/simulation results 04-2024/community_composition/',
                            clustering,'_',sim_nr, '_',frag, 
                            '_', mutation_rate,'_', max_mutation,'similar_dispersal.txt', sep='')
        } else {
          disp_type <- 'Different dispersal strategies' 
          max_mutation <- 0.05
          filename <- paste('./results/simulation results 04-2024/community_composition/',
                            clustering,'_',sim_nr, '_',frag, 
                            '_', mutation_rate,'_', max_mutation,'.txt', sep='')
        }
      
        m <- as.matrix(read.table(filename))
        area_size <- ncol(m)
        n_species <- length(unique(as.vector(m)))
        
        df2 <- data.frame(clustering = clustering,
                         sim_nr = sim_nr,
                         frag = frag,
                         mutation_rate = mutation_rate,
                         max_mutation = max_mutation,
                         area_size = area_size,
                         n_species = n_species,
                         dispersal_type = disp_type,
                         static_dynamic = 'dynamic')
        df <- rbind(df, df2)
        
        if (frag == 0){
          # static species loss:
          filename2 <- paste('./results/simulation results 04-2024/landscapes/landscapes_',
                             sim_nr, '_', clustering, '.txt', sep='')
          landscape <- as.matrix(read.table(filename2))
          
          for (i in 1:20){
            # static species loss:
            area_size <- ncol(m[,landscape[,i] == 1])
            n_species <- length(unique(as.vector(m[,landscape[,i] == 1])))
            
            df2 <- data.frame(clustering = clustering,
                              sim_nr = sim_nr,
                              frag = frag,
                              mutation_rate = mutation_rate,
                              max_mutation = max_mutation,
                              area_size = area_size,
                              n_species = n_species,
                              dispersal_type = disp_type,
                              static_dynamic = 'static')
            df <- rbind(df, df2)
            
            # backwards SAR:
            area_size <- i * 20
            n_species <- length(unique(as.vector(m[,1:area_size])))
            
            df2 <- data.frame(clustering = 0,
                              sim_nr = sim_nr,
                              frag = frag,
                              mutation_rate = mutation_rate,
                              max_mutation = max_mutation,
                              area_size = area_size,
                              n_species = n_species,
                              dispersal_type = disp_type,
                              static_dynamic = c('static', 'dynamic'))
            df <- rbind(df, df2)
            
          }
        }
      }
    }
  }
}

sdf <- summarySE(df, measurevar="n_species", 
                 groupvars=c("clustering","area_size", "dispersal_type", "static_dynamic"))

sdf$mu2 <- 'Random'
sdf$mu2[sdf$clustering == 3] <- 'Fractal'
sdf$mu2[sdf$clustering == 5] <- 'Clustered'
sdf$mu2[sdf$clustering == 0] <- 'Backwards SAR'

sdf$static_dynamic2 <- 'Static'
sdf$static_dynamic2[sdf$static_dynamic == 'dynamic'] <- 'Dynamic'
sdf$static_dynamic2 <- factor(sdf$static_dynamic2, levels=c('Static', 'Dynamic'))

sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal_type == 'Different dispersal strategies'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

pd <- position_dodge(0.01) 
ggplot(sdf, aes(x=area_size, y=n_species, color=mu2)) + 
  geom_errorbar(aes(ymin=n_species-ci, ymax=n_species+ci), width=0, position=pd) + 
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
