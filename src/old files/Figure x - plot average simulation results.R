# plot the data per simulation:
library(ggplot2)
library(ggpubr)
source('./src 03-2024/summarySE.R')

filename <- 'results/simulation results 04-2024/simulation_data/simulation_data_1.txt'
df       <- read.table(filename, header=T)

for (i in 2:10){
  filename <- paste('results/simulation results 04-2024/simulation_data/simulation_data_',
                    i, '.txt', sep='')
  df2      <- read.table(filename, header=T)
  df       <- rbind(df, df2)
}

#df$mu2             <- "μ = 1"
#df$mu2[df$mu == 2] <- "μ = 2"
#df$mu2[df$mu == 3] <- "μ = 3"
#df$mu2[df$mu == 4] <- "μ = 4"

df <- df[df$mu != 2,]
df$mu2 <- 'Random'
df$mu2[df$mu == 3] <- 'Fractal'
df$mu2[df$mu == 5] <- 'Clustered'

sdf <- summarySE(df, measurevar="nspecies", 
                 groupvars=c("mu2","f_hab_loss"))

pd <- position_dodge(0.01) # move them .05 to the left and right

p1 <- ggplot(sdf, aes(x=f_hab_loss, y=nspecies, color=mu2)) + 
  geom_errorbar(aes(ymin=nspecies-se, ymax=nspecies+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab('Total # species') +
  scale_y_continuous(trans='log10') + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf2 <- summarySE(df, measurevar="morans_I", 
                  groupvars=c("mu2","f_hab_loss"))

p2 <- ggplot(sdf2, aes(x=f_hab_loss, y=morans_I, color=mu2)) + 
  geom_errorbar(aes(ymin=morans_I-se, ymax=morans_I+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Moran's I") + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf3 <- summarySE(df, measurevar="m_patch_size", 
                  groupvars=c("mu2","f_hab_loss"))

p3 <- ggplot(sdf3, aes(x=f_hab_loss, y=m_patch_size, color=mu2)) + 
  geom_errorbar(aes(ymin=m_patch_size-se, ymax=m_patch_size+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Mean patch size") +
  scale_y_continuous(trans='log10') + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf4 <- summarySE(df, measurevar="n_patches", 
                  groupvars=c("mu2","f_hab_loss"))

p4 <- ggplot(sdf4, aes(x=f_hab_loss, y=n_patches, color=mu2)) + 
  geom_errorbar(aes(ymin=n_patches-se, ymax=n_patches+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Number of patches") + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf5 <- summarySE(df, measurevar="m_nspecies", 
                  groupvars=c("mu2","f_hab_loss"))

p5 <- ggplot(sdf5, aes(x=f_hab_loss, y=m_nspecies, color=mu2)) + 
  geom_errorbar(aes(ymin=m_nspecies-se, ymax=m_nspecies+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Mean # species") + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf6 <- summarySE(df, measurevar="m_METE_fit", 
                  groupvars=c("mu2","f_hab_loss"))

p6 <- ggplot(sdf6, aes(x=f_hab_loss, y=m_METE_fit, color=mu2)) + 
  geom_errorbar(aes(ymin=m_METE_fit-se, ymax=m_METE_fit+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Mean METE fit to SAD") + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf7 <- summarySE(df, measurevar="m_Shannon", 
                  groupvars=c("mu2","f_hab_loss"))

p7 <- ggplot(sdf7, aes(x=f_hab_loss, y=m_Shannon, color=mu2)) + 
  geom_errorbar(aes(ymin=m_Shannon-se, ymax=m_Shannon+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Mean Shannon Index") + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

sdf8 <- summarySE(df, measurevar="m_corr_length", 
                  groupvars=c("mu2","f_hab_loss"))

p8 <- ggplot(sdf8, aes(x=f_hab_loss, y=m_corr_length, color=mu2)) + 
  geom_errorbar(aes(ymin=m_corr_length-se, ymax=m_corr_length+se), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  xlab('') + 
  ylab("Mean correlation length") + 
  theme_bw() + 
  theme(legend.position = 'top') + 
  theme(legend.title=element_blank())

figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4, nrow=2, 
          common.legend = TRUE, align = "hv", 
          labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'))

annotate_figure(figure,
                bottom = text_grob("% Habitat loss", vjust=-1))


figure2 <- ggarrange(p2, p3, p4, ncol=3, nrow=1, common.legend = TRUE, 
                     labels = c('A', 'B', 'C'))
annotate_figure(figure2,
                bottom = text_grob("% Habitat loss", vjust=-1))

figure3 <- ggarrange(p5, p6, p7, p8, ncol=2, nrow=2, common.legend = TRUE, 
                     labels = c('A', 'B', 'C', 'D'), align='hv')
annotate_figure(figure3,
                bottom = text_grob("% Habitat loss", vjust=-1))

