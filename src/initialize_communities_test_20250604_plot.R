library(ggplot2)
source('./src/summarySE.R')

df <- read.table('results/dispersal_capacity/test_initial_dispersal_distribution1_1.txt', header=TRUE)
df <- df[df$mu == 3,]

for (i in 1:5){
  for (j in 1:3){
    filename <- paste('./results/dispersal_capacity/test_initial_dispersal_distribution',
                     j,'_', i,'.txt', sep='')
    df <- rbind(df, read.table(filename, header=TRUE)) 
  }
}


ggplot(df, aes(x=disp_cap, y=perc_indiv, color=factor(disp_distribution))) + 
  geom_point() + 
  geom_line(aes(linetype=factor(sim_nr)))

sdf <- summarySE(df, measurevar="perc_indiv", 
                 groupvars=c("disp_distribution","disp_cap"))

sdf$type <- 'initial peak at λ = 0.1'
sdf$type[sdf$disp_distribution == 1] <- 'uniform initial distribution'
sdf$type[sdf$disp_distribution == 3] <- 'initial peak at λ = 0.9'

cbp1 <- c("#009E73", "#E69F00", "#999999")
pd <- position.dodge(0.1)
p1 <- ggplot(sdf, aes(x=disp_cap, y=perc_indiv, color=factor(type))) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin=perc_indiv-se, ymax=perc_indiv+se), 
                width=0, position=pd) + 
  scale_color_manual(values=cbp1, name='') + 
  ylab('% Individuals') + 
  xlab('Dispersal capacity λ') + 
  theme_bw() + 
  theme(legend.position = 'top')
p1

# 6x4
tiff(filename = 'figures/initialize communities test 20250604.tif', 
     width = 6, height = 4, units = 'in', res = 600)
p1
dev.off()


# checking with the old simulations:
df2 <- read.table('results/dispersal_capacity/fragmented_dispersal capacity_data_0.txt', header=TRUE)
summary(df2)

sel <- df2$dispersal == 'different'
ggplot(df2[sel,], aes(disp_cap, y=perc_indiv, color=factor(mu))) + 
  geom_point() + 
  geom_line(aes(linetype=factor(sim_nr)))

unique(df2$disp_cap)
unique(df2$dispersal)
