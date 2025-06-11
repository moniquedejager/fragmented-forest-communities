library(ggplot2)
library(ggpubr)

filename <- 'results/subcommunity_data/fragmented_subcommunity_data_0.25.txt'
df <- read.table(filename, header = TRUE)

for (i in c(0.5, 0.75, 0.95)){
  filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    df <- rbind(df, m)
    
    #df <- m[m$dispersal == 'different',]
    #write.table(df, filename, append = FALSE, row.names = FALSE, col.names = TRUE)
  }
}
df <- df[df$dispersal == 'different',]
df$maxLambda <- 1
df$dispKernel <- 'Original'

for (i in c(0.25, 0.5, 0.75, 0.95)){
  for (j in c(0.5, 1.5)){
    filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '_maxLambda=', j,'.txt', sep='')
    if (file.exists(filename)){
      m  <- read.table(filename, header = TRUE)
      m$maxLambda <- j
      m$dispKernel <- 'Original'
      df <- rbind(df, m)
    }
  }
}

for (i in c(0.25, 0.5, 0.75, 0.95)){
  for (j in c(0.5, 1, 1.5)){
    filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '_maxLambda=', j,'_dispKernel=exponential.txt', sep='')
    if (file.exists(filename)){
      m  <- read.table(filename, header = TRUE)
      m$maxLambda <- j
      m$dispKernel <- 'Exponential'
      df <- rbind(df, m)
    }
  }
}

for (i in c(0.25, 0.5, 0.75, 0.95)){
  for (j in c('pareto', 'gaussian')){
    filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '_maxLambda=', 1,'_dispKernel=', j, '.txt', sep='')
    if (file.exists(filename)){
      m  <- read.table(filename, header = TRUE)
      m$maxLambda <- 1
      m$dispKernel <- j
      df <- rbind(df, m)
    }
  }
}

df$type <- paste(df$dispKernel, '(0 < λ ≤ ',df$maxLambda, ')',sep='')
sel <- (df$mu != 3)
ggplot(df[sel,], aes(x=factor(f_loss), y=n_species, color=type)) + 
  geom_boxplot() + 
  facet_wrap(vars(mu))

# to show that the size of the dispersal neighbourhood does not affect the 
# results when an exponential distribution is used: 
cbp1 <- c("#009E73", "#E69F00", "#999999")
df$dispNhood <- '21 x 21 cells'
df$dispNhood[df$dispKernel == 'Original'] <- '11 x 11 cells'
df$mu2 <- 'Random habitat destruction'
df$mu2[df$mu == 5] <- 'Clustered'
df$maxLambda2 <- paste('0 < λ ≤ ', df$maxLambda, sep='')
sel <- (df$mu != 3)&(df$dispKernel %in% c('Original', 'Exponential'))
p1 <- ggplot(df[sel,], aes(x=factor(f_loss), y=n_species, fill=dispNhood)) + 
  geom_boxplot() + 
  facet_grid(cols=vars(mu2), rows=vars(maxLambda2)) +
  scale_fill_manual(values=cbp1[1:2], name='') + 
  xlab('% Habitat loss') + 
  ylab('# Species per subcommunity') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

# 6x6:
tiff(filename = 'figures/plot sensitivity analysis 20250605 A.tif', 
     width = 6, height = 6, units = 'in', res = 600)
p1
dev.off()



# in a barplot:
sel <- (df$mu != 3)
df <- df[sel,]
source('./src/summarySE.R')
sdf <- summarySE(df, measurevar='n_species',
                 groupvars=c('mu', 'f_loss', 'maxLambda', 'dispKernel'))
cbp1 <- c("#009E73", "#E69F00", "#999999")
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')

# First, show that there is no effect of the dispersal neighbourhood in 
# the simulations with the original parameter combinations:
sel <- (sdf$maxLambda == 1)&(sdf$dispKernel %in% c('Original', 'Exponential'))
sdf$dispNhood <- '21 x 21 cells'
sdf$dispNhood[sdf$dispKernel == 'Original'] <- '11 x 11 cells'

ggplot(sdf[sel,], aes(x=f_loss*100, y=n_species, fill=dispNhood)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_wrap(vars(mu2)) + 
  geom_errorbar(aes(ymin=n_species-sd, ymax=n_species+sd), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('# Species per subcommunity') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

m1 <- sdf[sel,]

# average % ancestors from elsewhere:
df$Pm <- (1- df$Pm)*100
sdf <- summarySE(df, measurevar='Pm',
                 groupvars=c('mu', 'f_loss', 'maxLambda', 'dispKernel'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')
sel <- (sdf$maxLambda == 1)&(sdf$dispKernel %in% c('Original', 'Exponential'))
sdf$dispNhood <- '21 x 21 cells'
sdf$dispNhood[sdf$dispKernel == 'Original'] <- '11 x 11 cells'

m2 <- sdf[sel,]

sdf$ymins <- sdf$Pm-sdf$sd
sdf$ymins[sdf$ymins < 0] <- 0
ggplot(sdf[sel,], aes(x=f_loss*100, y=Pm, fill=dispNhood)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_wrap(vars(mu2)) + 
  geom_errorbar(aes(ymin=ymins, ymax=Pm+sd), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('% Ancestors from elsewhere') + 
  ylim(c(0, 100)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

# average dissimilarity:
filename <- 'results/dissimilarity/dissimilarity_data_1.txt'
df3 <- read.table(filename, header = TRUE)
sel <- (df3$mu != 3)&(df3$f_hab_loss %in% c(0.25, 0.5, 0.75, 0.95))&(df3$dispersal == 'different')
df3 <- df3[sel,]
df3$maxLambda <- 1
df3$dispKernel <- 'Original'

filename <- 'results/dissimilarity/dissimilarity_data_2_maxLambda=1_dispKernel=exponential.txt'
df3a <- read.table(filename, header = TRUE)
df3a$maxLambda <- 1
df3a$dispKernel <- 'Exponential'
df3 <- rbind(df3, df3a)

names(df3)[names(df3) == 'f_hab_loss'] <- 'f_loss'

sdf <- summarySE(df3, measurevar='dissimilarity',
                 groupvars=c('mu', 'f_loss', 'maxLambda', 'dispKernel'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')
sdf$dispNhood <- '21 x 21 cells'
sdf$dispNhood[sdf$dispKernel == 'Original'] <- '11 x 11 cells'

m3 <- sdf

sdf$ymaxs <- sdf$dissimilarity + sdf$sd 
sdf$ymaxs[sdf$ymaxs > 1] <- 1
ggplot(sdf, aes(x=f_loss*100, y=dissimilarity, fill=dispNhood)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_wrap(vars(mu2)) + 
  geom_errorbar(aes(ymin=dissimilarity-sd, ymax=ymaxs), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  ylim(c(0, 1)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

# add dispersal capacity? 
start <- TRUE
for (i in c(0.25, 0.5, 0.75, 0.95)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m        <- read.table(filename, header = TRUE)
    m        <- m[m$mutation_rate == 0.0003,]
    
    if (start == TRUE){
      df <- m
      start <- FALSE
    } else {
      df    <- rbind(df, m)
    }
  }
}
df <- df[df$dispersal == 'different',]
df <- df[(df$sim_nr == 1)&(df$mu != 3),]
df$maxLambda <- 1
df$dispKernel <- 'Original'

for (i in c(0.25, 0.5, 0.75, 0.95)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, 
                    '_maxLambda=1_dispKernel=exponential.txt', sep='')
  if (file.exists(filename)){
    m        <- read.table(filename, header = TRUE)
    m$maxLambda <- 1
    m$dispKernel <- 'Exponential'
    df    <- rbind(df, m)
  }
}

# plot the average dispersal strategy:
df$y <- df$perc_indiv/10 * df$disp_cap

sdf <- summarySE(df, measurevar="y", 
                 groupvars=c("mu","f_loss", 'maxLambda', 'dispKernel'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')
sdf$dispNhood <- '21 x 21 cells'
sdf$dispNhood[sdf$dispKernel == 'Original'] <- '11 x 11 cells'

m4 <- sdf

# plot them together
names(m1)[names(m1) == 'n_species'] <- 'y'
names(m2)[names(m2) == 'Pm'] <- 'y'
names(m3)[names(m3) == 'dissimilarity'] <- 'y'

m1$ylab <- '# Species per subcommunity'
m2$ylab <- '% Ancestors from elsewhere'
m3$ylab <- 'Bray-Curtis dissimilarity'
m4$ylab <- 'Dispersal capacity (λ)'

sdf <- rbind(m1, m2, m3, m4)

sdf$ymins <- sdf$y - sdf$sd
sdf$ymins[sdf$ymins < 0] <- 0
sdf$ymaxs <- sdf$y + sdf$sd
sdf$ymaxs[(sdf$ymaxs > 1)&(sdf$ylab == 'Bray-Curtis dissimilarity')] <- 1
sdf$ylab <- factor(sdf$ylab, levels=unique(sdf$ylab))
p1 <- ggplot(sdf, aes(x=f_loss*100, y=y, fill=dispNhood)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_grid(cols=vars(mu2), rows=vars(ylab), scales='free_y', switch='y') + 
  geom_errorbar(aes(ymin=ymins, ymax=ymaxs), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
p1

pd <- position_dodge(3)
p2 <- ggplot(sdf, aes(x=f_loss*100, y=y, color=dispNhood)) + 
  geom_line(position=pd) + 
  geom_point(position=pd) + 
  scale_color_manual(values=cbp1, name='') + 
  facet_grid(cols=vars(mu2), rows=vars(ylab), scales='free_y', switch='y') + 
  geom_errorbar(aes(ymin=ymins, ymax=ymaxs), 
                width=0, position=pd) + 
  xlab('% Habitat loss') + 
  ylab('') + 
  scale_x_continuous(breaks=c(25,50, 75, 95)) +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
p2

# 6x8:
# tiff file 600 dpi:
tiff(filename = 'figures/plot sensitivity analysis 20250605 A.tif', 
     width = 6, height = 8, units = 'in', res = 600)
p2
dev.off()

###############################################################################
# second, plot the other parameter combinations:

ggplot(sdf, aes(x=f_loss*100, y=n_species, fill=maxLambda2)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_wrap(vars(mu2)) + 
  geom_errorbar(aes(ymin=n_species-sd, ymax=n_species+sd), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('# Species per subcommunity') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
  
m1 <- sdf

# average % ancestors from elsewhere:
df$Pm <- (1- df$Pm)*100
sdf <- summarySE(df, measurevar='Pm',
                 groupvars=c('mu', 'f_loss', 'maxLambda'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')

m2 <- sdf

sdf$ymins <- sdf$Pm-sdf$sd
sdf$ymins[sdf$ymins < 0] <- 0
ggplot(sdf, aes(x=f_loss*100, y=Pm, fill=maxLambda2)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_wrap(vars(mu2)) + 
  geom_errorbar(aes(ymin=ymins, ymax=Pm+sd), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('% Ancestors from elsewhere') + 
  ylim(c(0, 1)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

# average dissimilarity:
filename <- 'results/dissimilarity/dissimilarity_data_1.txt'
df3 <- read.table(filename, header = TRUE)
sel <- (df3$mu != 3)&(df3$f_hab_loss %in% c(0.25, 0.5, 0.75, 0.95))&(df3$dispersal == 'different')
df3 <- df3[sel,]
df3$maxLambda <- 1

filename <- 'results/dissimilarity/dissimilarity_data_1_maxLambda=0.5.txt'
df3b <- read.table(filename, header = TRUE)
df3b$maxLambda <- 0.5

filename <- 'results/dissimilarity/dissimilarity_data_1_maxLambda=1.5.txt'
df3c <- read.table(filename, header = TRUE)
df3c$maxLambda <- 1.5

df3 <- rbind(df3, df3b, df3c)
names(df3)[names(df3) == 'f_hab_loss'] <- 'f_loss'

sdf <- summarySE(df3, measurevar='dissimilarity',
                 groupvars=c('mu', 'f_loss', 'maxLambda'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')

m3 <- sdf

sdf$ymaxs <- sdf$dissimilarity + sdf$sd 
sdf$ymaxs[sdf$ymaxs > 1] <- 1
ggplot(sdf, aes(x=f_loss*100, y=dissimilarity, fill=maxLambda2)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_wrap(vars(mu2)) + 
  geom_errorbar(aes(ymin=dissimilarity-sd, ymax=ymaxs), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('Bray-Curtis dissimilarity') + 
  ylim(c(0, 1)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

# add dispersal capacity? 
start <- TRUE
for (i in c(0.25, 0.5, 0.75, 0.95)){
  filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m        <- read.table(filename, header = TRUE)
    m        <- m[m$mutation_rate == 0.0003,]
    
    if (start == TRUE){
      df <- m
      start <- FALSE
    } else {
      df    <- rbind(df, m)
    }
  }
}
df <- df[df$dispersal == 'different',]
df <- df[(df$sim_nr == 1)&(df$mu != 3),]
df$maxLambda <- 1

for (i in c(0.25, 0.5, 0.75, 0.95)){
  for (j in c(0.5, 1.5)){
    filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, 
                      '_maxLambda=',j,'.txt', sep='')
    if (file.exists(filename)){
      m        <- read.table(filename, header = TRUE)
      m$maxLambda <- j
      df    <- rbind(df, m)
    }
  }
}

# plot the average dispersal strategy:
df$y <- df$perc_indiv/10 * df$disp_cap

sdf <- summarySE(df, measurevar="y", 
                 groupvars=c("mu","f_loss", 'maxLambda'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')

m4 <- sdf

# plot them together
names(m1)[names(m1) == 'n_species'] <- 'y'
names(m2)[names(m2) == 'Pm'] <- 'y'
names(m3)[names(m3) == 'dissimilarity'] <- 'y'

m1$ylab <- '# Species per subcommunity'
m2$ylab <- '% Ancestors from elsewhere'
m3$ylab <- 'Bray-Curtis dissimilarity'
m4$ylab <- 'Dispersal capacity (λ)'

sdf <- rbind(m1, m2, m3, m4)

sdf$ymins <- sdf$y - sdf$sd
sdf$ymins[sdf$ymins < 0] <- 0
sdf$ymaxs <- sdf$y + sdf$sd
sdf$ymaxs[(sdf$ymaxs > 1)&(sdf$ylab == 'Bray-Curtis dissimilarity')] <- 1
sdf$ylab <- factor(sdf$ylab, levels=unique(sdf$ylab))
p1 <- ggplot(sdf, aes(x=f_loss*100, y=y, fill=maxLambda2)) + 
  geom_bar(stat='identity', position=position_dodge(10), 
           just=0.9, width=10, color='black') + 
  scale_fill_manual(values=cbp1, name='') + 
  facet_grid(cols=vars(mu2), rows=vars(ylab), scales='free_y', switch='y') + 
  geom_errorbar(aes(ymin=ymins, ymax=ymaxs), 
                width=0, position=position_dodge(10)) + 
  xlab('% Habitat loss') + 
  ylab('') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
p1

pd <- position_dodge(3)
p2 <- ggplot(sdf, aes(x=f_loss*100, y=y, color=maxLambda2)) + 
  geom_line(position=pd) + 
  geom_point(position=pd) + 
  scale_color_manual(values=cbp1, name='') + 
  facet_grid(cols=vars(mu2), rows=vars(ylab), scales='free_y', switch='y') + 
  geom_errorbar(aes(ymin=ymins, ymax=ymaxs), 
                width=0, position=pd) + 
  xlab('% Habitat loss') + 
  ylab('') + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
p2

# 6x8:
# tiff file 600 dpi:
tiff(filename = 'figures/plot sensitivity analysis 20250605.tif', 
     width = 6, height = 8, units = 'in', res = 600)
p2
dev.off()
