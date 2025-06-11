library(ggplot2)
library(ggpubr)

filename <- 'results/subcommunity_data/fragmented_subcommunity_data_0.25.txt'
df <- read.table(filename, header = TRUE)
df$maxLambda <- 1
df$dispKernel <- 'Original'
df <- df[df$mu == 0,]

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
df$dispKernel[df$dispKernel == 'gaussian'] <- 'Gaussian'
df$dispKernel[df$dispKernel == 'pareto'] <- 'Pareto'

df$type <- paste(df$dispKernel, ' (0 < λ ≤ ',df$maxLambda, ')',sep='')

mod <- lm(n_species ~ f_loss * mu * type, data=df)
summary(mod)

mod <- lm(Pm ~ f_loss * mu * type, data=df)
summary(mod)


ggplot(df, aes(x=factor(f_loss), y=n_species, color=type)) + 
  geom_violin() + 
  scale_y_continuous(trans='log10') + 
  facet_wrap(vars(mu))
ggplot(df, aes(x=factor(f_loss), y=Pm, color=type)) + 
  geom_boxplot() + 
  facet_wrap(vars(mu))


# average dissimilarity:
filename <- 'results/dissimilarity/dissimilarity_data_2_maxLambda=1_dispKernel=exponential.txt'
df3 <- read.table(filename, header = TRUE)
df3$maxLambda <- 1
df3$dispKernel <- 'Exponential'
df3 <- df3[df3$mu == 0,]

for (i in c(0.5, 1, 1.5)){
  for (j in c('exponential', 'gaussian', 'pareto')){
    for (k in 1:5){ #1:5
      filename <- paste('results/dissimilarity/dissimilarity_data_', k, 
                        '_maxLambda=', i, '_dispKernel=', j, '.txt', sep='')
      if (file.exists(filename)){
        df3b <- read.table(filename, header = TRUE)
        df3b$maxLambda <- i
        df3b$dispKernel <- j
        df3 <- rbind(df3, df3b)
      }
    }
  }
}
df3$dispKernel[df3$dispKernel == 'exponential'] <- 'Exponential'
df3$dispKernel[df3$dispKernel == 'gaussian'] <- 'Gaussian'
df3$dispKernel[df3$dispKernel == 'pareto'] <- 'Pareto'
df3$type <- paste(df3$dispKernel, ' (0 < λ ≤ ',df3$maxLambda, ')',sep='')

names(df3)[names(df3) == 'f_hab_loss'] <- 'f_loss'

ggplot(df3, aes(x=factor(f_loss), y=dissimilarity, color=type)) + 
  geom_boxplot() + 
  facet_wrap(vars(mu))

mod <- lm(dissimilarity~f_loss*mu*type, data=df3)
summary(mod)

# add dispersal capacity? 
filename <- './results/dispersal_capacity/fragmented_dispersal capacity_data_0.25_maxLambda=1_dispKernel=exponential.txt'
df <- read.table(filename, header=TRUE)
df$maxLambda <- 1
df$dispKernel <- 'Original'
df <- df[df$mu == 0,]

for (i in c(0.25, 0.5, 0.75, 0.95)){
  for (j in c('exponential', 'gaussian', 'pareto')){
    for (k in c(0.5, 1, 1.5)){
      filename <- paste('./results/dispersal_capacity/fragmented_dispersal capacity_data_', i, 
                        '_maxLambda=', k, '_dispKernel=', j, '.txt', sep='')
      if (file.exists(filename)){
        df3b <- read.table(filename, header = TRUE)
        df3b$maxLambda <- k
        df3b$dispKernel <- j
        df <- rbind(df, df3b)
      }
    }
  }
}
df$dispKernel[df$dispKernel == 'exponential'] <- 'Exponential'
df$dispKernel[df$dispKernel == 'gaussian'] <- 'Gaussian'
df$dispKernel[df$dispKernel == 'pareto'] <- 'Pareto'

# plot the average dispersal strategy:
df$y <- df$perc_indiv/10 * df$disp_cap

sdf <- summarySE(df, measurevar="y", 
                 groupvars=c("mu","f_loss", 'maxLambda', 'dispKernel'))
sdf$mu2 <- 'Random habitat destruction'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$maxLambda2 <- paste('0 < λ ≤ ', sdf$maxLambda, sep='')
sdf$type <- paste(sdf$dispKernel, ' (0 < λ ≤ ',sdf$maxLambda, ')',sep='')

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

pd <- position_dodge(7)
cbp1 <- c( "#307314", "#009E73",  "#999999","#E69F00", "#6E2C00")
p2 <- ggplot(sdf, aes(x=f_loss*100, y=y, color=type)) + 
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
        strip.background = element_blank()) + 
  guides(color = guide_legend(ncol = 2)) 
p2

# 6x8:
# tiff file 600 dpi:
tiff(filename = 'figures/plot sensitivity analysis 20250605 B sim2.tif', 
     width = 6, height = 9, units = 'in', res = 600)
p2
dev.off()
