library(ggplot2)
library(ggpubr)

filename <- 'results/subcommunity_data/fragmented_subcommunity_data_0.05.txt'
df <- read.table(filename, header = TRUE)

for (i in seq(0, 0.95, 0.05)){
  filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    df <- rbind(df, m)
    
    #df <- m[m$dispersal == 'different',]
    #write.table(df, filename, append = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

df$METE_fit[df$METE_fit < 0] <- 0
#df$Shannon <- df$Shannon / log(df$n_species)

# with geom_raster:
group <- paste(df$f_loss, df$mu, df$dispersal, sep='-')
df3   <- data.frame(frag     = vector(length = 0),
                    mu       = vector(length = 0),
                    disp     = vector(length = 0),
                    fit      = vector(length = 0),
                    fit_dens = vector(length = 0),
                    nspecies = vector(length = 0),
                    nspec_dens = vector(length = 0),
                    Pm       = vector(length = 0),
                    Pm_dens  = vector(length = 0))

brks <- (0:22)/20 - 0.05

for (i in unique(group)){
  
  fit_hist      <- hist(df$METE_fit[(group == i)], breaks = (brks * 1))
  if (df$dispersal[(group == i)][1] == 'similar'){
    nspecies_hist <- hist(df$n_species[(group == i)], breaks = brks * 450)
  } else {
    nspecies_hist <- hist(df$n_species[(group == i)], breaks = brks * 100)
  }
  
  Pm_hist       <- hist(df$Pm[(group == i)], breaks = brks * 1)
 
  df3b   <- data.frame(frag    = df$f_loss[(group == i)][1],
                      mu       = df$mu[(group == i)][1],
                      disp     = df$dispersal[(group == i)][1],
                      fit      = fit_hist$mids,
                      fit_dens = fit_hist$counts / 
                        sum(fit_hist$counts) * 100,
                      nspecies = nspecies_hist$mids,
                      nspec_dens = nspecies_hist$counts / 
                        sum(nspecies_hist$counts) * 100,
                      Pm       = Pm_hist$mids,
                      Pm_dens  = Pm_hist$counts / 
                        sum(Pm_hist$counts) * 100)
  df3 <- rbind(df3, df3b)
}

# all plots together:
df4 <- data.frame(frag = rep(df3$frag, 3),
                  mu   = rep(df3$mu, 3),
                  disp = rep(df3$disp, 3),
                  y    = c(df3$fit, df3$nspecies, 1-df3$Pm),
                  z    = c(df3$fit_dens, df3$nspec_dens, df3$Pm_dens),
                  type = c(rep("METE's fit to SAD", length(df3$frag)),
                           rep('# Species per subcommunity', length(df3$frag)),
                           rep("% Ancestors from elsewhere", length(df3$frag))))

df4$z2             <- df4$z
df4$z2[df4$z == 0] <- NA

df5 <- df4[df4$frag > 0,]
df5 <- df5[(df5$mu != 2)&(df5$type != 'Evenness'),]
df5$clustering <- 'Random'
df5$clustering[df5$mu == 3] <- 'Fractal'
df5$clustering[df5$mu == 5] <- 'Clustered'

sel <- df5$disp == 'similar'
p1 <- ggplot(df5[sel,], aes(x=frag*100, y=y, fill=z2/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering), 
             rows = vars(type), 
             scales = 'free', 
             switch = 'y') + 
  scale_fill_gradientn(colors = c('yellow','seagreen', 
                                  'darkcyan', 'midnightblue', 'black', 'black'),
                       labels = scales::percent_format(),
                       na.value = 'transparent',
                       #trans = 'log10',
                       name="% subcommunities",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('% Habitat loss') + 
  ylab("") + 
  ggtitle("A. Same dispersal strategy") + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 15)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
p1

sel <- (df5$disp == 'different')
p2 <- ggplot(df5[sel,], aes(x=frag*100, y=y, fill=z2/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering), 
             rows = vars(type), 
             scales = 'free', 
             switch = 'y') + 
  scale_fill_gradientn(colors = c('yellow','seagreen', 
                                  'darkcyan', 'midnightblue', 'black', 'black'),
                       labels = scales::percent_format(),
                       na.value = 'transparent',
                       #trans = 'log10',
                       name="% subcommunities",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('% Habitat loss') + 
  ylab("") + 
  ggtitle("B. Different dispersal strategies") + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 15)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
p2

ggarrange(p1, p2, align='hv', common.legend = TRUE)

df5$disp2 <- 'Different dispersal'
df5$disp2[df5$disp == 'similar'] <- 'Same dispersal'
df5$disp2 <- factor(df5$disp2, levels=unique(df5$disp2)[2:1])
sel <- df5$type == '# Species per subcommunity'
sel <- df5$type == '% Ancestors from elsewhere'
sel <- df5$type == "METE's fit to SAD"
ggplot(df5[sel,], aes(x=frag*100, y=y, fill=z2/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering), 
             rows = vars(disp2), 
             scales = 'free') + 
  scale_fill_gradientn(colors = c('yellow','seagreen', 
                                  'darkcyan', 'midnightblue', 'black', 'black'),
                       labels = scales::percent_format(),
                       na.value = 'transparent',
                       #trans = 'log10',
                       name="% subcommunities",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('% Habitat loss') + 
  ylab(df5$type[sel][1]) + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 15)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

##### use averages instead:
source('./src/summarySE.R')

# average number of species per subcommunity:
sdf <- summarySE(df, measurevar='n_species',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

pd <- position_dodge(1) 
ggplot(sdf, aes(x=f_loss*100, y=n_species, color=mu2))+
  geom_errorbar(aes(ymin=n_species-ci, ymax=n_species+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  facet_grid(rows=vars(disp_type), scales='free_y') + 
  xlab('% Habitat loss') + 
  ylab('Average no. species per subcommunity') +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())

# average % ancestors from elsewhere:
sdf <- summarySE(df, measurevar='Pm',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

pd <- position_dodge(1) 
ggplot(sdf, aes(x=f_loss*100, y=1 - Pm, color=mu2))+
  geom_errorbar(aes(ymin=(1 - Pm)-ci, ymax=(1 - Pm)+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  facet_grid(rows=vars(disp_type), scales='free_y') + 
  xlab('% Habitat loss') + 
  ylab('Average % ancestors from elsewhere') +
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())

# average METE fit to SAD per subcommunity:
sdf <- summarySE(df, measurevar='METE_fit',
                 groupvars=c('mu', 'f_loss', 'dispersal'))
sdf$mu2 <- 'Random'
sdf$mu2[sdf$mu == 3] <- 'Fractal'
sdf$mu2[sdf$mu == 5] <- 'Clustered'
sdf$disp_type <- 'Same dispersal'
sdf$disp_type[sdf$dispersal == 'different'] <- 'Different dispersal'
sdf$disp_type <- factor(sdf$disp_type, levels = c('Same dispersal', 'Different dispersal'))

pd <- position_dodge(1) 
ggplot(sdf, aes(x=f_loss*100, y=METE_fit, color=mu2))+
  geom_errorbar(aes(ymin=METE_fit-ci, ymax=METE_fit+ci), width=0, position=pd) + 
  geom_line(position=pd) +
  geom_point(position=pd) + 
  facet_grid(rows=vars(disp_type)) + 
  xlab('% Habitat loss') + 
  ylab('METE fit to SAD') +
  ylim(c(0, 1)) + 
  theme_bw() + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank(),
        legend.title=element_blank())




# To calculate differences between clustering levels:
sel3 <- (df$mu == 5)
tapply((df$n_species[sel3] > 10), as.factor(df$f_hab_loss[sel3]), mean)

sel3 <- (df$mu == 1.1)
tapply((df$Pm[sel3]), as.factor(df$f_hab_loss[sel3]), mean)

sel3 <- (df$mu == 4)
tapply((df$m_dist[sel3]), as.factor(df$f_hab_loss[sel3]), mean)

sel3 <- (df$mu == 4)
tapply((df$METE_fit[sel3] == 0), as.factor(df$f_hab_loss[sel3]), mean)

sel3 <- (sel)&(df$mu == 4)&(df$METE_fit != 0)
tapply((df$METE_fit[sel3] ), as.factor(df$f_hab_loss[sel3]), mean)

sel3 <- (df$METE_fit != 0)
ggplot(df[sel3,], aes(x=Frag, y=fit, color=as.factor(mu))) + 
  geom_smooth()

sel3 <- (df2$mu == 1.1)
tapply((df2$Shannon[sel3]), as.factor(df2$f_hab_loss[sel3]), mean)

sel3 <- (df2$mu == 4)
tapply((df2$Shannon[sel3] <= 1), as.factor(df2$f_hab_loss[sel3]), mean)


# all together
group  <- paste(df$f_hab_loss, df$mu, sep='-')
df5    <- data.frame(frag  = tapply(df$f_hab_loss, group, mean),
                    mu    = tapply(df$mu, group, mean),
                    f_n10 = tapply((df$n_species > 10), group, mean),
                    f_fit = tapply((df$METE_fit > 0.5), group, mean),
                    f_H   = tapply((df$Shannon > 1), group, mean),
                    f_Pm  = tapply(df$Pm, group, mean))

n <- length(df5$frag)
df6 <- data.frame(frag = rep(df5$frag, 4),
                  mu   = rep(df5$mu, 4),
                  z    = c(df5$f_n10, df5$f_Pm, df5$f_fit, df5$f_H),
                  type = c(rep('% Subcommunities with > 10 species', n),
                           rep('% Ancestors from elsewhere', n),
                           rep('% Subcommunities with Shannon index > 1', n), 
                           rep("% Subcommunities with METE fit > 0.5", n)))

df6$type <- factor(df6$type, labels = unique(df6$type))

ggplot(df6, aes(x=frag, y=z*100, color = as.factor(mu))) + 
  geom_point(aes(shape = as.factor(mu))) + 
  geom_smooth(aes(linetype = as.factor(mu)), se = FALSE) + 
  xlab('% Habitat loss') + 
  ylab('') + 
  ylim(c(0,100)) + 
  facet_wrap(vars(type), ncol=2, nrow=2) + 
  theme_bw() + 
  guides(linetype = 'none',
         shape    = 'none',
         color    = guide_legend(title = "Clustering parameter μ", 
                                 override.aes = list(shape = c(16,17, 15)))) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())

df5 <- df6[(df6$mu != 2)&(df6$type != '% Subcommunities with Shannon index > 1'),]
df5$clustering <- 'Random'
df5$clustering[df5$mu == 3] <- 'Fractal'
df5$clustering[df5$mu == 5] <- 'Clustered'

ggplot(df5, aes(x=frag, y=z*100, color = clustering)) + 
  geom_point(aes(shape = as.factor(mu))) + 
  geom_smooth(aes(linetype = as.factor(mu)), se = FALSE) + 
  xlab('% Habitat loss') + 
  ylab('') + 
  ylim(c(0,100)) + 
  facet_wrap(vars(type), ncol=3, nrow=1) + 
  theme_bw() + 
  guides(linetype = 'none',
         shape    = 'none',
         color    = guide_legend(title = "", 
                                 override.aes = list(shape = c(16,17, 15)))) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())
