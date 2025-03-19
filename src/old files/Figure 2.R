library(ggplot2)
library(ggpubr)

filename <- 'results/simulation results 04-2024/subcommunity_data/subcommunity_data_1.txt'
df <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('results/simulation results 04-2024/subcommunity_data/subcommunity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    df <- rbind(df, m)
  }
}

df$METE_fit[df$METE_fit < 0] <- 0
#df$Shannon <- df$Shannon / log(df$n_species)

# with geom_raster:
group <- paste(df$f_hab_loss, df$mu, sep='-')
df3   <- data.frame(frag     = vector(length = 0),
                    mu       = vector(length = 0),
                    fit      = vector(length = 0),
                    fit_dens = vector(length = 0),
                    nspecies = vector(length = 0),
                    nspec_dens = vector(length = 0),
                    m_dist   = vector(length = 0),
                    m_dist_dens = vector(length = 0),
                    Pm       = vector(length = 0),
                    Pm_dens  = vector(length = 0),
                    H        = vector(length = 0),
                    H_dens   = vector(length = 0))

brks <- (0:22)/20 - 0.05

for (i in unique(group)){
  
  fit_hist      <- hist(df$METE_fit[(group == i)], breaks = (brks * 1))
  nspecies_hist <- hist(df$n_species[(group == i)], breaks = brks * 260)
  m_dist_hist   <- hist(df$m_dist[(group == i)], breaks = brks * 15)
  Pm_hist       <- hist(df$Pm[(group == i)], breaks = brks * 1)
  H_hist        <- hist(df$Shannon[(group == i)], breaks = brks * 7)

  df3b   <- data.frame(frag    = df$f_hab_loss[(group == i)][1],
                      mu       = df$mu[(group == i)][1],
                      fit      = fit_hist$mids,
                      fit_dens = fit_hist$counts / 
                        sum(fit_hist$counts) * 100,
                      nspecies = nspecies_hist$mids,
                      nspec_dens = nspecies_hist$counts / 
                        sum(nspecies_hist$counts) * 100,
                      m_dist   = m_dist_hist$mids,
                      m_dist_dens = m_dist_hist$counts / 
                        sum(m_dist_hist$counts) * 100,
                      Pm       = Pm_hist$mids,
                      Pm_dens  = Pm_hist$counts / 
                        sum(Pm_hist$counts) * 100,
                      H        = H_hist$mids,
                      H_dens   = H_hist$counts / 
                        sum(H_hist$counts) * 100)
  df3 <- rbind(df3, df3b)
}

# all plots together:
df4 <- data.frame(frag = rep(df3$frag, 5),
                  mu   = rep(df3$mu, 5),
                  y    = c(df3$fit, df3$nspecies, df3$m_dist, 
                           df3$Pm, df3$H),
                  z    = c(df3$fit_dens, df3$nspec_dens, df3$m_dist_dens, 
                           df3$Pm_dens, df3$H_dens),
                  type = c(rep("METE's fit to SAD", length(df3$frag)),
                           rep('# Species per subcommunity', length(df3$frag)),
                           rep("Average distance from ancestor",  length(df3$frag)),
                           rep("% Ancestors from elsewhere", length(df3$frag)),
                           rep("Evenness",length(df3$frag))))

df4$mu2              <- "μ = 1"
df4$mu2[df4$mu == 3] <- "μ = 3"
df4$mu2[df4$mu == 5] <- "μ = 5"

df4$z2             <- df4$z
df4$z2[df4$z == 0] <- NA
sel <- (df4$type != 'Average distance from ancestor')&(!is.na(df4$z2))
ggplot(df4[sel,], aes(x=frag, y=y, fill=z2/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(mu2), 
             rows = vars(type), 
             scales = 'free', 
             switch = 'y') + 
  scale_fill_gradientn(colors = c('yellow','seagreen', 
                                  'darkcyan', 'midnightblue', 'black', 'black'),
                       labels = scales::percent_format(),
                       na.value = '',
                       #trans = 'log10',
                       name="% subcommunities",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('% Habitat loss') + 
  ylab("") + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 15)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())


df5 <- df4[sel,]
df5 <- df5[(df5$mu != 2)&(df5$type != 'Evenness'),]
df5$clustering <- 'Random'
df5$clustering[df5$mu == 3] <- 'Fractal'
df5$clustering[df5$mu == 5] <- 'Clustered'

ggplot(df5, aes(x=frag, y=y, fill=z2/100)) + 
  geom_raster() + 
  facet_grid(cols = vars(clustering), 
             rows = vars(type), 
             scales = 'free', 
             switch = 'y') + 
  scale_fill_gradientn(colors = c('yellow','seagreen', 
                                  'darkcyan', 'midnightblue', 'black', 'black'),
                       labels = scales::percent_format(),
                       na.value = '',
                       #trans = 'log10',
                       name="% subcommunities",
                       limits = c(0, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('% Habitat loss') + 
  ylab("") + 
  theme_bw() + 
  guides(fill = guide_colorbar(barwidth = 15)) + 
  theme(legend.position = 'top',
        strip.placement = "outside", 
        strip.background = element_blank())



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
