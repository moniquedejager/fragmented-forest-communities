library(ggplot2)
library(ggpubr)

filename <- 'results/subcommunity_data/fragmented_subcommunity_data_0.5.txt'
df <- read.table(filename, header = TRUE)

ix <- 1:length(df$mu)
group <- paste(df$mu, df$x, df$y, sep='-')
ix <- tapply(ix, group, max)
df2 <- df[ix,]

for (i in seq(0.1, 0.95, 0.05)){
  filename <- paste('results/subcommunity_data/fragmented_subcommunity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    
    ix <- 1:length(m$mu)
    group <- paste(m$mu, m$x, m$y, sep='-')
    ix <- tapply(ix, group, max)
    df <- rbind(df, m[ix,])
  }
}

#df$METE_fit[df$METE_fit < 0] <- 0
df$Pm <- 1 - df$Pm
#df$Shannon <- df$Shannon / log(df$n_species)

# with geom_raster:
group <- paste(df$f_loss, df$mu, sep='-')
df3   <- data.frame(frag     = vector(length = 0),
                    mu       = vector(length = 0),
                    fit      = vector(length = 0),
                    fit_dens = vector(length = 0),
                    nspecies = vector(length = 0),
                    nspec_dens = vector(length = 0),
                    Pm       = vector(length = 0),
                    Pm_dens  = vector(length = 0))

brks <- (0:22)/20 - 0.05

for (i in unique(group)){
  
  fit_hist      <- hist(df$METE_fit[(group == i)], breaks = (brks * 12 - 11.8))
  nspecies_hist <- hist(df$n_species[(group == i)], breaks = brks * 500)
  Pm_hist       <- hist(df$Pm[(group == i)], breaks = brks * 1)

  df3b   <- data.frame(frag    = df$f_loss[(group == i)][1],
                      mu       = df$mu[(group == i)][1],
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
                  y    = c(df3$fit, df3$nspecies,  
                           df3$Pm),
                  z    = c(df3$fit_dens, df3$nspec_dens,
                           df3$Pm_dens),
                  type = c(rep("METE's fit to SAD", length(df3$frag)),
                           rep('# Species per subcommunity', length(df3$frag)),
                           rep("% Ancestors from elsewhere", length(df3$frag))))

df4$mu2              <- "μ = 1"
df4$mu2[df4$mu == 3] <- "μ = 3"
df4$mu2[df4$mu == 5] <- "μ = 5"

df4$z2             <- df4$z
df4$z2[df4$z == 0] <- NA
df4 <- df4[!is.na(df4$z2),]

df4$clustering <- 'Random'
df4$clustering[df4$mu == 3] <- 'Fractal'
df4$clustering[df4$mu == 5] <- 'Clustered'

ggplot(df4, aes(x=frag, y=y, fill=z2/100)) + 
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
