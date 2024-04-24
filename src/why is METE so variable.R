# why are there such massive differences in METE fit within the same level of
# habitat loss?

library(ggplot2)
library(ggpubr)

filename <- 'results/simulation results 03-2024/subcommunity_data/subcommunity_data_1.txt'
df <- read.table(filename, header = TRUE)

for (i in 2:10){
  filename <- paste('results/simulation results 03-2024/subcommunity_data/subcommunity_data_', i, '.txt', sep='')
  if (file.exists(filename)){
    m  <- read.table(filename, header = TRUE)
    df <- rbind(df, m)
  }
}

df$dist_to_center <- sqrt((df$x - 10)^2 + (df$y - 10)^2)
sel <- (df$f_hab_loss == 0) #&(abs(df$x - 25) < 15)&(abs(df$y - 25) < 15)
ggplot(df[sel,], aes(x=as.factor(mu), y=METE_fit)) + 
  geom_violin()

ggplot(df[sel,], aes(x=n_species, y=METE_fit, color=as.factor(max_mutation))) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mutation_rate)) 

ggplot(df[sel,], aes(x=Pm, y=METE_fit, color=as.factor(max_mutation))) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mutation_rate)) 

ggplot(df[sel,], aes(x=Pm, y=n_species, color=as.factor(max_mutation))) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mutation_rate)) 

ggplot(df[sel,], aes(x=dist_to_center, y=METE_fit, 
                     color=as.factor(max_mutation))) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mutation_rate)) 

ggplot(df[sel,], aes(x=dist_to_center, y=n_species,
                     color=as.factor(max_mutation))) + 
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(mutation_rate)) 

df$rx <- round(df$x, -1)
ggplot(df[sel,], aes(x=dist_to_center, y=n_species, color=as.factor(mu))) + 
  geom_point() +
  geom_smooth()
df$ry <- round(df$y, -1)

ggplot(df[sel,], aes(x=dist_to_center, y=n_species, color=as.factor(rx))) + 
  geom_point() +
  geom_smooth(se=FALSE, method='gam') + facet_grid(rows = vars(sim_nr), cols=vars(mu)) + 
  ylim(c(0, 30))

ggplot(df[sel,], aes(x=dist_to_center, y=METE_fit, color=as.factor(rx))) + 
  geom_point() +
  facet_grid(rows = vars(sim_nr), cols=vars(mu)) 

ggplot(df[sel,], aes(x=dist_to_center, y=n_species)) + 
  geom_point() +
  facet_grid(rows = vars(sim_nr), cols=vars(mu)) 

# at mu = 4, there are no patches near the center, only patches close to the 
# edge. This results in the significant effect of mu on METE fit! 

# center communities have fewer species and more dispersal between patches
# (as the local dispersing species have been outcompeted and are not 
# replenished by the metacommunity). 

# A different method to deal with this problem is to simulate the subcommunities
# without a bordering metacommunity, and use random mutations to obtain new 
# species (as other models do). This would mean that I have to redo the 
# simulations, after adjusting the model, of course! 
# boundary conditions should then become continuous. 

# what should be the mutation rate? 






# The proportion of t-50 ancestors from outside of the current subcommunity
# depends on the dispersal capacity of the species. Only few species are able
# to always disperse outside the parental subcommunity (Pm = 1). Hence, when Pm
# is high, the number of species is low, as the more locally dispersing species
# are not present.

df$dist <- 'Local dispersal'
df$dist[df$m_dist > 3] <- 'Intermediate dispersal'
df$dist[df$m_dist > 5] <- 'Far dispersal'
df$dist <- factor(df$dist, levels = unique(df$dist)[3:1])

ggplot(df[sel,], aes(x=n_species, y=METE_fit)) + 
  geom_point() +
  facet_wrap(vars(dist)) + 
  xlab('Number of species') + 
  ylab('METE fit to SAD') + 
  scale_x_continuous(breaks=c(150,200, 250, 300, 350), limits=c(145, 355)) + 
  scale_y_continuous(limits=c(-0.5, 1)) + 
  theme_bw() + 
  theme(strip.placement = "outside", 
        strip.background = element_blank())


ggplot(df[sel,], aes(x=Shannon, y=METE_fit)) + 
  geom_point() +
  facet_wrap(vars(dist)) + 
  xlab('Shannon Index') + 
  ylab('METE fit to SAD') + 
  scale_y_continuous(limits=c(-0.5, 1)) + 
  theme_bw() + 
  theme(strip.placement = "outside", 
        strip.background = element_blank())

df$evenness <- df$Shannon / log(df$n_species)
ggplot(df[sel,], aes(x=evenness, y=METE_fit)) + 
  geom_point() +
  facet_wrap(vars(dist)) + 
  xlab('Evenness') + 
  ylab('METE fit to SAD') + 
  scale_y_continuous(limits=c(-0.5, 1)) + 
  theme_bw() + 
  theme(strip.placement = "outside", 
        strip.background = element_blank())

df$logit_evenness <- log(df$evenness/(1 - df$evenness))
ggplot(df[sel,], aes(x=logit_evenness, y=METE_fit)) + 
  geom_point() +
  facet_wrap(vars(dist)) + 
  xlab('logit Evenness') + 
  ylab('METE fit to SAD') + 
  scale_y_continuous(limits=c(-0.5, 1)) + 
  theme_bw() + 
  theme(strip.placement = "outside", 
        strip.background = element_blank())

summary(lm(METE_fit~n_species*mu,data=df[sel,]))
# mu zou niet uit moeten maken, aangezien het landschap nog ongefragmenteerd is!

df$Pm_logit <- log(df$Pm/(1 - df$Pm))
ggplot(df[sel,], aes(x=Pm_logit, y=METE_fit)) + 
  geom_point() + 
  facet_wrap(vars(mu))

ggplot(df[sel,], aes(x=Pm_logit, y=n_species)) + 
  geom_point() + 
  facet_wrap(vars(mu))

ggplot(df[sel,], aes(x=m_dist, y=METE_fit)) + 
  geom_point()

ggplot(df[sel,], aes(x=m_dist, y=n_species)) + 
  geom_point()

ggplot(df[sel,], aes(x=m_dist, y=Pm)) + 
  geom_point()


# let's examine the bad fits more closely...
df[(df$METE_fit < -0.25)&(df$f_hab_loss == 0),]

df2 <- df[(df$mu == 1)&(df$f_hab_loss == 0)&(df$sim_nr == 3),]

filename <- './results/simulation results 03-2024/community_composition/1_0_3.txt'
m <-as.matrix(read.table(filename))

x <- rep(1:50, 50)
y <- sort(x)
group <- paste(x, y, sep='-')
group2 <- paste(df2$x, df2$y, sep='-')

i <- 1:2500

m2 <- m[,i[(x == 11)&(y == 8)]]
SAD <- sort(tapply(m2, as.factor(m2), length), decreasing = TRUE)
plot(SAD)

m3 <- m[,i[(x == 1)&(y == 19)]]
SAD2 <- sort(tapply(m3, as.factor(m3), length), decreasing = TRUE)
lines(SAD2)

mean(m2)
mean(m3)


m4 <- m[,group %in% group2]
group3 <- group[group %in% group2]
m_spec_nr <- colMeans(m )
df2$m_spec_nr <- 0
for (i in group2){
  print(paste(i, ', ', m_spec_nr[group3 == i]))
  df2$m_spec_nr[group2 == i] <- m_spec_nr[group3 == i]
}

ggplot(df2, aes(x=m_spec_nr, y=METE_fit)) + geom_point()
ggplot(df2, aes(x=m_spec_nr, y=m_dist)) + geom_point()

ggplot(df2, aes(x=x, y=y, fill=m_dist)) + geom_raster()



