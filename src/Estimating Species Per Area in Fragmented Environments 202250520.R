
library(ggplot2)
library(ggpubr)
library(meteR)

filename <- 'Fragmented-forest-communities/x64/Release/results/SARinPristineLandscape.txt'
df1 <- read.table(filename, header=TRUE)

ggplot(df1, aes(x=A, y=S, color=simnr)) + geom_point()

params_continuous <- data.frame(A0 = max(df1$A),
                                N0 = max(df1$A) * 1000,
                                S0 = max(df1$S))

sar1 <- meteSAR(Amin=1, 
                A0=params_continuous$A0, 
                S0=params_continuous$S0, 
                N0=params_continuous$N0) 

plot(sar1$pred$A, sar1$pred$S, ylim=c(0, 5000))

sar2 <- meteSAR(Amin=1, 
                A0=2025, 
                S0=4000, 
                N0=2025000) 

points(sar2$pred$A, sar2$pred$S, col='blue', pch=16)


sel <- (df2$simnr ==4)&(df2$clustering == 2.5)&(df2$f_loss == 0.9)
ggplot(df2[sel,], aes(x=A, y=s)) + geom_point()



sar3 <- meteSAR(Amin=1, 
                A0=2025, 
                S0=2000, 
                N0=3200) 

plot(df2a$A, df2a$S, col='orange', pch=16, ylim=c(0, 5000))
points(sar3$pred$A, sar3$pred$S, col='red', pch=16)


# voor een heleboel combinaties van S0 en N0 sar's maken en dan tegen de 
# sar's van de simulaties fitten om zo te zien hoe S0 en N0 aangepast moeten 
# worden per moransI en habitat cover:

# range N0: 2000 - 2,000,000
# range S0: 2000 - 10,000

params <- expand.grid(N0 = c(seq(2000, 10000, 100), seq(10000, 100000, 1000),
                             seq(100000, 2000000, 100000)),
                      S0 = seq(2000, 10000, 100))
params <- params[params$N0 > params$S0,]

sars <- data.frame(N0 = 0,
                   S0 = 0, 
                   A  = 0,
                   S  = 0)

for (i in 987:nrow(params)){
  sar <- meteSAR(Amin=1, 
                  A0=2025, 
                  S0=params$S0[i], 
                  N0=params$N0[i]) 
  sars_add <- data.frame(N0 = params$N0[i],
                         S0 = params$S0[i],
                         A  = sar$pred$A,
                         S  = sar$pred$S)
  sars <- rbind(sars, sars_add)
}

# als dit klaar is, sars wegschrijven zodat je het niet nog eens moet uitrekenen!! 
write.table(sars, 'Fragmented-forest-communities/x64/Release/results/SARs.txt', row.names = FALSE, col.names = TRUE)

sel <- (sars$N0 %in% c(3000, 10000, 50000, 100000, 500000, 2000000))&
  (sars$S0 %in% c(2000, 4000, 6000, 8000, 10000))

ggplot(sars[sel,], aes(x=A, y=S, color=factor(S0))) + 
  geom_line() + 
  facet_wrap(vars(N0))

# record the best fit (N0, S0, and fit) for each simulation, together with the 
# simulation parameters (f_loss, simnr, clustering, cover, and moransI)
df2 <- data.frame(simnr = 0, f_loss = 0, clustering = 0, 
                  cover = 0, moransI = 0, 
                  N0 = 0, S0 = 0, fit = 0)

for (clust in c('1.00', '1.50', '2.01','2.50', '3.00', '3.50', '4.00', '4.50', '5.00')){
  for (loss in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)){
    for (simnr in 1:5){
      filename <- paste('composition', clust,'_', simnr, 
                        '.00_', loss, '0_5500.00_2.00_6.50.txt', sep='')
      
      #filename <- 'composition1.00_3.00_0.10_5500.00_2.00_6.50.txt'
      
      # derive the simulation number, clustering, and habitat loss from the filename:
      dat        <- gsub("composition", "", filename)
      dat        <- unlist(strsplit(dat, '_'))
      
      clustering <- as.numeric(dat[1])
      sim_nr     <- as.numeric(dat[2])
      f_loss     <- as.numeric(dat[3])
      
      filename2 <- paste('Fragmented-forest-communities/x64/Release/community composition/', filename, sep='')
      df        <- read.table(filename2)
      names(df) <- c('x', 'y', 'species', 'n')
      df        <- df[df$x != -1,]
      df$patch  <- paste(df$x, df$y, sep='-')
      uPatch    <- unique(df$patch)

      #####################################################################
      # Calculating Moran's I:
      h <- rep(0, 2025)
      h[patchID %in% uPatch] <- 1
      
      # create vectors of the habitat values of each neighbor:
      zmean      <- mean(h)
      d0         <- h - zmean
      
      x2         <- x - 1
      x2[x2 < 1] <- n
      ID2        <- (y - 1)*n + x2
      n1         <- h[ID2]
      d1         <- n1 - zmean
      
      x2         <- x + 1
      x2[x2 > n] <- 1
      ID2        <- (y - 1)*n + x2
      n2         <- h[ID2]
      d2         <- n2 - zmean
      
      y2         <- y - 1
      y2[y2 < 1] <- n
      ID2        <- (y2 - 1)*n + x
      n3         <- h[ID2]
      d3         <- n3 - zmean
      
      y2         <- y + 1
      y2[y2 > n] <- 1
      ID2        <- (y2 - 1)*n + x
      n4         <- h[ID2]
      d4         <- n4 - zmean
      
      I <- 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2)) 
      #######################################################################
      
      # we need to add the zero's back in, if we want to do a proper job... 
      areas <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
      n     <- 45
      x     <- rep(1:n, n)
      y     <- sort(x)
      patchID <- paste(y-1, x-1, sep='-')
      df_add <- data.frame(x = x,
                           y = y,
                           species = 0,
                           n = 0,
                           patch = patchID)
      
      df <- rbind(df, df_add)

      uPatch <- unique(df$patch)
      S <- vector(length=0)
      for (i in 1:2025){
       S <- c(S, length(unique(df$species[df$patch == uPatch[i]])))  
      }
      mean(S)

      df2a <- data.frame(S = mean(S), 
                         A = 1)
      for (j in areas){
        S <- vector(length=0)
        for (i in 1:100){
          patches <- sample(uPatch, j, replace = FALSE)
          S <- c(S, length(unique(df$species[df$patch %in% patches])))
        }
        mean(S)
        df2a <- rbind(df2a, list(mean(S), j))
      }

      #ggplot(df2a, aes(x=A, y=S)) + geom_point()

      # which sar fits best? Calculate Goodness-of-fit per parameter combination:
      fit_data <- data.frame(S0 = 0,
                             N0 = 0,
                             fit = 0)
      
      group  <- paste(sars$N0, sars$S0, sep='-')
      ugroup <- unique(group)
      S_obs  <- df2a$S
      
      for (i in 1:length(ugroup)){
        sel    <- group == ugroup[i]
        S_pred <- sars$S[sel]
        fit    <- 1 - sum((S_obs - S_pred)^2)/sum((S_obs - mean(S_pred))^2)
        
        fit_data <- rbind(fit_data, list(sars$S0[sel][1], sars$N0[sel][1], fit))
      }
      fit_data <- fit_data[fit_data$S0 > 0,]

      #ggplot(fit_data, aes(x=N0, y=fit, color=factor(S0))) + geom_line() +
        scale_x_continuous(trans='log10')

      best_fit <- fit_data[fit_data$fit == max(fit_data$fit),]
      best_fit

      df2_add <- data.frame(simnr = sim_nr, 
                            f_loss = f_loss, 
                            clustering = clustering, 
                            cover = mean(h), 
                            moransI = I, 
                            N0 = best_fit$N0/1000, # keep in mind that I did this!!  
                            S0 = best_fit$S0, 
                            fit = best_fit$fit)
      df2 <- rbind(df2, df2_add)
    }
  }
}
df2 <- df2[df2$simnr > 0,]

write.table(df2, 'Fragmented-forest-communities/x64/Release/results/best_fits.txt', row.names = FALSE, col.names = TRUE)

ggplot(df2, aes(x=moransI, y=N0, color=S0)) + 
  geom_point() + 
  facet_wrap(vars(round(cover, 1)))

ggplot(df2, aes(x=moransI, y=S0, color=N0)) + 
  geom_point() + 
  facet_wrap(vars(round(cover, 1)), scales='free_y')

ggplot(df2, aes(x=moransI, y=fit, color=N0)) + 
  geom_point() + 
  facet_wrap(vars(round(cover, 1))) 

p1 <- ggplot(df2, aes(x=round(cover,1), y=round(moransI, 1), fill=N0)) + 
  geom_tile() + 
  scale_fill_viridis_c(name='N0 (x 1000)') + 
  xlab('% Habitat cover') + 
  ylab("Moran's I") + 
  theme_bw() + 
  theme(legend.position = 'top')
p1

p2 <- ggplot(df2, aes(x=round(cover,1), y=round(moransI, 1), fill=S0/1000)) + 
  geom_raster() + 
  scale_fill_viridis_c(name='S0 (x 1000)') + 
  xlab('% Habitat cover') + 
  ylab("Moran's I") + 
  theme_bw() + 
  theme(legend.position = 'top')

ggarrange(p1, p2)

mod1 <- lm(N0~cover*moransI*I(moransI^2)*I(cover^2),data=df2)
mod1 <- lm(formula = N0 ~ cover + moransI + I(moransI^2) + I(cover^2) + 
             cover:moransI + cover:I(moransI^2) + moransI:I(moransI^2) + 
             cover:I(cover^2) + moransI:I(cover^2) + I(moransI^2):I(cover^2) + 
             cover:moransI:I(moransI^2) + cover:moransI:I(cover^2) + cover:I(moransI^2):I(cover^2) + 
             moransI:I(moransI^2):I(cover^2) + 
             cover:moransI:I(moransI^2):I(cover^2), data = df2)
summary(mod1)
# R2 = 0.964

mod1b <- lm(formula = log(N0) ~ cover + moransI + I(cover^2) + 
              cover:moransI + moransI:I(moransI^2) + 
              cover:I(cover^2) + moransI:I(cover^2) + I(moransI^2):I(cover^2) + 
              cover:moransI:I(moransI^2) + cover:moransI:I(cover^2) + cover:I(moransI^2):I(cover^2) + 
              moransI:I(moransI^2):I(cover^2) + 
              cover:moransI:I(moransI^2):I(cover^2), data = df2)
summary(mod1b)
# R2 = 0.995

df2$fitted_N0 <- exp(mod1b$fitted.values)
p1 <- ggplot(df2, aes(x=N0, y=fitted_N0, color=cover)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  scale_color_viridis_c(name='% Habitat cover') + 
  geom_abline(slope=1, intercept=0) + 
  theme_bw() + 
  theme(legend.position = 'top')

mod2 <- lm(S0~cover*moransI*I(moransI^2)*I(cover^2),data=df2)
summary(mod2)
step(mod2b)

mod2b <- lm(formula = log(S0) ~ cover + I(cover^2) + 
              cover:moransI + cover:I(moransI^2) + moransI:I(moransI^2) + 
              cover:I(cover^2) + moransI:I(cover^2) + 
              cover:moransI:I(moransI^2) + cover:moransI:I(cover^2) + cover:I(moransI^2):I(cover^2) + 
              cover:moransI:I(moransI^2):I(cover^2), data = df2)
summary(mod2b)
# R2 = 0.994

df2$fitted_S0 <- exp(mod2b$fitted.values)
p2 <- ggplot(df2, aes(x=S0, y=fitted_S0, color=cover)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  scale_color_viridis_c(name='% Habitat cover') + 
  geom_abline(slope=1, intercept=0) + 
  theme_bw() + 
  theme(legend.position = 'top')
p2

ggarrange(p1, p2, common.legend = TRUE)

# N0 and S0 in the continuous lanscape were 2025000 and 5500, 
# what fraction of this is used in the best METE fit, and how
# do Moran's I and % cover affect these deviations?

p1 <- ggplot(df2, aes(x=round(cover,1), y=round(moransI, 1), fill=N0/2025)) + 
  geom_tile() + 
  scale_fill_viridis_c(name='best fitting N0 / original N0') + 
  xlab('% Habitat cover') + 
  ylab("Moran's I") + 
  theme_bw() + 
  theme(legend.position = 'top')
p1

p2 <- ggplot(df2, aes(x=round(cover,1), y=round(moransI, 1), fill=S0/5500)) + 
  geom_raster() + 
  scale_fill_viridis_c(name='best fitting S0 / original S0') + 
  xlab('% Habitat cover') + 
  ylab("Moran's I") + 
  theme_bw() + 
  theme(legend.position = 'top')
p2

ggarrange(p1, p2)

df2$dN0 <- df2$N0/2025
mod1 <- lm(formula = log(dN0) ~ cover + moransI + I(cover^2) + 
             cover:moransI + moransI:I(moransI^2) + 
             cover:I(cover^2) + moransI:I(cover^2) + I(moransI^2):I(cover^2) + 
             cover:moransI:I(moransI^2) + cover:moransI:I(cover^2) + cover:I(moransI^2):I(cover^2) + 
             moransI:I(moransI^2):I(cover^2) + 
             cover:moransI:I(moransI^2):I(cover^2), data = df2)
summary(mod1)
# R2 = 0.995

df2$fitted_dN0 <- exp(mod1$fitted.values)
p1 <- ggplot(df2, aes(x=dN0, y=fitted_dN0, color=cover)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  xlab('observed best fitting N0 / original N0') + 
  ylab('estimated best fitting N0 / original N0') + 
  scale_color_viridis_c(name='% Habitat cover') + 
  geom_abline(slope=1, intercept=0) + 
  theme_bw() + 
  theme(legend.position = 'top')
p1

df2$dS0 <- df2$S0/5500
mod2 <- lm(formula = log(dS0) ~ cover + I(cover^2) + 
             cover:moransI + cover:I(moransI^2) + moransI:I(moransI^2) + 
             cover:I(cover^2) + moransI:I(cover^2) + 
             cover:moransI:I(moransI^2) + cover:moransI:I(cover^2) + cover:I(moransI^2):I(cover^2) + 
             cover:moransI:I(moransI^2):I(cover^2), data = df2)
summary(mod2)
# R2 = 0.994

df2$fitted_dS0 <- exp(mod2$fitted.values)
p2 <- ggplot(df2, aes(x=dS0, y=fitted_dS0, color=cover)) + 
  geom_point() + 
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans='log10') + 
  xlab('observed best fitting S0 / original S0') + 
  ylab('estimated best fitting S0 / original S0') + 
  scale_color_viridis_c(name='% Habitat cover') + 
  geom_abline(slope=1, intercept=0) + 
  theme_bw() + 
  theme(legend.position = 'top')
p2

ggarrange(p1, p2, common.legend = TRUE)


ggplot(df2, aes(x=clustering, y=moransI)) + 
  geom_point()

ggplot(df2, aes(x=f_loss, y=cover)) + 
  geom_point()

