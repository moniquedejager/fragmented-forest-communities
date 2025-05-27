# dinsdag 19 mei: code weer aanpassen, want nu klopt het weer niet meer. 
# de +1-en uit de KL berekeningen halen;
# ipv de opgeslagen SAD van de originele continue habitat gebruiken, 
# alleen SAD berekenen over de patches die ook daadwerkelijk in de nieuwe
# simulatie zitten (die met loss en fragmentatie).

# update: na een nachtje slapen denk ik dat het beter is om bij een habitat
# patch te beginnen, en dan gewoon dezelfde methode aan te houden als hiervoor. 
# hoeven we de +1 ook niet te gebruiken. En aangezien het continuous boundaries 
# zijn, hoeven we daar ook niet over in te zitten.



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

plot(sar1$pred$A, sar1$pred$S)

# now, we can use the data from the pristine environment to calculate the 
# biodiversity loss (in terms of numbers of species) per area size, for all
# simulations:

areas <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
n     <- 45
x     <- rep(1:n, n)
y     <- sort(x)
patchID <- paste(y-1, x-1, sep='-')

df2 <- data.frame(simnr      = 0, # simulation number
                  clustering = 0, # clustering parameter (1 is highly fragmented, 5 is clustered)
                  f_loss     = 0, # fraction of habitat lost
                  A          = 0, # Area size (in patches)
                  S          = 0, # No. species in undestructed forest
                  a          = 0, # No. habitat patches in area of size S
                  s          = 0, # No. species in case of habitat loss and fragmentation
                  KL_A       = 0, # Adapted Kullback-Leibner divergence of A
                  KL_S       = 0, # Adapted Kullback-Leibner divergence of S
                  MoransI    = 0, # Moran's I
                  h          = 0) # Fraction of area that consists of habitat

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
      
      df$patch <- paste(df$x, df$y, sep='-')
      uPatch   <- unique(df$patch)
      
      # reorder to start at the first habitat patch:
      i <- (1:2025)[patchID == uPatch[1]]
      patchID2 <- patchID[c(i:2025, 1:(i-1))]
      
      a <- length(unique(df$patch[df$patch %in% patchID2[1]]))
      s <- length(unique(df$species[df$patch %in% patchID2[1]]))
      A <- 1
      S <- sar1$pred$S[1]
      KL_A1 <- (a - A)^2  #log(a/A)
      KL_S1 <- (s - S)^2  #log(s/S)
      KL_A  <- KL_A1
      KL_S  <- KL_S1
      
      for (i in areas){
        a1 <- length(unique(df$patch[df$patch %in% patchID2[1:i]]))
        a  <- c(a, a1)
        s1 <- length(unique(df$species[df$patch %in% patchID2[1:i]]))
        s  <- c(s, s1)
        A1 <- i  #(1:2025)[patchID == patchID[i]]
        A  <- c(A, A1)
        S1 <- sar1$pred$S[sar1$pred$A == i]
        S  <- c(S, S1)
        
        KL_A1 <- KL_A1 + (a1 - A1)^2 #log(a1/A1)
        KL_A  <- c(KL_A, KL_A1)
        KL_S1 <- KL_S1 + (s1 - S1)^2 #log(s1/S1)
        KL_S  <- c(KL_S, KL_S1)
      }
      
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

      df2a <- data.frame(simnr      = sim_nr,
                         clustering = clustering,
                         f_loss     = f_loss,
                         A          = A,
                         S          = S,
                         a          = a,
                         s          = s,
                         KL_A       = KL_A,
                         KL_S       = KL_S,
                         MoransI    = I,
                         h          = mean(h))
      
      df2 <- rbind(df2, df2a)
    }
  }
}
df2 <- df2[df2$simnr > 0,]

# logit transform MoransI: 
df2$MoransI_logit <- log(abs(df2$MoransI) / (1 - abs(df2$MoransI)))
df2$h_logit <- log(df2$h/(1 - df2$h))

ggplot(df2, aes(x=a, y=s, color=A)) + 
  geom_point() +
  facet_grid(cols=vars(round(h, 1)), rows=vars(round(MoransI, 1))) + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10') + 
  scale_color_viridis_c(trans='log10')

ggplot(df2, aes(x=A, y=KL_A)) + 
  geom_point() + 
  facet_grid(cols=vars(round(h, 1)), rows=vars(round(MoransI, 1))) 

sel <- (df2$simnr == 1)&(df2$clustering %in% c(1,5))&(df2$f_loss %in% c(0.1, 0.9))
ggplot(df2[sel,], aes(x=A, y=KL_A, color=factor(f_loss))) + 
  geom_line(aes(linetype=factor(clustering))) + 
  geom_point(aes(shape=factor(clustering)))
ggplot(df2[sel,], aes(x=A, y=KL_S, color=factor(f_loss))) + 
  geom_line(aes(linetype=factor(clustering))) + 
  geom_point(aes(shape=factor(clustering)))


ggplot(df2, aes(x=A, y=KL_S)) + 
  geom_line(aes(linetype=factor(simnr))) +
  geom_point(aes(color=a+1)) + 
  scale_color_viridis_c(trans='log10') + 
  facet_grid(rows = vars(clustering), cols=vars(f_loss))


ggplot(df2, aes(x=sqrt(KL_A), y=sqrt(KL_S))) + 
  geom_point() + 
  facet_grid(cols=vars(round(h_logit, 1)), rows=vars(round(MoransI_logit))) + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans='log10') + 
  scale_x_continuous(trans='log10')


sel <- (df2$KL_A > 0)&(df2$KL_S > 0)
mod <- lm(log(sqrt(KL_S) + 1) ~ 0 + log(sqrt(KL_A) + 1) * MoransI_logit * h_logit * I((log(sqrt(KL_A) + 1)^2)), data=df2[sel,])
summary(mod)

df2$fitted <- NA
df2$fitted[sel] <- mod$fitted.values

ggplot(df2, aes(x=log(sqrt(KL_S) + 1), y=fitted)) + geom_point() + 
  geom_abline(slope=1, intercept=0)

sel <- (df2$clustering %in% c(1, 3, 5)) & (df2$f_loss %in% c(0.1, 0.5, 0.9))

p1 <- ggplot(df2[sel,], aes(x=A, y=KL_A, color=factor(simnr))) + 
  geom_line() + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss), scales='free_y') + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()

p2 <- ggplot(df2[sel,], aes(x=A, y=KL_S, color=factor(simnr))) + 
  geom_line() + 
  geom_line(aes(x=A, y=fitted, color=factor(simnr)), linetype=2) + 
  facet_grid(cols=vars(clustering), rows=vars(f_loss), scales='free_y') + 
  #scale_y_continuous(trans='log10') + 
  theme_bw()

ggarrange(p1, p2, nrow=1, ncol=2, common.legend = TRUE)



