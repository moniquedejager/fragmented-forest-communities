library(ggplot2)
setwd('C:/Users/jager015/OneDrive - Universiteit Utrecht/Documents/Edwin/Fragmentation, restoration, dispersal and communities')
source('./src 03-2024/summarySE.R')

df           <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c('A', 'S', 'type', 'mu', 'sim_nr')

for (mu in 1:4){
  for (sim_nr in 1:10){
    filename <- paste('./results/simulation results 03-2024/community_composition/', mu, '_', sim_nr,'_95.txt', sep='')
    if (file.exists(filename)){
      filename <- paste('./results/simulation results 03-2024/community_composition/', mu, '_', sim_nr,'_0.txt', sep='')
      m <- read.table(filename)
      m <- as.matrix(m)
      
      # first, create the species-area relationship with the contiguous envirnonment:
      A <- vector(length = 0)
      S <- vector(length = 0)
      
      for (i in seq(2500, 125, -125)){
        for (j in 1:100){
          ids <- sample(1:2500, i, replace = FALSE)
          A <- c(A, i)
          S <- c(S, length(unique(as.vector(m[,ids]))))
        }
      }
      
      df2 <- data.frame(A = A,
                       S = S,
                       type = 'backwards SAR',
                       mu = mu,
                       sim_nr = sim_nr)
      
      df <- rbind(df, df2)
      
      # What happens when habitat loss occurs? First, let's see what happens
      # in terms of geometric/static habitat loss:
      filename  <- paste('./results/simulation results 03-2024/landscapes/landscapes_', sim_nr, '_', mu, '.txt', sep='')
      landscape <- read.table(filename)
      
      A <- colSums(landscape)
      S <- sapply(1:20, function(ix){
        S <- length(unique(as.vector(m[,landscape[,ix] == 1])))
        return(S)
      })
      
      df2 <- data.frame(A = A,
                        S = S,
                        type = 'Geometric/Static',
                        mu   = mu, 
                        sim_nr = sim_nr)
      
      df <- rbind(df, df2)
      
      # Now, let's find out the extinction debt in case of the demographic/dynamic
      # part:
      A <- A
      S <- sapply(seq(0, 95, 5), function(ix){
        filename <- paste('./results/simulation results 03-2024/community_composition/', mu, '_', sim_nr,'_', ix,'.txt', sep='')
        if (file.exists(filename)){
          m <- read.table(filename)
          m <- as.matrix(m)
          S <- (length(unique(as.vector(m))))
        } else { S <- NA }
        return(S)
      })
      
      df2 <- data.frame(A = A,
                        S = S,
                        type = 'Demographic/Dynamic',
                        mu   = mu, 
                        sim_nr = sim_nr)
      
      df <- rbind(df, df2)
    }
  }
}

sel <- (df$type == 'backwards SAR')
ggplot(df[sel,], aes(x=A, y=S, color = as.factor(mu))) + 
  geom_jitter(width=50, alpha=0.3)

sdf <- summarySE(df, measurevar="S", groupvars=c("A","mu", "type"))
ggplot(sdf, aes(x=A, y=S, colour=as.factor(mu), shape = as.factor(mu))) + 
  geom_errorbar(aes(ymin=S-ci, ymax=S+ci), width=0) +
  geom_line() +
  geom_point() + 
  facet_wrap(vars(type)) + 
  scale_x_continuous(trans='log10', limits=c(100, 2500)) + 
  scale_y_continuous(trans='log10') + 
  xlab('Area') + 
  ylab('# Species')

# Calculate extinction debt: --> moet anders nu we heel veel samples van de backwards SAR hebben!
df2 <- data.frame(A = df$A[df$type == 'backwards SAR'],
                  E = df$S[df$type == 'Demographic/Dynamic'] - 
                    df$S[df$type == 'backwards SAR'],
                  type = 'Demographic/Dynamic',
                  mu   = df$mu[df$type == 'backwards SAR'], 
                  sim_nr = df$sim_nr[df$type == 'backwards SAR'])
df3 <- data.frame(A = df$A[df$type == 'backwards SAR'],
                  E = df$S[df$type == 'Geometric/Static'] - 
                    df$S[df$type == 'backwards SAR'],
                  type = 'Geometric/Static',
                  mu   = df$mu[df$type == 'backwards SAR'], 
                  sim_nr = df$sim_nr[df$type == 'backwards SAR'])
df2 <- rbind(df2, df3)

sdf <- summarySE(df2, measurevar="E", groupvars=c("A","mu", "type"))
pd  <- position_dodge(0)
ggplot(sdf, aes(x=A, y=E, colour=as.factor(mu), shape = as.factor(mu))) +
  geom_hline(yintercept = 0, linetype='dashed', lwd=1) + 
  geom_errorbar(aes(ymin=E-se, ymax=E+se), width=0, position=pd) +
  geom_line(linewidth=1, position=pd) +
  geom_point(size=3, position=pd) + 
  scale_x_continuous(trans='log10') + 
  facet_wrap(vars(type)) + 
  xlab('Area') + 
  ylab('Extinction debt')

# difference in extinction debt between geometric and demographic:
df3 <- data.frame(A = df2$A[df2$type == 'Geometric/Static'],
                  dE = df2$E[df2$type == 'Geometric/Static'] - 
                    df2$E[df2$type == 'Demographic/Dynamic'],
                  mu   = df2$mu[df2$type == 'Geometric/Static'], 
                  sim_nr = df2$sim_nr[df2$type == 'Geometric/Static'])
df3$hab_loss <- 100 - df3$A/25
sdf <- summarySE(df3, measurevar="dE", groupvars=c("hab_loss","mu"))

ggplot(sdf, aes(x=hab_loss, y=dE, colour=as.factor(mu), shape = as.factor(mu))) +
  geom_hline(yintercept = 0, linetype='dashed', lwd=1) + 
  geom_errorbar(aes(ymin=dE-se, ymax=dE+se), width=0) +
  geom_line(linewidth=1) +
  geom_point(size=3) + 
  xlab('% Habitat loss') + 
  ylab('Difference in extinction debt') + 
  annotate(geom="text", x=50, y=30, 
           label="Geometric extinction debt is higher") + 
  annotate(geom="text", x=50, y=-25, 
           label="Demographic extinction debt is higher") + 
  theme_classic()
