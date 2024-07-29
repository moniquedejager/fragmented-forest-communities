# Assuming that the relation between the inverse cumulative 
# probability of migration and distance follows an exponential distribution, 
# we were able to extraplolate the probability that an individual comes from 
# the metapopulation from real datasets to our model case. We hereby use the 
# estimates of the probabilities of local, adjacent, and metacommunity 
# migrations given in table 2 of Pos et al. 2017. 

# data from table 2 of Pos et al. 2017:
m.local <- c(0.79, 0.79, 0.75, 0.75, 0.69, 0.69, 0.65, 0.59, 0.59, 0.55, 0.49,
             0.49, 0.45, 0.39, 0.29, 0.19)
m.adj <- c(0.2, 0.01, 0.05, 0.2, 0.3, 0.01, 0.3, 0.4, 0.01, 0.4, 0.5, 
           0.01, 0.5, 0.01, 0.01, 0.01)
m.meta <- c(0.01, 0.2, 0.2, 0.05, 0.01, 0.3, 0.05, 0.01, 0.4, 0.05, 0.01, 
            0.5, 0.05, 0.6, 0.7, 0.8)

# Per study, we estimate lambda and thereafter the probability that an
# individual comes from the metapopulation in our simulations (i.e. not 
# from any of the 400 subcommunities holding 1000 individuals each):

Pm <- vector(length=0)
for (i in 1:16){
  m1 <- m.adj[i] + m.meta[i]
  m2 <- m.meta[i]
  x1 <- rep(1, length(m1))
  x2 <- rep(2, length(m1))
  x0 <- rep(0, length(m1))
  m0 <- x0 + 1
  x <- c(x0, x1, x2)
  y <- c(m0, m1, m2)
  plot(x,y, log='y')
  ly <- log(y)
  plot(x,ly)
  smod <- summary(lm(ly~0+x))
  
  l <- 16
  Pm <- c(Pm, exp(smod$coefficients[1]*l))
}

boxplot(Pm, log='y')
abline(h=0.0003)

range(Pm)
median(Pm)
