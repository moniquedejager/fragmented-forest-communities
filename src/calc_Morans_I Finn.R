# calculate Moran's I
  
# libraries:
library(GGally)    
 
# dit deel is uit mijn eigen script, dus daar moet de omgevingsdata van jou in 
# komen te staan. n is het aantal cellen (ervan uitgaande dat we met een vierkant 
# werken), x is een vector met alle x-coordinaten, y een vector met alle y-
# coordinaten, en h een vector met nullen en enen (nul = geen habitat cell, een
# = wel habitat cell). 
n <- rv$nx
x <- rep(1:n, n)
y <- sort(x)
h <- rv$comm_type == 'sub'
  
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

# Moran's I = 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2))  
I <- 1 / 4 * (sum(d0*(d1 + d2 + d3 + d4)) / sum(d0^2)) 
