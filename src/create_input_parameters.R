# create the input file for the simulations in C++:
params <- expand.grid(simNumber = 1:5, 
                      f_loss = 0:9/10,
                      clustering = seq(1, 5, 0.5))
params$initialization = 0
params$initialization[params$f_loss == 0] = 1
params$n_species = 5500
params$mutation_rate = 0.002
params$H = 6.5

params <- params[order(params$f_loss),]

end_file <- data.frame(simNumber = -1, 
                       f_loss = -1, 
                       clustering = -1,
                       initialization = 0,
                       n_species = -1,
                       mutation_rate = -1,
                       H = -1)
params <- rbind(params, end_file)

write.table(params, 'results/parameter_values.txt', 
            row.names = FALSE, 
            col.names = FALSE,
            append = FALSE)


# input file for calibration of mutation rate (= probability of receiving a
# random offspring from the larger metacommunity) & H:

x1 <- expand.grid(H = seq(4.5, 6.5, 0.5),
                 mutation_rate = seq(0.001, 0.002, 0.0002),
                 n_species = seq(5500, 10000, 500))

x2 <- expand.grid(H = seq(4.5, 7.5, 0.5),
                  mutation_rate = seq(0.001, 0.0024, 0.0002),
                  n_species = seq(4500, 10000, 500))

group1 <- paste(x1$H, x1$mutation_rate, x1$n_species, sep='-')
group2 <- paste(x2$H, x2$mutation_rate, x2$n_species, sep='-')
x <- x2[!(group2 %in% group1),]

x <- x[!(x$H %in% c(5, 6)),]
x <- x[x$mutation_rate > 0.0017,]
x <- x[x$n_species < 6500,]

params <- data.frame(simNumber = 1,
                     f_loss = 0,
                     clustering = 1,
                     initialization = 1,
                     n_species = x$n_species,
                     mutation_rate = x$mutation_rate,
                     H = x$H)
                     
                     #mutation_rate = c(0.00001, 0.00003, 0.0001, 0.0003, 
                      #                 0.001, 0.003, 0.01, 0.03))

end_file <- data.frame(simNumber = -1, 
                       f_loss = -1, 
                       clustering = -1,
                       initialization = 0,
                       n_species = -1,
                       mutation_rate = -1,
                       H = -1)
params <- rbind(params, end_file)

write.table(params, 'results/parameter_values.txt', 
            row.names = FALSE, 
            col.names = FALSE,
            append = FALSE)

params$mutation_rate*1000

### for Fleur:

params <- data.frame(simNumber = 1:5,
                     f_loss = 0,
                     clustering = 1,
                     initialization = 1,
                     n_species = 5500,
                     mutation_rate = 0.002,
                     H = 6.5)

#mutation_rate = c(0.00001, 0.00003, 0.0001, 0.0003, 
#                 0.001, 0.003, 0.01, 0.03))

end_file <- data.frame(simNumber = -1, 
                       f_loss = -1, 
                       clustering = -1,
                       initialization = 0,
                       n_species = -1,
                       mutation_rate = -1,
                       H = -1)
params <- rbind(params, end_file)

write.table(params, 'results/parameter_values.txt', 
            row.names = FALSE, 
            col.names = FALSE,
            append = FALSE)
