# create the input file for the simulations in C++:
params <- expand.grid(simNumber = 1:5, 
                      f_loss = 0:9/10,
                      clustering = seq(1, 5, 0.5))
params$initialization = 0
params$initialization[params$f_loss == 0] = 1
params$n_species = 1000

params <- params[order(params$f_loss),]

end_file <- data.frame(simNumber = -1, 
                       f_loss = -1, 
                       clustering = -1,
                       initialization = 0,
                       n_species = -1)
params <- rbind(params, end_file)

write.table(params, 'results/parameter_values.txt', 
            row.names = FALSE, 
            col.names = FALSE,
            append = FALSE)


# input file for calibration of n_species:
params <- data.frame(simNumber = 1,
                     f_loss = 0,
                     clustering = 1,
                     initialization = 1,
                     n_species = (3:10)*5000)

end_file <- data.frame(simNumber = -1, 
                       f_loss = -1, 
                       clustering = -1,
                       initialization = 0,
                       n_species = -1)
params <- rbind(params, end_file)

write.table(params, 'results/parameter_values.txt', 
            row.names = FALSE, 
            col.names = FALSE,
            append = FALSE)
