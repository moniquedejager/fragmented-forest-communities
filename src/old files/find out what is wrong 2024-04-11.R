

patchID  <- sort(rep(1:10, 20))
spec     <- 1:200

m_spec1 <- m_spec2 <- vector(length = 0)

for (i in 1:100){

newpatch <- sample(1:10, 200, replace=TRUE)
newind   <- sample(1:20, 200, replace=TRUE)
sel      <- (patchID > 1)

a        <- sample(1:sum(sel), sum(sel), replace=FALSE)

spec_new <- spec
spec_new[(newpatch[sel][a] - 1)*20 + newind[sel][a]] <- spec[sel][a]
m_spec1 <- c(m_spec1, mean(spec_new))

spec_new <- spec
ind_id <- ((newpatch - 1)*20 + newind)[sel]
spec_new[ind_id[a]] <- spec[sel][a]
m_spec2 <- c(m_spec2, mean(spec_new))

}

boxplot(m_spec1, m_spec2)
abline(h = 100.5)

