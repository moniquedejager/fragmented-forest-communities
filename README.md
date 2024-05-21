# fragmented-forest-communities: restoration

## Methods
In summary, we simulated subcommunity dynamics in a 2-dimensional, semi-spatial, near-neutral, individual-based model. 
The environment consisted of 20 x 20 cells; in each cell, a subcommunity consisting of 1,000 individuals could exist. 
This restoration part of the model is preceded by a fragmentation part. 
Initially, all cells were habitable and contained a subcommunity in this first part of the model. 
These pristine landscapes underwent increasing habitat loss in steps of 5%.
Each time, community dynamics needed to be stabilized before habitat destruction continued. 
Habitat destruction was simulated at three fragmentation levels: (i) clustered, (ii) fractal, and (iii) random habitat destruction. 
Each simulation was run 10 times to account for the stochastic nature of the model.
 
Each fragmented landscape resulting from part 1 was used as the basis for the restoration simulations.
For each of these landscapes, 5, 10, 15, etc. percent of the total area was restored, until all patches contained habitat.
Habitat patches were added using three different spatial configurations: (i) clustered, (ii) fractal, and (iii) random habitat restoration. 
Each simulation was run until the average number of species showed no trend in 50 iterations for 5 times.  

All model simulations and analyses of the resulting data were performed in R version xx (ref to Zenodo). 

### Creating initial communities 
First, the subcommunities were initialized in the modelled landscape. We created a 2-dimensional lattice of 20 x 20 cells with continuous boundaries, where each cell contained a subcommunity consisting of 1,000 individuals. Each subcommunity started with one individual per species of species S1 to S1,000. Each timestep, every individual was able to disperse within its own or to a neighboring subcommunity and replace another individual. The Euclidean distance between cells determined the probability that an individual from subcommunity j could replace an individual in subcommunity i (Fig. S1). To limit matrix sizes and thereby increase computation speed, dispersal was limited to the closest 121 cells (within the ranges -5 ≤ Δxij ≤ 5 and -5 ≤ Δyij ≤ 5, where Δxij and Δyij are the distances (in numbers of patches) in x and y directions between the focal subcommunity i and its neighbor j). For each of these 121 cells surrounding (and including) a focal subcommunity, the probability to disperse here was calculated using a 2-dimensional exponential probability density function: 

 P_(i,j)=(2πλ^2 )^(-1)∙e^(-l⁄λ),							(eq. 1)

where λ is the scaling exponent and l is the Euclidian distance between subcommunities i and j. Other probability density functions can also be used, to describe different dispersal kernels. We ran two different simulation blocks, one in which all species dispersed with λ = 0.3, hereafter called same-dispersal simulations, and one in which species differed in their dispersal capacity (0 < λ < 0.5), hereafter called different-dispersal simulations. 
Every individual in each subcommunity randomly selected a subcommunity to disperse to, based on its dispersal kernel (weighted random selection on Pi,j). The individual to replace from that subcommunity was randomly selected as well (without any differences in weights). Because such random replacements may overlap (i.e. an individual was replaced twice in the same iteration) and R always uses the last value as the resulting replacement, we randomly ordered the replacements. Replacement of the individuals in the subcommunities was carried out iteratively, until there were 5 sequences of 50 timesteps in which no significant relation (p < 0.05) was to be found between timestep number and the average number of species per subcommunity. The resulting subcommunities were recorded for further use in subsequent simulations.

### Fragmenting the environment
We simulated increasing habitat loss in steps of 5% habitat destruction (resulting in 20 sequential environments with different fragmentation levels), and with 3 different levels of clustering of habitat destruction. The first non-habitat cell was randomly selected from all cells; all other cells remained habitat. While the number of non-habitat cells created was smaller than the desired number of non-habitat cells (which depends on the fraction of the environment that should contain habitat), non-habitat cells are sequentially selected by means of weighted random selection, where weights of cells were attributable to their distance from destructed, non-habitat cells. Weights (Wi) were calculated from distances using a 2-D pareto distribution: 

 W_i=  ∑_(j=1)^n▒〖1/2π∙(2- μ)/(l_max^(2-μ)- l_min^(2-μ) )∙l_ij^(-μ) 〗,						(eq. 2)

where lij is the Euclidian distance between subcommunity i and non-habitat cell j, n is the number of non-habitat cells, lmin and lmax are the minimum and maximum distances to consider (l¬min = 1, lmax = 50), and µ is the scaling exponent (µ > 1). When µ is small, the non-habitat cells are distributed more homogeneously; non-habitat cells become more clustered with increasing µ. The different spatial configurations of habitat destruction were created using µ = 1 for random, µ = 3 for fractal, and µ = 5 for clustered habitat destruction (Fig. 1).
The spatial configuration of the environments with clustered, fractal, and random habitat destruction differ substantially from each other (Fig. S2). Under random habitat destruction, Moran’s I, an indicator of spatial autocorrelation, is close to zero, the average patch size decreases rapidly with habitat loss, and the number of patches shows a large peak at 75% habitat loss. In contrast, under clustered habitat destruction, spatial autocorrelation is high, mean patch size decreases less rapidly with habitat loss, and the number of patches shows remains much more constant across the different levels of habitat loss.  

### Community dynamics  
Simulations are initialized using the communities that were recorded at the end of the previous simulation round (e.g. the 5% fragmentation simulation starts with the communities as they were at the end of the initialization round, the 10% fragmentation simulation starts with the communities resulting from the 5% fragmentation simulation, etc.). During a simulation round, all individuals iteratively reproduce and replace other individuals, as described above in the section Creating initial communities. New species could form with a mutation rate of 0.0001. These new species could maximally differ Δλ = 0.05 in their dispersal strategy from their origin species in the different-species simulations, and Δλ = 0 in the same-species simulations.    

### Data analysis
For each simulation, we recorded – per subcommunity – (i) the number of species, (ii) the fraction of individuals originating from ancestors that inhabited the current patch 50 timesteps earlier, and (iii) the fit of the prediction of Maximum Entropy Theory in Ecology (METE) to the species abundance distribution (SAD). We estimated SAD using the sad function from the package meteR (Rominger & Merow, 2016) on the sampled data. For the given number of species present, we generated 1,000 SAD’s and calculated the average abundance per rank. Using this average estimated SAD and the actual simulated SAD, we calculated the goodness of fit (GOF¬METE):

 〖GOF〗_METE= 1- ∑_(i=1)^S▒〖(n_(o,i)-n_(e,i) )^2 〗/∑_(i=1)^S▒〖(n_(o,i)-(n_o ) ̅)^2 〗,			(eq. 3)

where S is the number of species present in the local community, n0,i and ne,i are, respectively, the number of individuals observed and expected (by the METE estimation) for rank i, and (n_o ) ̅ is the average number of individuals per species. Values of GOFMETE¬ close to 1 indicate that the prediction by METE of the SAD is quite well and, hence, the subcommunity appears to be in a good state. For each simulation, we furthermore recorded the entire metacommunity structure (x- and y-coordinates, and species number per individual), as these were needed in the habitat restoration simulations. We furthermore calculated Bray-Curtis dissimilarity, after each simulation round, between all possible combinations of the 20 subcommunities that were left after 95% habitat loss.  
