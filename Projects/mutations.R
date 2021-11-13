source("../simu.R")

# Mutations introduce two new parameters:
# rate.mut is the mutation rate per haplotype
# var.mut  is the variance of a mutational effect
#
# (warning: var.mut is NOT the mutational variance, 
#  varM = 2*rate.mut*var.mut) 

############## Here is the only significant change: 
# If a mutation occurs during gametogenesis (with rate rate.mut), 
# a random locus is modified by a Gaussian random numbber of variance var.mut

make.gamete <- function(indiv, rate.mut, var.mut) {
	gam <- indiv$genotype[cbind(1:nrow(indiv$genotype), sample(c(1,2), nrow(indiv$genotype), replace=TRUE))]
	if (runif(1) < rate.mut) {
		mut.loc <- sample(seq_along(gam), 1)
		gam[mut.loc] <- rnorm(1, mean=gam[mut.loc], sd=sqrt(var.mut))
	}
	gam
}

############## The rest is just copy/paste with the two new parameters transmitted down from simulation() to make.gamete()

make.offspring <- function(mother, father, var.env, rate.mut, var.mut) {
	# Makes an individual out of two parents. 
	genotype <- cbind(make.gamete(mother, rate.mut, var.mut), make.gamete(father, rate.mut, var.mut))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}

reproduction <- function(population, pop.size, var.env, rate.mut, var.mut) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(n=pop.size, 
		expr=make.offspring(
				mother=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE), 
				father=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE),
				var.env=var.env ,
				rate.mut=rate.mut, var.mut=var.mut),
				simplify = FALSE)
}

simulation <- function(generations=20, pop.size = 100, num.loci = 5, var.init = 1, var.env = 1, sel.strength = 1, sel.optimum = 0, rate.mut=0.0001, var.mut=1) {
	# Runs a simulation
	pop <- init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env)
	summ <- data.frame()
	for (gg in 1:generations) {
		pop <- update.fitness(pop, sel.strength, sel.optimum)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(pop, pop.size=pop.size, var.env=var.env, rate.mut=rate.mut, var.mut=var.mut)
	}
	summ
}

################################################################################################################################
##########################################
# Exploration of the effect of mutations #
##########################################

# Effect on long-term directional selection

ngen <- 20 # easier to see the effect with long simulations
mutrates <- c(0, 0.02, 0.1)
N        <- 1000
var.mut  <- 1
n.loci   <- 10

pdf("mutations-fig.pdf", width=10, height=15)
	layout(1:3)
	par(cex=1)
	
	plot(NULL, xlim=c(0, ngen), ylim=c(0, 10), xlab="Generations", ylab="Mean phenotype")
	for (i in seq_along(mutrates)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=10, rate.mut=mutrates[i], var.mut=var.mut)$phen.mean, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(mutrates)), legend=rev(paste0("mu = ", mutrates)))
	mtext("Mutations help the population to achieve an optimum that was \nnot reachable from the standing genetic variation", 3)
	
	plot(NULL, xlim=c(0, ngen), ylim=c(0.4, 0.7), xlab="Generations", ylab="Mean fitness")
	for (i in seq_along(mutrates)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=0, rate.mut=mutrates[i], var.mut=var.mut)$fit.mean, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(mutrates)), legend=rev(paste0("mu = ", mutrates)))
	mtext("Mutations generate genetic load when the population is at the fitness optimum (=0)", 3)

	plot(NULL, xlim=c(0, ngen), ylim=c(0, 1.5), xlab="Generations", ylab="Genetic variance")
	for (i in seq_along(mutrates)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=0, rate.mut=mutrates[i], var.mut=var.mut)$gen.var, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(mutrates)), legend=rev(paste0("mu = ", mutrates)))
	mtext("Mutations contribute to the equilibrium genetic variance (optimum=0)", 3)

dev.off()
