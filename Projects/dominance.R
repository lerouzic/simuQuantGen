source("../simu.R")

# Dominance requires an additional parameter: the dominance coefficient.
# dom.coef:
#           < -1: underdominance
#           = -1: the "lower" allele is completely dominant
#           < 0 : partial dominance of the lower allele
#           = 0 : no dominance
#           > 0 : partial dominance of the larger allele
#           = 1 : total dominance
#           > 1 : overdominance 

############## The significant change happens here: 

GPmap <- function(genotype, dom.coef) {
	# Returns the genotypic value (mean phenotype) corresponding to a genotype
	diffgen <- apply(genotype, 1, diff)
	sum(2*rowMeans(genotype) + dom.coef*abs(diffgen))
}

############## The rest is just copy/paste with the two new parameters transmitted down from simulation() to GPmap()

get.phenotype <- function(genotype, var.env, dom.coef) {
	rnorm(1, mean=GPmap(genotype, dom.coef), sd=sqrt(var.env))
}

init.individual <- function(var.init, num.loci, var.env, dom.coef) {
	genotype <- matrix(
		rnorm(2*num.loci, mean=0, sd=sqrt(var.init/2/num.loci)), 
		ncol=2)
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype, dom.coef),
		phenotype = get.phenotype(genotype, var.env, dom.coef),
		fitness   = 1
	)
}

init.population <- function(pop.size, var.init, num.loci, var.env, dom.coef) {
	replicate(pop.size, init.individual(var.init, num.loci, var.env, dom.coef), simplify=FALSE)
}

make.offspring <- function(mother, father, var.env, dom.coef) {
	# Makes an individual out of two parents. 
	genotype <- cbind(make.gamete(mother), make.gamete(father))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype, dom.coef),
		phenotype = get.phenotype(genotype, var.env, dom.coef),
		fitness   = 1
	)
}

reproduction <- function(population, pop.size, var.env, dom.coef) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(n=pop.size, 
		expr=make.offspring(
				mother=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE), 
				father=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE),
				var.env=var.env, 
				dom.coef=dom.coef),
				simplify = FALSE)
}

simulation <- function(generations=20, pop.size = 100, num.loci = 5, var.init = 1, var.env = 1, sel.strength = 1, sel.optimum = 0, dom.coef = 0) {
	# Runs a simulation
	pop <- init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env, dom.coef=dom.coef)
	summ <- data.frame()
	for (gg in 1:generations) {
		pop <- update.fitness(pop, sel.strength, sel.optimum)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(pop, pop.size=pop.size, var.env=var.env, dom.coef=dom.coef)
	}
	summ
}


################################################################################################################################
###################################################
# Exploration of the effect of dominance          #
###################################################

ngen <- 20 # easier to see the effect with long simulations
dom <- c(0, 0.5, 1, 1.5)
N        <- 1000
n.loci   <- 5
sel.strength <- 2

pdf("dominance-fig.pdf", width=10, height=15)
	layout(1:3)
	par(cex=1)
	
	plot(NULL, xlim=c(0, ngen), ylim=c(-10, 10), xlab="Generations", ylab="Mean phenotype")
	for (i in seq_along(dom)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=10, dom.coef=dom[i], sel.strength=sel.strength)$phen.mean, col=i)
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=-10, dom.coef=dom[i], sel.strength=sel.strength)$phen.mean, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(dom)), legend=rev(paste0("h = ", dom)))
	mtext("Dominance makes selection response asymmetric", 3)


	plot(NULL, xlim=c(0, ngen), ylim=c(0, 3), xlab="Generations", ylab="Genetic variance")
	for (i in seq_along(dom)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=10, dom.coef=dom[i], sel.strength=sel.strength)$gen.var, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(dom)), legend=rev(paste0("h = ", dom)))
	mtext("Overdominance prevents the fixation of beneficial alleles (optimum=10)", 3)


	plot(NULL, xlim=c(0, ngen), ylim=c(0, 1), xlab="Generations", ylab="Mean fitness")
	for (i in seq_along(dom)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=5, dom.coef=dom[i], sel.strength=sel.strength)$fit.mean, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(dom)), legend=rev(paste0("h = ", dom)))
	mtext("Adaptation through overdominance has a cost (optimum=5)", 3)
dev.off()


