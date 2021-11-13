source("../simu.R")

# Recombinations require an additional parameter, the recombination rate.
# Complex vs simple genetic maps come almost "for free": if rate.rec is a vector, 
# then it encodes all recombination rates. 

############## The significant change happens here: 

make.gamete <- function(indiv, rate.rec) {
	recs <- cumsum(runif(length(rate.rec)+1) < c(0.5, rate.rec))
	indiv$genotype[cbind(1:nrow(indiv$genotype), 1+(recs %% 2))]
}

############## The rest is just copy/paste with the two new parameters transmitted down from simulation() to make.gamete()

make.offspring <- function(mother, father, var.env, rate.rec) {
	# Makes an individual out of two parents. 
	genotype <- cbind(make.gamete(mother, rate.rec), make.gamete(father, rate.rec))
	list(
		genotype  = genotype, 
		genot.value= GPmap(genotype),
		phenotype = get.phenotype(genotype, var.env),
		fitness   = 1
	)
}

reproduction <- function(population, pop.size, var.env, rate.rec) {
	# Returns the next generation
	fitnesses <- sapply(population, "[[", "fitness")
	replicate(n=pop.size, 
		expr=make.offspring(
				mother=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE), 
				father=unlist(sample(population, 1, prob=fitnesses), recursive=FALSE),
				var.env=var.env ,
				rate.rec=rate.rec),
				simplify = FALSE)
}

simulation <- function(generations=20, pop.size = 100, num.loci = 5, var.init = 1, var.env = 1, sel.strength = 1, sel.optimum = 0, rate.rec=0.5) {
	# Runs a simulation
	if (length(rate.rec) == 1)
		rate.rec <- rep(rate.rec, num.loci-1)
	stopifnot(length(rate.rec) == num.loci-1)
	
	pop <- init.population(pop.size=pop.size, var.init=var.init, num.loci=num.loci, var.env=var.env)
	summ <- data.frame()
	for (gg in 1:generations) {
		pop <- update.fitness(pop, sel.strength, sel.optimum)
		summ <- rbind(summ, summary.population(pop))
		if (gg < generations)
			pop <- reproduction(pop, pop.size=pop.size, var.env=var.env, rate.rec=rate.rec)
	}
	summ
}


################################################################################################################################
###################################################
# Exploration of the effect of recombination rate #
###################################################

ngen <- 20 # easier to see the effect with long simulations
recrates <- c(0, 0.1, 0.5)
N        <- 1000
n.loci   <- 10
sel.strength <- 2

pdf("recombination-fig.pdf", width=10, height=10)
	layout(1:2)
	par(cex=1)
	
	plot(NULL, xlim=c(0, ngen), ylim=c(0, 10), xlab="Generations", ylab="Mean phenotype")
	for (i in seq_along(recrates)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=10, rate.rec=recrates[i], sel.strength=sel.strength)$phen.mean, col=i)
	}
	legend("bottomright", lty=1, col=rev(seq_along(recrates)), legend=rev(paste0("r = ", recrates)))
	mtext("Recombinations helps adaptation", 3)


	plot(NULL, xlim=c(0, ngen), ylim=c(0, 2), xlab="Generations", ylab="Genetic variance")
	for (i in seq_along(recrates)) {
		lines(simulation(generations=ngen, pop.size=N, num.loci=n.loci, sel.optimum=0, rate.rec=recrates[i], sel.strength=sel.strength)$gen.var, col=i)
	}
	legend("topright", lty=1, col=rev(seq_along(recrates)), legend=rev(paste0("r = ", recrates)))
	mtext("Recombinations maintains the genetic variance (optimum = 0)", 3)
	

dev.off()

